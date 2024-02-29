/*
 * Copyright (c) 2021 Anthony J. Greenberg
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Summarize variant tables by hashing
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2021 Anthony J. Greenberg
 * \version 0.5
 *
 * Implementation of classes that take binary variant files and generate lossy summaries with hashing.
 *
 */

#include <chrono>
#include <ctime>
#include <iomanip>  // for put_time
#include <cstring>
#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include <ios>
#include <vector>
#include <unordered_map>
#include <array>
#include <utility>  // for std::pair
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <limits>
#include <future>
#include <thread>
#include <immintrin.h>

#include "gvarHash.hpp"
#include "vashFunctions.hpp"
#include "random.hpp"
#include "similarityMatrix.hpp"

using namespace BayesicSpace;

// GenoTableBin methods
constexpr size_t   GenoTableBin::nMagicBytes_    = 3;                // number of leading bytes for .bed files
constexpr uint8_t  GenoTableBin::oneBit_         = 0b00000001;       // One set bit for masking
constexpr uint8_t  GenoTableBin::byteSize_       = 8;                // Size of one byte in bits
constexpr uint8_t  GenoTableBin::bedGenoPerByte_ = 4;                // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableBin::llWordSize_     = 8;                // 64 bit word size in bytes

// Constructors
GenoTableBin::GenoTableBin(const std::string &inputFileName, const uint32_t &nIndividuals, std::string logFileName, const size_t &nThreads)
															: nIndividuals_{nIndividuals}, nThreads_{nThreads}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime = std::time(nullptr);
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype binarization from a .bed file started on " + logStream.str() + "\n";
	logStream.clear();
	if (nIndividuals <= 1) {
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting\n";
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nIndividuals > std::numeric_limits<size_t>::max() / nIndividuals ) { // a square will overflow
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too big to make a square relationship matrix; aborting\n";
		throw std::string("ERROR: the number of individuals (") + std::to_string(nIndividuals) + 
			std::string( ") is too big to make a square relationship matrix in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nThreads_ = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_ = std::max(nThreads_, 1UL);
	const uint32_t nBedBytesPerLocus = nIndividuals_ / bedGenoPerByte_ + static_cast<uint32_t>( (nIndividuals_ % bedGenoPerByte_) > 0);
	std::fstream inStr;
	// Start by measuring file size
	inStr.open(inputFileName, std::ios::in | std::ios::binary | std::ios::ate);
	if ( inStr.fail() ) {
		logMessages_ += "ERROR: failed to open file " + inputFileName + "\n";
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const uint64_t endPosition = static_cast<uint64_t>( inStr.tellg() );
	if (endPosition <= nMagicBytes_) {
		logMessages_ += "ERROR: no genotype records in file " + inputFileName + "\n";
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t nTotalBedBytes{static_cast<uint64_t>(endPosition) - nMagicBytes_};
	inStr.close();
	if ( nTotalBedBytes > std::numeric_limits<uint32_t>::max() ) {
		logMessages_ += "ERROR: .bed file (" + inputFileName + ") too large\n";
		throw std::string("ERROR: there must be fewer than 2^32 bytes in the .bed file ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nLoci_ = static_cast<uint32_t>(nTotalBedBytes) / nBedBytesPerLocus;
	logMessages_ += "Number of individuals: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: "        + std::to_string(nLoci_) + "\n";
	logMessages_ += "Number of threads: "     + std::to_string(nThreads_) + "\n";

	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStr.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	binLocusSize_ = nIndividuals_ / byteSize_ + static_cast<size_t>( (nIndividuals_ % byteSize_) > 0 );
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);
	BedDataStats locusGroupAttributes{};
	const size_t ramSize                = getAvailableRAM() / 2UL;                                                 // measuring here, after all the major allocations; use half to leave resources for other operations
	locusGroupAttributes.nLociToRead    = std::min( ramSize / nBedBytesPerLocus, static_cast<size_t>(nLoci_) );    // number of .bed loci to read at a time
	locusGroupAttributes.nMemChunks     = nLoci_ / locusGroupAttributes.nLociToRead;
	const size_t remainingLoci          = nLoci_ % locusGroupAttributes.nLociToRead;
	const size_t remainingBytes         = remainingLoci * nBedBytesPerLocus;
	locusGroupAttributes.nBytesToRead   = std::min( locusGroupAttributes.nLociToRead * nBedBytesPerLocus,
													static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) );
	locusGroupAttributes.nLociPerThread = std::max(locusGroupAttributes.nLociToRead / nThreads_, 1UL);
	locusGroupAttributes.nBytesPerLocus = nIndividuals_ / bedGenoPerByte_ + static_cast<size_t>(nIndividuals_ % bedGenoPerByte_ > 0);
	logMessages_                       += "RAM available for reading the .bed file: " + std::to_string(ramSize) + " bytes\n";
	logMessages_                       += ".bed file will be read in " + std::to_string(locusGroupAttributes.nMemChunks) + " chunk(s)\n";
	assert( ( remainingBytes < std::numeric_limits<std::streamsize>::max() ) //NOLINT
			&& "ERROR: remainingBytes larger than maximum streamsize in GenoTableBin constructor");

	locusGroupAttributes.firstLocusIdx = 0;
	locusGroupAttributes.firstLocusIdx = bed2bin_(locusGroupAttributes, inStr);
	if (remainingLoci > 0) {
		assert( ( remainingBytes < std::numeric_limits<std::streamsize>::max() ) // NOLINT
				&& "ERROR: remainingBytes exceeds maximum streamsize in GenoTableBin constructor" );
		locusGroupAttributes.nLociPerThread = std::max(remainingLoci / nThreads_, 1UL);
		locusGroupAttributes.nBytesToRead   = remainingBytes;
		locusGroupAttributes.nLociToRead    = remainingLoci;
		locusGroupAttributes.nMemChunks     = 1;
		bed2bin_(locusGroupAttributes, inStr);
	}
	inStr.close();
}

GenoTableBin::GenoTableBin(const std::vector<int> &maCounts, const uint32_t &nIndividuals, std::string logFileName, const size_t &nThreads)
							: nIndividuals_{nIndividuals}, nLoci_{static_cast<uint32_t>( maCounts.size() / static_cast<size_t>(nIndividuals) )}, nThreads_{nThreads}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime = std::time(nullptr);
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype binarization from minor allele count vector started on " + logStream.str() + "\n";
	logStream.clear();
	if ( ( maCounts.size() / static_cast<size_t>(nIndividuals) ) > std::numeric_limits<uint32_t>::max() ) {
		logMessages_ += "ERROR: too many loci \n";
		throw std::string("ERROR: there must be fewer than 2^32 loci in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nIndividuals <= 1) {
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting\n";
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( (maCounts.size() % nIndividuals) > 0 ) {
		logMessages_ += "ERROR: length of allele count vector (" + std::to_string( maCounts.size() ) + " is not divisible by the provided number of individuals (" +
			std::to_string(nIndividuals) + "\n";
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(nIndividuals) + std::string( ") in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( maCounts.empty() ) {
		logMessages_ += "ERROR: empty vector of minor allele counts\n";
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nThreads_     = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_     = std::max(nThreads_, 1UL);
	logMessages_ += "Number of individuals: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: " + std::to_string(nLoci_) + "\n";
	logMessages_ += "Number of threads: " + std::to_string(nThreads_) + "\n";
	binLocusSize_ = nIndividuals_ / byteSize_ + static_cast<size_t>( (nIndividuals_ % byteSize_) > 0 );
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);
	const size_t ranVecSize     = nIndividuals_ / llWordSize_ + static_cast<size_t>( (nIndividuals_ % llWordSize_) > 0 );
	const size_t nLociPerThread = nLoci_ / nThreads_;
	if (nLociPerThread > 0) {
		CountAndSize threadCounts{0, 0};
		threadCounts.count = nThreads_;
		threadCounts.size  = nLociPerThread;
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
		threadRanges.back().second = nLoci_;
		std::vector< std::future<void> > tasks;
		tasks.reserve(nThreads_);
		for (const auto &eachTR : threadRanges) {
			tasks.emplace_back(
				std::async([this, &maCounts, &eachTR, ranVecSize]{
					mac2binBlk_(maCounts, eachTR, ranVecSize);
				})
			);
		}
	} else {
		const std::pair<size_t, size_t> fullRange{0, nLoci_};
		mac2binBlk_(maCounts, fullRange, ranVecSize);
	}
}

void GenoTableBin::saveGenoBinary(const std::string &outFileName) const {
	std::fstream out;
	assert( ( binGenotypes_.size() < std::numeric_limits<std::streamsize>::max() ) // NOLINT
			&& "ERROR: binGenotypes_ size exceeds maximum streamsize in GenoTableBin.saveGenoBinary()");
	out.open(outFileName, std::ios::out | std::ios::binary | std::ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), static_cast<std::streamsize>( binGenotypes_.size() ) ); // OK because we are casting to const char*
	out.close();
}

void GenoTableBin::allJaccardLD(const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks) const {
	logMessages_ += "Calculating all pairwise LD using full Jaccard similarity estimates\n";
	std::vector<std::string> locusNames{};
	if ( !bimAndLDnames.inputFileName.empty() ) {
		std::fstream bimExistenceTest(bimAndLDnames.inputFileName, std::ios::in);
		const bool bimExists = bimExistenceTest.good();
		bimExistenceTest.close();
		if (bimExists) {
			logMessages_ += "Getting locus names from the " + bimAndLDnames.inputFileName + " .bim file\n";
			locusNames    = getLocusNames(bimAndLDnames.inputFileName);
		}
		assert( (locusNames.size() == nLoci_) // NOLINT
				&& "ERROR: number of loci in the .bim file not the same as nLoci_");
	}

	const SimilarityMatrix emptyMatrix;
	const size_t maxInRAM = getAvailableRAM() / ( static_cast<size_t>(2) * emptyMatrix.elementSize() );      // use half to leave resources for other operations
	const size_t nPairs   = nLoci_ * ( nLoci_ - static_cast<size_t>(1) ) / static_cast<size_t>(2);
	const size_t nChunks  = std::max(nPairs / maxInRAM, suggestNchunks);
	std::vector<size_t> chunkSizes{makeChunkSizes(nPairs, nChunks)};

	logMessages_ += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "; ";
	logMessages_ += "calculating in " + std::to_string(nChunks) + " chunk(s)\n";

	// set up the header
	std::fstream output;
	output.open(bimAndLDnames.outputFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\n";
	output.close();
	size_t cumChunkIdx{0};
	for (const auto &eachChunkSize : chunkSizes) {
		const size_t nThreads{std::min(eachChunkSize, nThreads_)};
		LocationWithLength currStartAndSize{};
		currStartAndSize.start  = cumChunkIdx;
		currStartAndSize.length = eachChunkSize;
		std::vector< std::pair<RowColIdx, RowColIdx> > threadRanges{makeChunkRanges(currStartAndSize, nThreads)};

		SimilarityMatrix result{jaccardThreaded_(threadRanges)};
		result.save(bimAndLDnames.outputFileName, nThreads);
		cumChunkIdx += eachChunkSize;
	}
}

void GenoTableBin::saveLogFile() const {
	std::fstream outLog;
	outLog.open(logFileName_, std::ios::out | std::ios::trunc);
	outLog << logMessages_;
	outLog.close();
}

void GenoTableBin::bed2binBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const LocationWithLength &locusSpan) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	RanDraw locPRNG;
	size_t begByte{locusSpan.start * binLocusSize_};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus) {
		LocationWithLength bedWindow{0, 0};
		bedWindow.start  = iBedLocus * locusSpan.length;
		bedWindow.length = locusSpan.length;
		LocationWithLength binWindow{0, 0};
		binWindow.start  = begByte;
		binWindow.length = binLocusSize_;
		binarizeBedLocus(bedWindow, bedData, nIndividuals_, binWindow, binGenotypes_);
		begByte += binLocusSize_;
	}
}

size_t GenoTableBin::bed2binThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const LocationWithLength &locusSpan) {
	size_t locusInd = locusSpan.start;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges) {
		tasks.emplace_back(
			std::async([this, &bedData, &eachTR, locusInd, &locusSpan]{
				LocationWithLength currentLocusSpan{0, 0};
				currentLocusSpan.start  = locusInd;
				currentLocusSpan.length = locusSpan.length;
				bed2binBlk_(bedData, eachTR, currentLocusSpan);
			})
		);
		locusInd += eachTR.second - eachTR.first;
	}
	return locusInd;
}

size_t GenoTableBin::bed2bin_(const BedDataStats &locusGroupStats, std::fstream &bedStream) {
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = locusGroupStats.nLociPerThread;
	size_t locusInd    = locusGroupStats.firstLocusIdx;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	assert( (locusGroupStats.nLociToRead >= threadRanges.back().second) // NOLINT
								&& "ERROR: nLociToRead smaller than threadRanges.back().second in bed2bin_" );
	const size_t excessLoci    = locusGroupStats.nLociToRead - threadRanges.back().second;
	threadRanges.back().second = locusGroupStats.nLociToRead;
	std::vector<char> bedChunkToRead(locusGroupStats.nBytesToRead, 0);
	for (size_t iChunk = 0; iChunk < locusGroupStats.nMemChunks; ++iChunk) {
		assert( ( locusGroupStats.nBytesToRead < std::numeric_limits<std::streamsize>::max() ) // NOLINT
								&& "ERROR: nBedBytesToRead exceeds maximum streamsize in bed2bin_" );
		bedStream.read( bedChunkToRead.data(), static_cast<std::streamsize>(locusGroupStats.nBytesToRead) );
		LocationWithLength currentLocusSpan{0, 0};
		currentLocusSpan.start  = locusInd;
		currentLocusSpan.length = locusGroupStats.nBytesPerLocus;
		locusInd                = bed2binThreaded_(bedChunkToRead, threadRanges, currentLocusSpan);
		locusInd               += excessLoci;
	}
	return locusInd;
}

void GenoTableBin::mac2binBlk_(const std::vector<int> &macData, const std::pair<size_t, size_t> &locusIndRange, const size_t &randVecLen) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	constexpr uint8_t middleMask{0b10000011};
	constexpr uint8_t endTwoBitMask{0b00000011};
	RanDraw locPRNG;
	auto remainderInd          = static_cast<uint8_t>(binLocusSize_ * byteSize_ - nIndividuals_);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
	remainderInd               = byteSize_ - remainderInd;
	// Create a vector to store random bytes for stochastic heterozygote resolution
	std::vector<uint64_t> rand(randVecLen);
	auto *randBytes    = reinterpret_cast<uint8_t*>( rand.data() );
	for (size_t iLocus = locusIndRange.first; iLocus < locusIndRange.second; ++iLocus) {
		// Fill the random byte vector
		for (auto &randValue : rand) {
			randValue = locPRNG.ranInt();
		}
		size_t iIndiv         = 0;
		size_t iRB            = 0;                                                                 // randBytes index
		const size_t begIndiv = iLocus * nIndividuals_;
		const size_t begByte  = iLocus * binLocusSize_;
		std::vector<uint8_t> missMasks(binLocusSize_, 0);
		size_t i0Byte = 0;                                                                         // to index the missMasks vector
		for (size_t iByte = begByte; iByte < begByte + binLocusSize_ - 1; ++iByte) {                // treat the last byte separately
			uint8_t binByte = 0;
			for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte) {
				auto curIndiv            = static_cast<uint8_t>(macData[begIndiv + iIndiv]);       // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= middleMask;                                            // mask everything in the middle
				const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
				missMasks[i0Byte]        |= static_cast<uint8_t>(missingMask << iInByte);
				curIndiv                 &= endTwoBitMask;
				const uint8_t randMask    = (randBytes[iRB] >> iInByte) & oneBit_;                 // 0b00000000 or 0b00000001 with equal chance
				uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
				curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
				binByte                  |= static_cast<uint8_t>(curBitMask << iInByte);
				++iIndiv;
			}
			// should be safe: each thread accesses different vector elements
			binGenotypes_[iByte] = binByte;
			++i0Byte;
			++iRB;
		}
		// now deal with the last byte in the individual
		uint8_t lastBinByte = 0;
		for (uint8_t iRem = 0; iRem < remainderInd; ++iRem) {
			auto curIndiv             = static_cast<uint8_t>(macData[begIndiv + iIndiv]);          // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= middleMask;                                                // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
			missMasks.back()         |= static_cast<uint8_t>(missingMask << iRem);
			curIndiv                 &= endTwoBitMask;                                                
			const uint8_t randMask    = (randBytes[binLocusSize_ - 1] >> iRem) & oneBit_;          // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                   // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                              // zero it out if missing value is set

			lastBinByte              |= static_cast<uint8_t>(curBitMask << iRem);
			++iIndiv;
		}
		// should be safe: each thread accesses different vector elements
		binGenotypes_[begByte + binLocusSize_ - 1] = lastBinByte;
		const LocationWithLength genoWindow{begByte, binLocusSize_};
		float maf = static_cast<float>( countSetBits(binGenotypes_, genoWindow) ) / static_cast<float>(nIndividuals_);
		if (maf > 0.5F) { // always want the alternative to be the minor allele
			i0Byte = 0;
			for (size_t i = begByte; i < begByte + binLocusSize_; ++i) {
				// should be safe: each thread accesses different vector elements
				binGenotypes_[i] = (~binGenotypes_[i]) & (~missMasks[i0Byte]);
				++i0Byte;
			}
			// should be safe: each thread accesses different vector elements
			binGenotypes_[begByte + binLocusSize_ - 1] &= lastByteMask; // unset the remainder bits
		}
	}
}

SimilarityMatrix GenoTableBin::jaccardBlock_(const std::pair<RowColIdx, RowColIdx> &blockRange) const {
	SimilarityMatrix result;
	uint32_t iRow{blockRange.first.iRow};
	if (blockRange.first.iRow != blockRange.second.iRow) {
		for (uint32_t jCol = blockRange.first.jCol; jCol < iRow; ++jCol) { // first, possibly incomplete, row
			RowColIdx localRC{};
			localRC.iRow = iRow;
			localRC.jCol = jCol;
			JaccardPair localJP{makeJaccardPair_(localRC)};
			result.insert(localRC, localJP);
		}
		++iRow;
		while ( iRow < (blockRange.second.iRow) ) { // complete triangle
			for (uint32_t jCol = 0; jCol < iRow; ++jCol) {
				RowColIdx localRC{};
				localRC.iRow = iRow;
				localRC.jCol = jCol;
				JaccardPair localJP{makeJaccardPair_(localRC)};
				result.insert(localRC, localJP);
			}
			++iRow;
		}
	}
	for (uint32_t jColRem = 0; jColRem < blockRange.second.jCol; ++jColRem) { // last, possibly incomplete, row
		RowColIdx localRC{};
		localRC.iRow = iRow;
		localRC.jCol = jColRem;
		JaccardPair localJP{makeJaccardPair_(localRC)};
		result.insert(localRC, localJP);
	}
	return result;
}

SimilarityMatrix GenoTableBin::jaccardThreaded_(const std::vector< std::pair<RowColIdx, RowColIdx> > &indexPairs) const {
	std::vector<SimilarityMatrix> threadResults( indexPairs.size() );
	std::vector< std::future<void> > tasks;
	tasks.reserve( indexPairs.size() );
	size_t iThread{0};
	for (const auto &eachRange : indexPairs) {
		tasks.emplace_back(
			std::async([this, eachRange, iThread, &threadResults]{
				threadResults.at(iThread) = jaccardBlock_(eachRange);
			})
		);
		++iThread;
	}

	for (const auto &eachThread : tasks) {
		eachThread.wait();
	}

	std::for_each(
		threadResults.begin() + 1,
		threadResults.end(),
		[&threadResults](SimilarityMatrix &eachMatrix) {
			threadResults.at(0).merge( std::move(eachMatrix) );
		}
	);

	return threadResults.at(0);
}

JaccardPair GenoTableBin::makeJaccardPair_(const RowColIdx &rowColumn) const {
	std::vector<uint8_t> locus(binLocusSize_);
	JaccardPair localJP{};
	const size_t rowBin = rowColumn.iRow * binLocusSize_;
	const size_t colBin = rowColumn.jCol * binLocusSize_;
	for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc) {
		locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] & binGenotypes_[colBin + iBinLoc];
	}
	localJP.nIntersect = countSetBits(locus);
	for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc) {
		locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] | binGenotypes_[colBin + iBinLoc];
	}
	localJP.nUnion = countSetBits(locus);
	return localJP;
}

// GenoTableHash methods
constexpr size_t   GenoTableHash::nMagicBytes_    = 3;                                    // number of leading bytes for .bed files
constexpr uint8_t  GenoTableHash::oneBit_         = 0b00000001;                           // One set bit for masking 
constexpr uint8_t  GenoTableHash::byteSize_       = 8;                                    // Size of one byte in bits 
constexpr uint8_t  GenoTableHash::bedGenoPerByte_ = 4;                                    // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableHash::llWordSize_     = 8;                                    // 64 bit word size in bytes 
constexpr uint32_t GenoTableHash::roundMask_      = 0xfffffff8;                           // mask for rounding down to nearest whole-byte value
constexpr uint64_t GenoTableHash::allBitsSet_     = std::numeric_limits<uint64_t>::max(); // 64-bit word with all bits set
constexpr size_t   GenoTableHash::wordSizeInBits_ = 64;                                   // 64-bit word size
constexpr uint16_t GenoTableHash::emptyBinToken_  = std::numeric_limits<uint16_t>::max(); // Value corresponding to an empty token 

// Constructors
GenoTableHash::GenoTableHash(const std::string &inputFileName, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, std::string logFileName)
								: kSketches_{indivSketchCounts.kSketches}, fSketches_{static_cast<float>(indivSketchCounts.kSketches)}, nLoci_{0}, nThreads_{nThreads}, emptyBinIdxSeed_{0}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime{std::time(nullptr)};
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype hashing from a .bed file started on " + logStream.str() + "\n";
	logStream.clear();
	if (indivSketchCounts.nIndividuals <= 1) {
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(indivSketchCounts.nIndividuals) + ") is too small; aborting\n";
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ < 3) {
		logMessages_ += "ERROR: number of sketches (" + std::to_string(kSketches_) + ") is too small; aborting\n";
		throw std::string("ERROR: sketch number must be at least three in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ > indivSketchCounts.nIndividuals) {
		logMessages_ += "ERROR: number of sketches (" + std::to_string(kSketches_) + ") is larger than the number of individuals; aborting\n";
		throw std::string("ERROR: sketch number must be smaller than the number of individuals in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	// Round up the number of individuals to nearest divisible by kSketches_
	sketchSize_   = indivSketchCounts.nIndividuals / kSketches_ + static_cast<uint32_t>( (indivSketchCounts.nIndividuals % kSketches_) > 0 );
	nIndividuals_ = sketchSize_ * kSketches_;
	if (indivSketchCounts.kSketches >= emptyBinToken_) {
		logMessages_ += "ERROR: sketch size (" + std::to_string(indivSketchCounts.kSketches) + ") is too big; aborting\n";
		throw std::string("ERROR: Number of sketches (") + std::to_string(indivSketchCounts.kSketches) + std::string(") implies sketch size (") +
			std::to_string(sketchSize_) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t nBedBytes{indivSketchCounts.nIndividuals / bedGenoPerByte_ + static_cast<size_t>( (indivSketchCounts.nIndividuals % bedGenoPerByte_) > 0 )};
	nThreads_     = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_     = std::max(nThreads_, 1UL);
	logMessages_ += "Number of threads used: " + std::to_string(nThreads_) + "\n";
	std::fstream inStream;
	// Start by measuring file size
	inStream.open(inputFileName, std::ios::in | std::ios::binary | std::ios::ate);
	if ( inStream.fail() ) {
		logMessages_ += "ERROR: failed to open file " + inputFileName + "; aborting\n";
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const auto endPosition{static_cast<size_t>( inStream.tellg() )};
	if (endPosition <= nMagicBytes_) {
		logMessages_ += "ERROR: no loci in the input .bed file " + inputFileName + "; aborting\n";
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	inStream.close();
	const size_t fileSize{endPosition - nMagicBytes_};
	const size_t tmpNloci{fileSize / nBedBytes};
	if ( tmpNloci > std::numeric_limits<uint32_t>::max() ) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(tmpNloci) + "\n";
		throw std::string("ERROR: there must be fewer than 2^32 loci in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nLoci_        = static_cast<uint32_t>(tmpNloci);
	logMessages_ += "Number of individuals: "         + std::to_string(indivSketchCounts.nIndividuals) + "\n";
	logMessages_ += "Number of individuals to hash: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: "                + std::to_string(nLoci_) + "\n";
	logMessages_ += "Hash size: "                     + std::to_string(kSketches_) + "\n";

	RanDraw prng;
	emptyBinIdxSeed_ = prng.ranInt();
	locusSize_       = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;                    // round up to the nearest multiple of 8
	nFullWordBytes_  = (nIndividuals_ - 1) / byteSize_;
	sketches_.resize(static_cast<size_t>(kSketches_) * nLoci_, emptyBinToken_);
	inStream.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStream.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	BedDataStats locusGroupAttributes{};
	locusGroupAttributes.nBytesPerLocus = indivSketchCounts.nIndividuals / bedGenoPerByte_ + static_cast<size_t>(indivSketchCounts.nIndividuals % bedGenoPerByte_ > 0);
	const size_t ramSize                = getAvailableRAM() / 2UL;                                    // measuring here, after all the major allocations; use half to leave resources for other operations
	locusGroupAttributes.nLociToRead    = std::min( ramSize / locusGroupAttributes.nBytesPerLocus, static_cast<size_t>(nLoci_) );       // number of .bed loci to read at a time
	const size_t remainingLoci          = nLoci_ % locusGroupAttributes.nLociToRead;
	const size_t remainingBytes         = remainingLoci * locusGroupAttributes.nBytesPerLocus;
	locusGroupAttributes.nMemChunks     = nLoci_ / locusGroupAttributes.nLociToRead;
	locusGroupAttributes.nBytesToRead   = std::min( locusGroupAttributes.nLociToRead * locusGroupAttributes.nBytesPerLocus,
													static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) );
	locusGroupAttributes.nLociPerThread = locusGroupAttributes.nLociToRead / nThreads_;

	logMessages_ += "RAM available for reading the .bed file: " + std::to_string(ramSize) + " bytes\n";
	logMessages_ += ".bed file will be read in "                + std::to_string(locusGroupAttributes.nMemChunks) + " chunk(s)\n";

	// Sample with replacement additional individuals to pad out the total
	std::vector< std::pair<size_t, size_t> > addIndv;
	for (size_t iAddIndiv = indivSketchCounts.nIndividuals; iAddIndiv < nIndividuals_; ++iAddIndiv) {
		addIndv.emplace_back( iAddIndiv, prng.sampleInt(indivSketchCounts.nIndividuals) );
	}
	if ( !addIndv.empty() ) {
		std::string addIndexes;
		for (const auto &eachIdx : addIndv) {
			addIndexes += std::to_string(eachIdx.second) + " ";
		}
		logMessages_ += "Re-sampled individuals: " + addIndexes + "\n";
	}
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{prng.fyIndexesUp(nIndividuals_)};

	locusGroupAttributes.firstLocusIdx = 0;
	locusGroupAttributes.firstLocusIdx = bed2oph_(locusGroupAttributes, inStream, ranInts, addIndv);
	if (remainingLoci > 0) {
		locusGroupAttributes.nLociPerThread = std::max(remainingLoci / nThreads_, 1UL);
		locusGroupAttributes.nBytesToRead   = remainingBytes;
		locusGroupAttributes.nLociToRead    = remainingLoci;
		locusGroupAttributes.nMemChunks     = 1;
		bed2oph_(locusGroupAttributes, inStream, ranInts, addIndv);
	}
	inStream.close();
}

GenoTableHash::GenoTableHash(const std::vector<int> &maCounts, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, std::string logFileName) 
		: nIndividuals_{indivSketchCounts.nIndividuals}, kSketches_{indivSketchCounts.kSketches}, fSketches_{static_cast<float>(indivSketchCounts.kSketches)},
				nLoci_{static_cast<uint32_t>(maCounts.size() / indivSketchCounts.nIndividuals)}, nThreads_{nThreads}, emptyBinIdxSeed_{0}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime = std::time(nullptr);
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype hashing from a minor allele count vector started on " + logStream.str() + "\n";
	logStream.clear();
	if (indivSketchCounts.nIndividuals <= 1) {
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(indivSketchCounts.nIndividuals) + ") is too small; aborting\n";
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( (maCounts.size() % indivSketchCounts.nIndividuals) > 0) {
		logMessages_ += "ERROR: minor allele vector size (" + std::to_string( maCounts.size() ) + ") is not evenly divisible by the number of individuals (" +
							std::to_string(nIndividuals_) + "); aborting\n";
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(indivSketchCounts.nIndividuals) + std::string( ") in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( maCounts.empty() ) {
		logMessages_ += "ERROR: minor allele count vector is empty; aborting\n";
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ < 3) {
		logMessages_ += "ERROR: sketch size (" + std::to_string(kSketches_) + ") is too small; aborting\n";
		throw std::string("ERROR: sketch size must be at least three in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ > indivSketchCounts.nIndividuals) {
		logMessages_ += "ERROR: number of sketches (" + std::to_string(kSketches_) + ") is larger than the number of individuals; aborting\n";
		throw std::string("ERROR: sketch number must be smaller than the number of individuals in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	RanDraw prng;
	emptyBinIdxSeed_ = prng.ranInt();

	nThreads_     = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_     = std::max(nThreads_, 1UL);
	sketchSize_   = indivSketchCounts.nIndividuals / kSketches_ + static_cast<uint16_t>( (indivSketchCounts.nIndividuals % kSketches_) > 0 );
	nIndividuals_ = sketchSize_ * kSketches_;
	if (indivSketchCounts.kSketches >= emptyBinToken_) {
		logMessages_ += "ERROR: sketch size (" + std::to_string(indivSketchCounts.kSketches) + ") is too small; aborting\n";
		throw std::string("ERROR: Number of sketches (") + std::to_string(kSketches_) + std::string(") implies sketch size (") +
			std::to_string(indivSketchCounts.kSketches) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	locusSize_      = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;                    // round up to the nearest multiple of 8
	nFullWordBytes_ = (nIndividuals_ - 1) / byteSize_;
	sketches_.resize(static_cast<size_t>(kSketches_) * nLoci_, emptyBinToken_);
	const size_t ranVecSize = locusSize_ / llWordSize_ + static_cast<size_t>( (locusSize_ % llWordSize_) > 0);
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{prng.fyIndexesUp(nIndividuals_)};

	logMessages_ += "Number of threads used: " + std::to_string(nThreads_) + "\n";
	logMessages_ += "Number of individuals: "  + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: "         + std::to_string(nLoci_) + "\n";
	logMessages_ += "Hash size: "              + std::to_string(kSketches_) + "\n";

	const size_t nLociPerThread = nLoci_ / nThreads_;
	if (nLociPerThread > 0) {
		CountAndSize threadCounts{0, 0};
		threadCounts.count = nThreads_;
		threadCounts.size  = nLociPerThread;
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
		threadRanges.back().second = nLoci_;
		std::vector< std::future<void> > tasks;
		tasks.reserve(nThreads_);
		for (const auto &eachTR : threadRanges) {
			tasks.emplace_back(
				std::async([this, &maCounts, &eachTR, ranVecSize, &ranInts]{
					LocationWithLength threadLoci{0, 0};
					threadLoci.start  = eachTR.first;
					threadLoci.length = eachTR.second - eachTR.first;
					mac2ophBlk_(maCounts, threadLoci, ranVecSize, ranInts);
				})
			);
		}
		for (const auto &eachTask : tasks) {
			eachTask.wait();
		}
	} else {
		LocationWithLength allLoci{0, 0};
		allLoci.length = nLoci_;
		mac2ophBlk_(maCounts, allLoci, ranVecSize, ranInts);
	}
}

void GenoTableHash::allHashLD(const float &similarityCutOff, const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks) const {
	std::vector<uint32_t> allLocusIndexes(nLoci_);
	std::iota(allLocusIndexes.begin(), allLocusIndexes.end(), 0);

	const SimilarityMatrix emptyMatrix;
	const size_t maxInRAM = getAvailableRAM() / ( static_cast<size_t>(2) * emptyMatrix.elementSize() );      // use half to leave resources for other operations
	const size_t nPairs   = static_cast<size_t>(nLoci_) * (static_cast<size_t>(nLoci_) - 1UL) / 2UL;
	const size_t nChunks  = std::max(nPairs / maxInRAM, suggestNchunks);
	std::vector<size_t> chunkSizes{makeChunkSizes(nPairs, nChunks)};

	logMessages_ += "Calculating all pairwise LD\n";
	logMessages_ += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "\n";
	logMessages_ += "calculating in " + std::to_string(nChunks) + " chunk(s)\n";

	std::fstream output;
	output.open(bimAndLDnames.outputFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\n";
	output.close();
	size_t cumChunkIdx{0};
	for (const auto &eachChunkSize : chunkSizes) {
		const size_t nThreads{std::min(eachChunkSize, nThreads_)};
		LocationWithLength currStartAndSize{};
		currStartAndSize.start  = cumChunkIdx;
		currStartAndSize.length = eachChunkSize;

		std::vector< std::pair<RowColIdx, RowColIdx> > threadRanges{makeChunkRanges(currStartAndSize, nThreads)};
		SimilarityMatrix result{hashJacThreaded_(threadRanges, allLocusIndexes, similarityCutOff)};
		result.save(bimAndLDnames.outputFileName, nThreads, bimAndLDnames.inputFileName);

		cumChunkIdx += eachChunkSize;
	}
}

std::vector<HashGroup> GenoTableHash::makeLDgroups(const size_t &nRowsPerBand) const {
	assert( (nRowsPerBand != 0) // NOLINT
			&& "ERROR: nRowsPerBand must not be 0 in makeLDgroups()" );
	assert( (nRowsPerBand < kSketches_) // NOLINT
			&& "ERROR: nRowsPerBand must be less than kSketches_ in makeLDgroups()" );
	const size_t nBands = kSketches_ / nRowsPerBand;                                                          // only using full-size bands because smaller ones permit inclusion of low-similarity pairs
	assert( ( nBands <= std::numeric_limits<uint16_t>::max() ) // NOLINT
			&& "ERROR: number of bands cannot exceed uint16_t max in makeLDgroups()" );

	logMessages_ += "Grouping loci\n";
	logMessages_ += "Number of rows per band: " + std::to_string(nRowsPerBand) + "\n";
	logMessages_ += "Number of bands: "         + std::to_string(nBands) + "\n";

	RanDraw prng;
	const auto sketchSeed = static_cast<uint32_t>( prng.ranInt() );
	std::unordered_map< uint32_t, std::vector<uint32_t> > ldGroups;                                           // the hash table

	for (size_t iLocus = 0; iLocus < nLoci_; ++iLocus) {
		size_t iSketch = 0;
		for (uint16_t iBand = 0; iBand < static_cast<uint16_t>(nBands); ++iBand) {
			std::vector<uint16_t> bandVec{iBand};                                                             // add the band index to the hash, so that only corresponding bands are compared

			auto firstSketchIt = sketches_.cbegin()
				+ static_cast<std::vector<uint16_t>::difference_type>(iSketch + iLocus * kSketches_);         // iSketch tracks band IDs
			auto lastSketchIt = firstSketchIt
				+ static_cast<std::vector<uint16_t>::difference_type>(nRowsPerBand);
			std::copy( firstSketchIt, lastSketchIt, std::back_inserter(bandVec) );

			LocationWithLength bandVecWindow{0, 0};
			bandVecWindow.start  = 0;
			bandVecWindow.length = bandVec.size();
			const uint32_t hash  = murMurHash(bandVec, bandVecWindow, sketchSeed);
			ldGroups[hash].push_back( static_cast<uint32_t>(iLocus) );
			iSketch += nRowsPerBand;
		}
	}
	std::vector< std::vector<uint32_t> > groups;
	for (auto &eachGrp : ldGroups) {
		// remove groups with one locus
		if (eachGrp.second.size() >= 2) {
			groups.emplace_back(eachGrp.second);
		}
	}
	// pre-sort the groups by position
	// this carries a ~30% overhead, but speeds the downstream pair sorting
	// enough that overall there is a ~10% speed-up and at least no overhead
	// it also enables processing by chunks if the whole sparse table does not fit in RAM
	std::sort( groups.begin(), groups.end(),
				[](const std::vector<uint32_t> &group1, const std::vector<uint32_t> &group2) {
					return (group1[0] == group2[0] ? group1[1] < group2[1] : group1[0] < group2[0]);
				}
			);
	// de-duplicate the groups
	logMessages_ += "Number of groups before de-duplication: " + std::to_string( groups.size() ) + "\n";
	auto lastUniqueIt = std::unique(
		groups.begin(),
		groups.end(),
		[sketchSeed](const std::vector<uint32_t> &first, const std::vector<uint32_t> &second) {
			return murMurHash(first, sketchSeed) == murMurHash(second, sketchSeed);
		}
	);
	groups.erase( lastUniqueIt, groups.end() );
	groups.shrink_to_fit();

	logMessages_ += "Number of groups after de-duplication: " + std::to_string( groups.size() ) + "\n";

	std::vector<HashGroup> indexedGroups;
	indexedGroups.reserve( groups.size() );
	uint64_t cumNpairs{0};
	for (auto &eachGroup : groups) {
		HashGroup currHG;
		cumNpairs              += eachGroup.size() * (eachGroup.size() - 1) / 2;
		currHG.cumulativeNpairs = cumNpairs;
		currHG.locusIndexes     = std::move(eachGroup);
		indexedGroups.emplace_back(currHG);
	}

	return indexedGroups;
}

void GenoTableHash::makeLDgroups(const size_t &nRowsPerBand, const InOutFileNames &bimAndGroupNames) const {
	std::vector<HashGroup> ldGroups{this->makeLDgroups(nRowsPerBand)};
	logMessages_ += "Saving group IDs only\n";
	std::vector<std::string> locusNames{};
	if ( !bimAndGroupNames.inputFileName.empty() ) {
		std::fstream bimExistenceTest(bimAndGroupNames.inputFileName, std::ios::in);
		const bool bimExists = bimExistenceTest.good();
		bimExistenceTest.close();
		if (bimExists) {
			logMessages_ += "Getting locus names from the " + bimAndGroupNames.inputFileName + " .bim file\n";
			locusNames    = getLocusNames(bimAndGroupNames.inputFileName);
		}
		assert( (locusNames.size() == nLoci_) // NOLINT
				&& "ERROR: number of loci in the .bim file not the same as nLoci_");
	}

	std::fstream out;
	out.open(bimAndGroupNames.outputFileName, std::ios::out | std::ios::trunc);
	out << "groupID\tlocusIdx\n";
	uint32_t groupID{1};
	if ( locusNames.empty() ) {
		for (const auto &eachGroup : ldGroups) {
			for (const auto &locusIdx : eachGroup.locusIndexes) {
				out << "G" << groupID << "\t" << locusIdx + 1 << "\n";
			}
			++groupID;
		}
		out.close();
		return;
	}
	for (const auto &eachGroup : ldGroups) {
		for (const auto &locusIdx : eachGroup.locusIndexes) {
			out << "G" << groupID << "\t" << locusNames[locusIdx] << "\n";
		}
		++groupID;
	}
	out.close();
}

void GenoTableHash::ldInGroups(const SparsityParameters &sparsityValues, const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks) const {
	std::vector<HashGroup> ldGroups{this->makeLDgroups(sparsityValues.nRowsPerBand)};
	
	const size_t totalPairNumber{ldGroups.back().cumulativeNpairs};                                                                                    // total number of pairs
	logMessages_ += "Estimating LD in groups\n";
	logMessages_ += "number of pairs in the hash table: " + std::to_string(totalPairNumber) + "\n";

	const SimilarityMatrix emptyMatrix;
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * emptyMatrix.elementSize() );      // use half to leave resources for other operations
	const size_t nChunks  = std::max(totalPairNumber / maxInRAM, suggestNchunks);

	logMessages_ += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "\n";
	logMessages_ += "calculating in " + std::to_string(nChunks) + " chunk(s)\n";

	std::fstream output;
	output.open(bimAndLDnames.outputFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\n";
	output.close();

	const std::vector<size_t> chunkSizes{makeChunkSizes( totalPairNumber, std::min(nChunks, totalPairNumber) )};
	BayesicSpace::HashGroupItPairCount startPair{};
	startPair.hgIterator = ldGroups.cbegin();
	startPair.pairCount  = 0;
	const size_t lastPairNumber{ldGroups.back().locusIndexes.size() * (ldGroups.back().locusIndexes.size() - 1) / 2};
	for (const auto &eachCS : chunkSizes) {
		SimilarityMatrix groupSimilarities;
		// actual matrix sizes may be smaller than expected because of common pairs among groups
		// so we keep adding until we run out of space to limit the number of saves and pair duplication
		// that can result from not being able to de-duplicate pairs that are already saved
		while (groupSimilarities.nElements() < eachCS) {
			const size_t currentChunkSize = eachCS - groupSimilarities.nElements();
			std::vector< std::pair<HashGroupItPairCount, HashGroupItPairCount> > groupRanges;
			const std::vector<size_t> threadSizes{makeChunkSizes( currentChunkSize, std::min(nThreads_, currentChunkSize) )};
			groupRanges.reserve( threadSizes.size() );
			for (const auto &eachThrSize : threadSizes) {
				groupRanges.emplace_back( makeGroupRanges(ldGroups, startPair, eachThrSize) );
				startPair = groupRanges.back().second;
			}
			groupSimilarities.merge( hashJacThreaded_(groupRanges, sparsityValues.similarityCutOff) );
			if ( ( startPair.hgIterator == std::prev( ldGroups.cend() ) ) && (startPair.pairCount == lastPairNumber) ) {
				groupSimilarities.save(bimAndLDnames.outputFileName, nThreads_, bimAndLDnames.inputFileName);
				return;
			}
		}
		groupSimilarities.save(bimAndLDnames.outputFileName, nThreads_, bimAndLDnames.inputFileName);
	}
	output.close();
}

void GenoTableHash::saveLogFile() const {
	std::fstream outLog;
	outLog.open(logFileName_, std::ios::out | std::ios::trunc);
	outLog << logMessages_;
	outLog.close();
}

void GenoTableHash::permuteBits_(const std::vector<size_t> &permutationIdx, std::vector<uint8_t> &binLocus) const {
	size_t iIndiv = 0;
	size_t iByte  = 0;
	while(iByte < nFullWordBytes_) {
		for (uint8_t iInLocusByte = 0; iInLocusByte < byteSize_; ++iInLocusByte) {
			auto bytePair            = static_cast<uint16_t>(binLocus[iByte]);
			const size_t perIndiv    = permutationIdx[iIndiv++];                                                           // post-increment to use current value for index first
			const size_t permByteInd = perIndiv / byteSize_;
			const auto permInByteInd = static_cast<uint8_t>( perIndiv - (perIndiv & roundMask_) );
			// Pair the current locus byte with the byte containing the value to be swapped
			// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
			bytePair                |= static_cast<uint16_t>(binLocus[permByteInd] << byteSize_);
			const auto mask          = static_cast<uint16_t>(oneBit_ << iInLocusByte);
			const auto perMask       = static_cast<uint8_t>(oneBit_ << permInByteInd);
			const auto shiftDistance = static_cast<uint16_t>( (byteSize_ - iInLocusByte) + permInByteInd );                // subtraction is safe b/c byteSize is the loop terminator
			const auto temp1         = static_cast<uint16_t>( ( bytePair ^ (bytePair >> shiftDistance) ) & mask );
			const auto temp2         = static_cast<uint16_t>(temp1 << shiftDistance);
			bytePair                ^= temp1 ^ temp2;
			// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus)
			// Must modify the current byte in each loop iteration because permutation indexes may fall into it
			binLocus[iByte]       ^= static_cast<uint8_t>( ( binLocus[iByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask) );
			binLocus[permByteInd] ^= static_cast<uint8_t>( ( binLocus[permByteInd] ^ static_cast<uint8_t>(bytePair >> byteSize_) ) & perMask );
 		}
		++iByte;
	}
	// Finish the individuals in the remaining partial byte, if any
	uint8_t iInLocusByte = 0;
	while (iIndiv < nIndividuals_ - 1) {
		auto bytePair            = static_cast<uint16_t>(binLocus[iByte]);
		const size_t perIndiv    = permutationIdx[iIndiv++];                                   // post-increment to use current value for index first
		const size_t permByteInd = perIndiv / byteSize_;
		const auto permInByteInd = static_cast<uint8_t>( perIndiv - (perIndiv & roundMask_) );
		// Pair the current locus byte with the byte containing the value to be swapped
		// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
		bytePair                |= static_cast<uint16_t>(binLocus[permByteInd] << byteSize_);
		const auto mask          = static_cast<uint16_t>(oneBit_ << iInLocusByte);
		const auto perMask       = static_cast<uint8_t>(oneBit_ << permInByteInd);
		const auto shiftDistance = static_cast<uint16_t>( (byteSize_ - iInLocusByte) + permInByteInd );           // subtraction is safe b/c byteSize is the loop terminator
		const auto temp1         = static_cast<uint16_t>( ( bytePair ^ (bytePair >> shiftDistance) ) & mask );
		const auto temp2         = static_cast<uint16_t>(temp1 << shiftDistance);
		bytePair                ^= temp1 ^ temp2;
		// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus)
		// Must modify the current byte in each loop iteration because permutation indexes may fall into it
		binLocus[iByte]       ^= static_cast<uint8_t>( ( binLocus[iByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask) );
		binLocus[permByteInd] ^= static_cast<uint8_t>( ( binLocus[permByteInd] ^ static_cast<uint8_t>(bytePair >> byteSize_) ) & perMask );
		++iInLocusByte;
	}
}

void GenoTableHash::locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint8_t> &binLocus) {
	// Start with a permutation to make OPH
	permuteBits_(permutation, binLocus);
	// Now make the sketches
	RanDraw prng(emptyBinIdxSeed_);
	std::vector<uint32_t> seeds{static_cast<uint32_t>( prng.ranInt() )};
	std::vector<size_t> filledIndexes;                                                       // indexes of the non-empty sketches
	size_t iByte{0};
	size_t sketchBeg{locusInd * kSketches_};
	size_t iSketch{0};
	uint64_t sketchTail{0};                                                                  // left over buts from beyond the last full byte of the previous sketch
	size_t locusChunkSize = (llWordSize_ > binLocus.size() ? binLocus.size() : llWordSize_);
	while ( iByte < binLocus.size() ) {
		uint64_t nWordUnsetBits{wordSizeInBits_};
		uint64_t nSumUnsetBits{0};
		while ( (nWordUnsetBits == wordSizeInBits_) && ( iByte < binLocus.size() ) ) {
			uint64_t locusChunk{allBitsSet_};
			assert( ( iByte < binLocus.size() ) // NOLINT
					&& "ERROR: iByte must be less than locus size in bytes in locusOPH_()" );
			const size_t nRemainingBytes{binLocus.size() - iByte};
			locusChunkSize = static_cast<size_t>(nRemainingBytes >= llWordSize_) * llWordSize_ + static_cast<size_t>(nRemainingBytes < llWordSize_) * nRemainingBytes;
			memcpy(&locusChunk, binLocus.data() + iByte, locusChunkSize);
			locusChunk    &= allBitsSet_ << sketchTail;
			nWordUnsetBits = _tzcnt_u64(locusChunk);
			nSumUnsetBits += nWordUnsetBits - sketchTail;
			sketchTail     = 0;
			iByte         += locusChunkSize;
		}
		iSketch += nSumUnsetBits / sketchSize_;
		if (iSketch >= kSketches_) {
			break;
		}
		filledIndexes.push_back(iSketch);
		sketches_[sketchBeg + iSketch] = static_cast<uint16_t>(nSumUnsetBits % sketchSize_);
		++iSketch;
		const uint64_t bitsDone{iSketch * sketchSize_};
		iByte      = bitsDone / byteSize_;
		sketchTail = bitsDone % byteSize_;
	}
	assert( (filledIndexes.size() <= kSketches_) // NOLINT
					&& "ERROR: filledIndexes.size() must not be greater than sketch number (kSketches_) in locusOPH_()" );
	size_t iSeed = 0;                                           // index into the seed vector
	if ( filledIndexes.empty() ) {                              // if the whole locus is monomorphic, pick a random index as filled
		filledIndexes.push_back( prng.sampleInt(kSketches_) );
	}
	size_t emptyCount = kSketches_ - filledIndexes.size();
	while (emptyCount > 0) {
		for (const auto eachFI : filledIndexes) {
			std::array<uint32_t, SIZE_OF_SIZET> key{};
			memcpy(key.data(), &eachFI, SIZE_OF_SIZET);
			auto newIdx = static_cast<uint32_t>(murMurHash(key, seeds[iSeed]) % kSketches_ + sketchBeg);
			// should be safe: each thread accesses different vector elements
			if (sketches_[newIdx] == emptyBinToken_) {
				sketches_[newIdx] = sketches_[eachFI + sketchBeg];
				--emptyCount;
				break;
			}
		}
		++iSeed;
		if ( iSeed == seeds.size() ) {
			seeds.push_back( static_cast<uint32_t>( prng.ranInt() ) );
		}
	}
}

void GenoTableHash::bed2ophBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t>&bedLocusIndRange,
					const LocationWithLength &bedLocusSpan, const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	size_t iLocus{bedLocusSpan.start};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus) {
		std::vector<uint8_t> binLocus(locusSize_, 0);
		LocationWithLength bedWindow{0, 0};
		bedWindow.start  = iBedLocus * bedLocusSpan.length;
		bedWindow.length = bedLocusSpan.length;
		LocationWithLength binWindow{0, 0};
		binWindow.length = locusSize_;
		binarizeBedLocus(bedWindow, bedData, nIndividuals_, binWindow, binLocus);
		// pad the locus to have a whole number of sketches 
		for (const auto &addI : padIndiv) {
			const size_t iLocByte    = addI.first / byteSize_;
			const auto iInLocByte    = static_cast<uint8_t>(addI.first % byteSize_);
			auto bytePair            = static_cast<uint16_t>(binLocus[iLocByte]);
			const size_t permByteInd = addI.second / byteSize_;
			const auto permInByteInd = static_cast<uint8_t>(addI.second % byteSize_);
			// Pair the current locus byte with the byte containing the value to be swapped
			// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
			bytePair                |= static_cast<uint16_t>(binLocus[permByteInd] << byteSize_);
			const auto mask          = static_cast<uint16_t>(1 << iInLocByte);
			const auto shiftDistance = static_cast<uint16_t>( (byteSize_ - iInLocByte) + permInByteInd );                        // subtraction is safe b/c iInLocByte is modulo byteSize
			const auto temp1         = static_cast<uint16_t>( ( bytePair ^ (bytePair >> shiftDistance) ) & mask );
			const auto temp2         = static_cast<uint16_t>(temp1 << shiftDistance);
			bytePair                ^= temp1 ^ temp2;
			// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus1)
			// Must modify the current byte in each loop iteration because permutation indexes may fall into it
			binLocus[iLocByte] ^= static_cast<uint8_t>( ( binLocus[iLocByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask) );
		}
		locusOPH_(iLocus, permutation, binLocus);
		++iLocus;
	}
}

size_t GenoTableHash::bed2ophThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges,
						const LocationWithLength &bedLocusSpan, const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv) {
	size_t locusInd = bedLocusSpan.start;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges) {
		LocationWithLength currentBedLocusSpan{0, 0};
		currentBedLocusSpan.start  = locusInd;
		currentBedLocusSpan.length = bedLocusSpan.length;
		tasks.emplace_back(
			std::async([this, &bedData, &eachTR, currentBedLocusSpan, &permutation, &padIndiv]{ // must pass currentLocusSpan by value, otherwise the threads see the same values
				bed2ophBlk_(bedData, eachTR, currentBedLocusSpan, permutation, padIndiv);
			})
		);
		locusInd += eachTR.second - eachTR.first;
	}
	for (const auto &eachTask : tasks) {
		eachTask.wait();
	}
	return locusInd;
}

size_t GenoTableHash::bed2oph_(const BedDataStats &locusGroupStats, std::fstream &bedStream, const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv) {
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = locusGroupStats.nLociPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	assert( (locusGroupStats.nLociToRead >= threadRanges.back().second) // NOLINT
				&& "ERROR: nLociToRead smaller than threadRanges.back().second in bed2oph_()");
	const size_t excessLoci = locusGroupStats.nLociToRead - threadRanges.back().second;
	threadRanges.back().second = locusGroupStats.nLociToRead;
	std::vector<char> bedChunkToRead(locusGroupStats.nBytesToRead, 0);
	size_t locusInd{locusGroupStats.firstLocusIdx};
	for (size_t iChunk = 0; iChunk < locusGroupStats.nMemChunks; ++iChunk) {
		assert( ( locusGroupStats.nBytesToRead < std::numeric_limits<std::streamsize>::max() ) // NOLINT
				&& "ERROR: amount to read larger than maximum streamsize in bed2oph_()");
		bedStream.read( bedChunkToRead.data(), static_cast<std::streamsize>(locusGroupStats.nBytesToRead) );
		LocationWithLength bedLocusSpan{0, 0};
		bedLocusSpan.start  = locusInd;
		bedLocusSpan.length = locusGroupStats.nBytesPerLocus;
		locusInd            = bed2ophThreaded_(bedChunkToRead, threadRanges, bedLocusSpan, permutation, padIndiv);
		locusInd           += excessLoci;
	}
	return locusInd;
}

void GenoTableHash::mac2ophBlk_(const std::vector<int> &macData, const LocationWithLength &locusBlock, const size_t &randVecLen, const std::vector<size_t> &permutation) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	RanDraw prng;
	constexpr uint8_t middleMask{0b10000011};
	constexpr uint8_t endTwoBitMask{0b00000011};
	auto remainderInd          = static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
	remainderInd               = byteSize_ - remainderInd;
	// Create a vector to store random bytes for stochastic heterozygote resolution
	std::vector<uint64_t> rand(randVecLen);
	auto *randBytes          = reinterpret_cast<uint8_t*>( rand.data() );
	const size_t endLocusInd = locusBlock.start + locusBlock.length;
	for (size_t iLocus = locusBlock.start; iLocus < endLocusInd; ++iLocus) {
		// Fill the random byte vector
		for (auto &randValue : rand) {
			randValue = prng.ranInt();
		}
		size_t iIndiv         = 0;
		const size_t begIndiv = iLocus * nIndividuals_;
		std::vector<uint8_t> missMasks(locusSize_, 0);
		std::vector<uint8_t> binLocus(locusSize_, 0);
		size_t i0Byte = 0;                                                                         // to index the missMasks vector
		for (size_t iByte = 0; iByte < locusSize_ - 1; ++iByte) {                                  // treat the last byte separately
			for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte) {
				auto curIndiv             = static_cast<uint8_t>(macData[begIndiv + iIndiv]);      // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= middleMask;                                            // mask everything in the middle
				const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
				missMasks[i0Byte]        |= static_cast<uint8_t>(missingMask << iInByte);
				curIndiv                 &= endTwoBitMask;
				const uint8_t randMask    = (randBytes[iByte] >> iInByte) & oneBit_;               // 0b00000000 or 0b00000001 with equal chance
				uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
				curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
				binLocus[iByte]          |= static_cast<uint8_t>(curBitMask << iInByte);
				++iIndiv;
			}
			++i0Byte;
		}
		// now deal with the last byte in the individual
		for (uint8_t iRem = 0; iRem < remainderInd; ++iRem) {
			auto curIndiv             = static_cast<uint8_t>(macData[begIndiv + iIndiv]);          // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= middleMask;                                                // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
			missMasks.back()         |= static_cast<uint8_t>(missingMask << iRem);
			curIndiv                 &= endTwoBitMask;                                                
			const uint8_t randMask    = (randBytes[locusSize_ - 1] >> iRem) & oneBit_;             // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                   // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                              // zero it out if missing value is set

			binLocus.back()          |= static_cast<uint8_t>(curBitMask << iRem);
			++iIndiv;
		}
		const float maf = static_cast<float>( countSetBits(binLocus) ) / static_cast<float>(nIndividuals_);
		if (maf > 0.5F) { // always want the alternative to be the minor allele
			for (size_t iBL = 0; iBL < locusSize_; ++iBL) {
				binLocus[iBL] = (~binLocus[iBL]) & (~missMasks[iBL]);
			}
			binLocus.back() &= lastByteMask; // unset the remainder bits
		}
		locusOPH_(iLocus, permutation, binLocus);
	}
}

SimilarityMatrix GenoTableHash::hashJacBlock_(const std::pair<RowColIdx, RowColIdx> &blockRange, const std::vector<uint32_t> &locusIndexes, const float &similarityCutOff) const {
	SimilarityMatrix result;
	uint32_t iRow{blockRange.first.iRow};
	if (blockRange.first.iRow != blockRange.second.iRow) {
		for (uint32_t jCol = blockRange.first.jCol; jCol < iRow; ++jCol) { // first, possibly incomplete, row
			RowColIdx localRC{};
			localRC.iRow = locusIndexes[iRow];
			localRC.jCol = locusIndexes[jCol];
			JaccardPair localJP{makeJaccardPair_(localRC)};
			if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
				result.insert(localRC, localJP);
			}
		}
		++iRow;
		while ( iRow < (blockRange.second.iRow) ) { // complete triangle
			for (uint32_t jCol = 0; jCol < iRow; ++jCol) {
				RowColIdx localRC{};
				localRC.iRow = locusIndexes[iRow];
				localRC.jCol = locusIndexes[jCol];
				JaccardPair localJP{makeJaccardPair_(localRC)};
				if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
					result.insert(localRC, localJP);
				}
			}
			++iRow;
		}
		for (uint32_t jColRem = 0; jColRem < blockRange.second.jCol; ++jColRem) {  // last, possibly incomplete, row
			RowColIdx localRC{};
			localRC.iRow = locusIndexes[iRow];
			localRC.jCol = locusIndexes[jColRem];
			JaccardPair localJP{makeJaccardPair_(localRC)};
			if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
				result.insert(localRC, localJP);
			}
		}
		return result;
	}
	for (uint32_t jCol = blockRange.first.jCol; jCol < blockRange.second.jCol; ++jCol) {  // last, possibly incomplete, row
		RowColIdx localRC{};
		localRC.iRow = locusIndexes[iRow];
		localRC.jCol = locusIndexes[jCol];
		JaccardPair localJP{makeJaccardPair_(localRC)};
		if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
			result.insert(localRC, localJP);
		}
	}

	return result;
}

SimilarityMatrix GenoTableHash::hashJacBlock_(const std::pair<HashGroupItPairCount, HashGroupItPairCount> &blockRange, const float &similarityCutOff) const {
	if (blockRange.first.hgIterator == blockRange.second.hgIterator) { // block falls entirely within a group
		std::pair<RowColIdx, RowColIdx> rowColumnPair{};
		rowColumnPair.first  = recoverRCindexes(blockRange.first.pairCount);
		rowColumnPair.second = recoverRCindexes(blockRange.second.pairCount);

		return hashJacBlock_(rowColumnPair, blockRange.first.hgIterator->locusIndexes, similarityCutOff);
	}
	std::pair<RowColIdx, RowColIdx> rowColumnPair{};
	rowColumnPair.first = recoverRCindexes(blockRange.first.pairCount);
	const size_t nPairs{blockRange.first.hgIterator->locusIndexes.size() * (blockRange.first.hgIterator->locusIndexes.size() - 1) / 2};
	rowColumnPair.second = recoverRCindexes(nPairs);
	SimilarityMatrix result{hashJacBlock_(rowColumnPair, blockRange.first.hgIterator->locusIndexes, similarityCutOff)};
	// process complete groups
	std::for_each(
		blockRange.first.hgIterator + 1,
		blockRange.second.hgIterator,
		[this, &result, &similarityCutOff](const HashGroup &eachGroup) {
			std::pair<RowColIdx, RowColIdx> localRCPair{};
			localRCPair.first.iRow = 1;
			localRCPair.first.jCol = 0;
			const size_t locNpairs{eachGroup.locusIndexes.size() * (eachGroup.locusIndexes.size() - 1) / 2};
			localRCPair.second = recoverRCindexes(locNpairs);
			// TODO: eliminate after debugging
			SimilarityMatrix tmp{hashJacBlock_(localRCPair, eachGroup.locusIndexes, similarityCutOff)};
			result.merge( std::move(tmp) );
		}
	);
	// last, possibly incomplete, group
	rowColumnPair.first.iRow = 1;
	rowColumnPair.first.jCol = 0;
	rowColumnPair.second     = recoverRCindexes(blockRange.second.pairCount);
	result.merge( hashJacBlock_(rowColumnPair, blockRange.second.hgIterator->locusIndexes, similarityCutOff) );

	return result;
}

SimilarityMatrix GenoTableHash::hashJacThreaded_(const std::vector< std::pair<RowColIdx, RowColIdx> > &indexPairs, const std::vector<uint32_t> &locusIndexes, const float &similarityCutOff) const {
	std::vector<SimilarityMatrix> threadResults( indexPairs.size() );
	std::vector< std::future<void> > tasks;
	tasks.reserve( indexPairs.size() );
	size_t iThread{0};
	for (const auto &eachRange : indexPairs) {
		tasks.emplace_back(
			std::async([this, eachRange, &locusIndexes, &similarityCutOff, iThread, &threadResults]{
				threadResults.at(iThread) = hashJacBlock_(eachRange, locusIndexes, similarityCutOff);
			})
		);
		++iThread;
	}

	for (const auto &eachThread : tasks) {
		eachThread.wait();
	}

	std::for_each(
		threadResults.begin() + 1,
		threadResults.end(),
		[&threadResults](SimilarityMatrix &eachMatrix) {
			threadResults.at(0).merge( std::move(eachMatrix) );
		}
	);

	return threadResults.at(0);
}

SimilarityMatrix GenoTableHash::hashJacThreaded_(const std::vector< std::pair<HashGroupItPairCount, HashGroupItPairCount> > &blockRanges, const float &similarityCutOff) const {
	std::vector<SimilarityMatrix> threadResults( blockRanges.size() );
	threadResults.at(0) = hashJacBlock_(blockRanges.at(1), similarityCutOff);
	/*
	std::vector< std::future<void> > tasks;
	tasks.reserve( blockRanges.size() );
	size_t iThread{0};
	for (const auto &eachRange : blockRanges) {
		tasks.emplace_back(
			std::async([this, eachRange, iThread, &threadResults, &similarityCutOff]{
				threadResults.at(iThread) = hashJacBlock_(eachRange, similarityCutOff);
			})
		);
		++iThread;
	}

	for (const auto &eachThread : tasks) {
		eachThread.wait();
	}

	std::for_each(
		threadResults.begin() + 1,
		threadResults.end(),
		[&threadResults](SimilarityMatrix &eachMatrix) {
			threadResults.at(0).merge( std::move(eachMatrix) );
		}
	);
	*/

	return threadResults.at(0);
}

JaccardPair GenoTableHash::makeJaccardPair_(const RowColIdx &rowColumn) const noexcept {
	const auto kSkDst{static_cast<std::vector<uint16_t>::difference_type>(kSketches_)};
	const auto start1 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(rowColumn.iRow) * kSketches_;
	const auto start2 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(rowColumn.jCol) * kSketches_;
	// count equal elements using the inner_product idiom
	const int simVal = std::inner_product( start1, start1 + kSkDst, start2, 0, std::plus<>(), std::equal_to<>() );
	JaccardPair localJP{};
	localJP.nIntersect = static_cast<uint32_t>(simVal);
	localJP.nUnion     = static_cast<uint32_t>(kSketches_);
	return localJP;
}
