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

#include <iostream>

#include "gvarHash.hpp"
#include "vashFunctions.hpp"
#include "random.hpp"

using namespace BayesicSpace;

// GenoTableBin methods
constexpr size_t   GenoTableBin::nMagicBytes_    = 3;                // number of leading bytes for .bed files
constexpr uint8_t  GenoTableBin::oneBit_         = 0b00000001;       // One set bit for masking
constexpr uint8_t  GenoTableBin::byteSize_       = 8;                // Size of one byte in bits
constexpr uint8_t  GenoTableBin::bedGenoPerByte_ = 4;                // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableBin::llWordSize_     = 8;                // 64 bit word size in bytes
constexpr size_t   GenoTableBin::maxNlocusPairs_ = 6074000999;       // approximate number of loci that does not overflow with n*(n-1)/2

// Constructors
GenoTableBin::GenoTableBin(const std::string &inputFileName, const size_t &nIndividuals, std::string logFileName, const size_t &nThreads)
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
	const size_t nBedBytesPerLocus = nIndividuals_ / bedGenoPerByte_ + static_cast<size_t>( (nIndividuals_ % bedGenoPerByte_) > 0);
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
	nLoci_ = nTotalBedBytes / nBedBytesPerLocus;
	if ( nLoci_ > std::numeric_limits<uint32_t>::max() ) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + "\n";
		throw std::string("ERROR: there must be fewer than 2^32 loci in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
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
	const size_t ramSize                = getAvailableRAM() / 2UL;                               // measuring here, after all the major allocations; use half to leave resources for other operations
	locusGroupAttributes.nLociToRead    = std::min(ramSize / nBedBytesPerLocus, nLoci_);         // number of .bed loci to read at a time
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

GenoTableBin::GenoTableBin(const std::vector<int> &maCounts, const size_t &nIndividuals, std::string logFileName, const size_t &nThreads)
							: nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals}, nThreads_{nThreads}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime = std::time(nullptr);
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype binarization from minor allele count vector started on " + logStream.str() + "\n";
	logStream.clear();
	if ( nLoci_ > std::numeric_limits<uint32_t>::max() ) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + "\n";
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

GenoTableBin::GenoTableBin(GenoTableBin &&toMove) noexcept : nIndividuals_{0}, nLoci_{0}, binLocusSize_{0}, nThreads_{0} {
	*this = std::move(toMove);
}

GenoTableBin& GenoTableBin::operator=(GenoTableBin &&toMove) noexcept {
	if (this != &toMove) {
		binGenotypes_ = std::move(toMove.binGenotypes_);
		nIndividuals_ = toMove.nIndividuals_;
		nLoci_        = toMove.nLoci_;
		binLocusSize_ = toMove.binLocusSize_;
		nThreads_     = toMove.nThreads_;
		logMessages_  = std::move(toMove.logMessages_);
		logFileName_  = std::move(toMove.logFileName_);

		toMove.nIndividuals_ = 0;
		toMove.nLoci_        = 0;
		toMove.binLocusSize_ = 0;
		toMove.nThreads_     = 0;
	}
	return *this;
}

void GenoTableBin::saveGenoBinary(const std::string &outFileName) const {
	std::fstream out;
	assert( ( binGenotypes_.size() < std::numeric_limits<std::streamsize>::max() ) // NOLINT
			&& "ERROR: binGenotypes_ size exceeds maximum streamsize in GenoTableBin.saveGenoBinary()");
	out.open(outFileName, std::ios::out | std::ios::binary | std::ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), static_cast<std::streamsize>( binGenotypes_.size() ) ); // OK because we are casting to const char*
	out.close();
}

std::vector<IndexedPairLD> GenoTableBin::allJaccardLD() const {
	if (nLoci_ > maxNlocusPairs_) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxNlocusPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	logMessages_ += "Calculating all pairwise LD using full Jaccard similarity estimates and passing the result to the calling function\n";
	const size_t nPairs{nLoci_ * ( nLoci_ - static_cast<size_t>(1) ) / static_cast<size_t>(2)};
	std::vector<IndexedPairLD> LDmat(nPairs);
	const size_t nLocusPairsPerThread = std::max( LDmat.size() / nThreads_, static_cast<size_t>(1) );
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = nLocusPairsPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	threadRanges.back().second = LDmat.size();
	jaccardThreaded_(threadRanges, 0, LDmat);
	return LDmat;
}

void GenoTableBin::allJaccardLD(const std::string &ldFileName) const {
	if (nLoci_ > maxNlocusPairs_) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxNlocusPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	const size_t maxInRAM             = getAvailableRAM() / ( static_cast<size_t>(2) * sizeof(IndexedPairLD) );      // use half to leave resources for other operations
	const size_t nPairs               = nLoci_ * ( nLoci_ - static_cast<size_t>(1) ) / static_cast<size_t>(2);
	const size_t nChunks              = std::max( nPairs / maxInRAM, static_cast<size_t>(1) );
	const size_t chunkSize            = std::min(nPairs, maxInRAM);
	const size_t remainingPairs       = nPairs % nChunks;
	const size_t nLocusPairsPerThread = std::max( chunkSize / nThreads_, static_cast<size_t>(1) );

	logMessages_ += "Calculating all pairwise LD using full Jaccard similarity estimates\n";
	logMessages_ += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "; ";
	logMessages_ += "calculating in " + std::to_string(nChunks) + " chunk(s)\n";

	std::fstream output;
	output.open(ldFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\trSq\n";
	size_t overallPairInd{0};
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = nLocusPairsPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	const size_t excessLoci    = chunkSize - threadRanges.back().second;
	threadRanges.back().second = chunkSize;
	for (size_t iChunk = 0; iChunk < nChunks; ++iChunk) {
		std::vector<IndexedPairLD> LDmatChunk(chunkSize);
		overallPairInd  = jaccardThreaded_(threadRanges, overallPairInd, LDmatChunk);
		overallPairInd += excessLoci;
		saveValues(LDmatChunk, output);
	}
	if (remainingPairs > 0) {
		std::vector<IndexedPairLD> LDmatChunk(remainingPairs);
		const size_t nRemainPairsPerThread = std::max( remainingPairs / nThreads_, static_cast<size_t>(1) );
		threadCounts.count                 = nThreads_;
		threadCounts.size                  = nRemainPairsPerThread;
		threadRanges                       = makeThreadRanges(threadCounts);
		threadRanges.back().second         = remainingPairs;
		overallPairInd                     = jaccardThreaded_(threadRanges, overallPairInd, LDmatChunk);
		saveValues(LDmatChunk, output);
	}
	output.close();
}

void GenoTableBin::allJaccardLD(const InOutFileNames &bimAndLDnames) const {
	if (nLoci_ > maxNlocusPairs_) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxNlocusPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	const size_t maxInRAM       = getAvailableRAM() / ( static_cast<size_t>(2) * sizeof(IndexedPairLD) );      // use half to leave resources for other operations
	const size_t nPairs         = nLoci_ * ( nLoci_ - static_cast<size_t>(1) ) / static_cast<size_t>(2);
	const size_t nChunks        = std::max( nPairs / maxInRAM, static_cast<size_t>(1) );
	const size_t chunkSize      = std::min(nPairs, maxInRAM);
	const size_t remainingPairs = nPairs % nChunks;

	logMessages_ += "Calculating all pairwise LD using full Jaccard similarity estimates\n";
	logMessages_ += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "; ";
	logMessages_ += "calculating in " + std::to_string(nChunks) + " chunk(s)\n";

	logMessages_ += "Getting locus names from the " + bimAndLDnames.inputFileName + " .bim file\n";
	std::vector<std::string> locusNames{getLocusNames(bimAndLDnames.inputFileName)};
	assert( (locusNames.size() == nLoci_) // NOLINT
			&& "ERROR: number of loci in the .bim file not the same as nLoci_");
	logMessages_ += "Read " + std::to_string( locusNames.size() ) + " locus names from the .bim file\n";

	std::fstream output;
	output.open(bimAndLDnames.outputFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\trSq\n";
	const size_t nLocusPairsPerThread = std::max( chunkSize / nThreads_, static_cast<size_t>(1) );
	size_t overallPairInd{0};
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = nLocusPairsPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	const size_t excessLoci    = maxInRAM - threadRanges.back().second;
	threadRanges.back().second = maxInRAM;
	for (size_t iChunk = 0; iChunk < nChunks; ++iChunk) {
		std::vector<IndexedPairLD> LDmatChunk(maxInRAM);
		overallPairInd  = jaccardThreaded_(threadRanges, overallPairInd, LDmatChunk);
		overallPairInd += excessLoci;
		saveValues(LDmatChunk, locusNames, output);
	}
	if (remainingPairs > 0) {
		std::vector<IndexedPairLD> LDmatChunk(remainingPairs);
		const size_t nRemainPairsPerThread = std::max( remainingPairs / nThreads_, static_cast<size_t>(1) );
		threadCounts.count                 = nThreads_;
		threadCounts.size                  = nRemainPairsPerThread;
		threadRanges                       = makeThreadRanges(threadCounts);
		threadRanges.back().second         = remainingPairs;
		overallPairInd                     = jaccardThreaded_(threadRanges, overallPairInd, LDmatChunk);
		saveValues(LDmatChunk, locusNames, output);
	}
	output.close();
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
	uint64_t locSeed{0};
	{
		std::lock_guard<std::mutex> lockGuard(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
	size_t begByte{locusSpan.start * binLocusSize_};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus) {
		LocationWithLength bedWindow{0, 0};
		bedWindow.start  = iBedLocus * locusSpan.length;
		bedWindow.length = locusSpan.length;
		LocationWithLength binWindow{0, 0};
		binWindow.start  = begByte;
		binWindow.length = binLocusSize_;
		binarizeBedLocus(bedWindow, bedData, nIndividuals_, locPRNG, binWindow, binGenotypes_);
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
	uint64_t locSeed{0};
	{
		std::lock_guard<std::mutex> lockGuard(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
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

void GenoTableBin::jaccardBlock_(const std::pair<size_t, size_t> &blockVecRange, const size_t &blockStartAll, std::vector<IndexedPairLD> &ldVec) const {
	const size_t nnLoci = nLoci_ * (nLoci_ - 1) / 2 - 1; // overflow checked in the calling function; do this here for encapsulation, may move to the calling function later
	std::vector<uint8_t> locus(binLocusSize_);
	size_t curJacMatInd{blockStartAll};
	const auto fIndiv{static_cast<float>(nIndividuals_)};
	for (size_t iVecInd = blockVecRange.first; iVecInd < blockVecRange.second; ++iVecInd) {
		// compute row and column indexes from the vectorized lower triangle index
		// got these expressions by combining various web sources and verifying
		const size_t kpIdx  = nnLoci - curJacMatInd;
		const size_t pIdx   = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kpIdx) ) ) - 1) / 2;
		const auto row      = static_cast<uint32_t>(nLoci_ - 2 - pIdx);
		const auto col      = static_cast<uint32_t>(nLoci_ - (kpIdx - pIdx * (pIdx + 1) / 2) - 1);
		const size_t rowBin = row * binLocusSize_;
		const size_t colBin = col * binLocusSize_;
		for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc) {
			locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] & binGenotypes_[colBin + iBinLoc];
		}
		const uint64_t isect = countSetBits(locus);
		for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc) {
			locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] | binGenotypes_[colBin + iBinLoc];
		}
		const LocationWithLength binRowWindow{rowBin, binLocusSize_};
		const LocationWithLength binColWindow{colBin, binLocusSize_};
		const uint64_t uni     = countSetBits(locus);
		const uint64_t locus1n = countSetBits(binGenotypes_, binRowWindow);
		const uint64_t locus2n = countSetBits(binGenotypes_, binColWindow);
		// should be safe: each thread accesses different vector elements
		// if isect is 0, union must also be 0, so we set the distance to 0.0
		ldVec[iVecInd].element1ind = row;
		ldVec[iVecInd].element2ind = col;
		const auto fsect{static_cast<float>(isect)}; 
		const float locus1p{static_cast<float>(locus1n) / fIndiv};                                  // locus 1 allele frequency
		const float locus2p{static_cast<float>(locus2n) / fIndiv};                                  // locus 2 allele frequency
		const float pApB{locus1p * locus2p};
		const float dValue{fsect / fIndiv - pApB};                                                  // D statistic
		ldVec[iVecInd].jaccard = (isect > 0 ? fsect / static_cast<float>(uni) : 0.0F);
		float tmpRsq{dValue * dValue / ( pApB * (1.0F - locus1p) * (1.0F - locus2p) )};
		// normalize; will set any inf values (due to locus being mis-called as monomorphic when all hets are set to '0')
		tmpRsq             = std::min(tmpRsq, 1.0F);
		ldVec[iVecInd].rSq = std::max(tmpRsq, 0.0F);                                                
		++curJacMatInd;
	}
}

size_t GenoTableBin::jaccardThreaded_(const std::vector< std::pair<size_t, size_t> > &pairIndRanges, const size_t &blockStartAll, std::vector<IndexedPairLD> &ldVec) const {
	size_t overallPairInd = blockStartAll;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : pairIndRanges) {
		tasks.emplace_back(
			std::async([this, &ldVec, &eachTR, overallPairInd]{
				jaccardBlock_(eachTR, overallPairInd, ldVec);
			})
		);
		overallPairInd += eachTR.second - eachTR.first;
	}
	for (const auto &eachThread : tasks) {
		eachThread.wait();
	}
	return overallPairInd;
}

// GenoTableHash methods
constexpr size_t   GenoTableHash::nMagicBytes_    = 3;                                    // number of leading bytes for .bed files
constexpr uint8_t  GenoTableHash::oneBit_         = 0b00000001;                           // One set bit for masking 
constexpr uint8_t  GenoTableHash::byteSize_       = 8;                                    // Size of one byte in bits 
constexpr uint8_t  GenoTableHash::bedGenoPerByte_ = 4;                                    // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableHash::llWordSize_     = 8;                                    // 64 bit word size in bytes 
constexpr size_t   GenoTableHash::maxPairs_       = 6074000999;                           // approximate maximum number that does not overflow with n*(n-1)/2
constexpr uint64_t GenoTableHash::roundMask_      = 0xfffffffffffffff8;                   // mask for rounding down to nearest whole-byte value
constexpr uint64_t GenoTableHash::allBitsSet_     = std::numeric_limits<uint64_t>::max(); // 64-bit word with all bits set
constexpr size_t   GenoTableHash::wordSizeInBits_ = 64;                                   // 64-bit word size
constexpr uint16_t GenoTableHash::emptyBinToken_  = std::numeric_limits<uint16_t>::max(); // Value corresponding to an empty token 

// Constructors
GenoTableHash::GenoTableHash(const std::string &inputFileName, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, std::string logFileName)
								: kSketches_{indivSketchCounts.kSketches}, fSketches_{static_cast<float>(indivSketchCounts.kSketches)}, nLoci_{0}, nThreads_{nThreads}, logFileName_{std::move(logFileName)} {
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
	sketchSize_   = indivSketchCounts.nIndividuals / kSketches_ + static_cast<size_t>( (indivSketchCounts.nIndividuals % kSketches_) > 0 );
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
	nLoci_ = fileSize / nBedBytes;
	if ( nLoci_ > std::numeric_limits<uint32_t>::max() ) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + "\n";
		throw std::string("ERROR: there must be fewer than 2^32 loci in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	logMessages_ += "Number of individuals: "         + std::to_string(indivSketchCounts.nIndividuals) + "\n";
	logMessages_ += "Number of individuals to hash: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: "                + std::to_string(nLoci_) + "\n";
	logMessages_ += "Hash size: "                     + std::to_string(kSketches_) + "\n";

	locusSize_      = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;                    // round up to the nearest multiple of 8
	nFullWordBytes_ = (nIndividuals_ - 1) / byteSize_;
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	inStream.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStream.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	BedDataStats locusGroupAttributes{};
	locusGroupAttributes.nBytesPerLocus = indivSketchCounts.nIndividuals / bedGenoPerByte_ + static_cast<size_t>(indivSketchCounts.nIndividuals % bedGenoPerByte_ > 0);
	const size_t ramSize                = getAvailableRAM() / 2UL;                                    // measuring here, after all the major allocations; use half to leave resources for other operations
	locusGroupAttributes.nLociToRead    = std::min(ramSize / locusGroupAttributes.nBytesPerLocus, nLoci_);       // number of .bed loci to read at a time
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
		addIndv.emplace_back( iAddIndiv, rng_.sampleInt(indivSketchCounts.nIndividuals) );
	}
	if ( !addIndv.empty() ) {
		std::string addIndexes;
		for (const auto &eachIdx : addIndv) {
			addIndexes += std::to_string(eachIdx.second) + " ";
		}
		logMessages_ += "Re-sampled individuals: " + addIndexes + "\n";
	}
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{rng_.fyIndexesUp(nIndividuals_)};
	std::vector<uint32_t> seeds{static_cast<uint32_t>( rng_.ranInt() )};

	locusGroupAttributes.firstLocusIdx = 0;
	locusGroupAttributes.firstLocusIdx = bed2oph_(locusGroupAttributes, inStream, ranInts, addIndv, seeds);
	std::cout << "remaining loci: " << remainingLoci << "\n";
	if (remainingLoci > 0) {
		locusGroupAttributes.nLociPerThread = std::max(remainingLoci / nThreads_, 1UL);
		locusGroupAttributes.nBytesToRead   = remainingBytes;
		locusGroupAttributes.nLociToRead    = remainingLoci;
		locusGroupAttributes.nMemChunks     = 1;
		bed2oph_(locusGroupAttributes, inStream, ranInts, addIndv, seeds);
	}
	inStream.close();
}

GenoTableHash::GenoTableHash(const std::vector<int> &maCounts, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, std::string logFileName) 
		: nIndividuals_{indivSketchCounts.nIndividuals}, kSketches_{indivSketchCounts.kSketches}, fSketches_{static_cast<float>(indivSketchCounts.kSketches)},
				nLoci_{maCounts.size() / indivSketchCounts.nIndividuals}, nThreads_{nThreads}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime = std::time(nullptr);
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype hashing from a minor allele count vector started on " + logStream.str() + "\n";
	logStream.clear();
	if ( nLoci_ > std::numeric_limits<uint32_t>::max() ) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + "\n";
		throw std::string("ERROR: there must be fewer than 2^32 loci in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
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
	nThreads_     = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_     = std::max(nThreads_, 1UL);
	sketchSize_   = indivSketchCounts.nIndividuals / kSketches_ + static_cast<size_t>( (indivSketchCounts.nIndividuals % kSketches_) > 0 );
	nIndividuals_ = sketchSize_ * kSketches_;
	if (indivSketchCounts.kSketches >= emptyBinToken_) {
		logMessages_ += "ERROR: sketch size (" + std::to_string(indivSketchCounts.kSketches) + ") is too small; aborting\n";
		throw std::string("ERROR: Number of sketches (") + std::to_string(kSketches_) + std::string(") implies sketch size (") +
			std::to_string(indivSketchCounts.kSketches) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	locusSize_      = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;                    // round up to the nearest multiple of 8
	nFullWordBytes_ = (nIndividuals_ - 1) / byteSize_;
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	const size_t ranVecSize = locusSize_ / llWordSize_ + static_cast<size_t>( (locusSize_ % llWordSize_) > 0);
	// Calculate the actual sketch number based on the realized sketch size
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{rng_.fyIndexesUp(nIndividuals_)};
	std::vector<uint32_t> seeds{static_cast<uint32_t>( rng_.ranInt() )};
	seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );

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
				std::async([this, &maCounts, &eachTR, ranVecSize, &ranInts, &seeds]{
					LocationWithLength threadLoci{0, 0};
					threadLoci.start  = eachTR.first;
					threadLoci.length = eachTR.second - eachTR.first;
					mac2ophBlk_(maCounts, threadLoci, ranVecSize, ranInts, seeds);
				})
			);
		}
		for (const auto &eachTask : tasks) {
			eachTask.wait();
		}
	} else {
		LocationWithLength allLoci{0, 0};
		allLoci.length = nLoci_;
		mac2ophBlk_(maCounts, allLoci, ranVecSize, ranInts, seeds);
	}
}

GenoTableHash::GenoTableHash(GenoTableHash &&toMove) noexcept : nIndividuals_{0}, kSketches_{0}, fSketches_{0.0}, sketchSize_{0}, nLoci_{0}, locusSize_{0}, nFullWordBytes_{0}, nThreads_{0} {
	*this = std::move(toMove);
}

GenoTableHash& GenoTableHash::operator=(GenoTableHash &&toMove) noexcept {
	if (this != &toMove) {
		sketches_       = std::move(toMove.sketches_);
		nIndividuals_   = toMove.nIndividuals_;
		kSketches_      = toMove.kSketches_;
		fSketches_      = toMove.fSketches_;
		sketchSize_     = toMove.sketchSize_;
		nLoci_          = toMove.nLoci_;
		locusSize_      = toMove.locusSize_;
		nFullWordBytes_ = toMove.nFullWordBytes_;
		nThreads_       = toMove.nThreads_;
		logFileName_    = std::move(toMove.logFileName_);
		logMessages_    = std::move(toMove.logMessages_);

		toMove.nIndividuals_   = 0;
		toMove.kSketches_      = 0;
		toMove.sketchSize_     = 0;
		toMove.nLoci_          = 0;
		toMove.locusSize_      = 0;
		toMove.nFullWordBytes_ = 0;
		toMove.nThreads_       = 0;
	}
	return *this;
}

std::vector<IndexedPairSimilarity> GenoTableHash::allHashLD() const {
	if (nLoci_ > maxPairs_) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	logMessages_ += "Calculating all pairwise LD and passing the result to the calling function\n";
	const size_t nPairs{nLoci_ * ( nLoci_ - static_cast<size_t>(1) ) / static_cast<size_t>(2)};
	std::vector<IndexedPairSimilarity> LDmat(nPairs);
	const size_t nLocusPairsPerThread = std::max( LDmat.size() / nThreads_, static_cast<size_t>(1) );
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = nLocusPairsPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	threadRanges.back().second = LDmat.size();
	hashJacThreaded_(threadRanges, 0, LDmat);
	return LDmat;
}

void GenoTableHash::allHashLD(const std::string &ldFileName) const {
	if (nLoci_ > maxPairs_) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t maxInRAM             = getAvailableRAM() / ( 2UL * sizeof(IndexedPairSimilarity) );      // use half to leave resources for other operations
	const size_t nPairs               = nLoci_ * (nLoci_ - 1) / 2;
	const size_t chunkSize            = std::min(nPairs, maxInRAM);
	const size_t nChunks              = std::max( nPairs / maxInRAM, static_cast<size_t>(1) );
	const size_t remainingPairs       = nPairs % nChunks;
	const size_t nLocusPairsPerThread = std::max( chunkSize / nThreads_, static_cast<size_t>(1) );

	logMessages_ += "Calculating all pairwise LD\n";
	logMessages_ += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "\n";
	logMessages_ += "calculating in " + std::to_string(nChunks) + " chunk(s)\n";

	size_t overallPairInd{0};
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = nLocusPairsPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	const size_t excessLoci    = chunkSize - threadRanges.back().second;
	threadRanges.back().second = chunkSize;

	std::fstream output;
	output.open(ldFileName, std::ios::trunc | std::ios::out);
	output << "groupID\tlocus1\tlocus2\tjaccLD\n";
	for (size_t iChunk = 0; iChunk < nChunks; ++iChunk) {
		std::vector<IndexedPairSimilarity> LDmatChunk(chunkSize);
		overallPairInd  = hashJacThreaded_(threadRanges, overallPairInd, LDmatChunk);
		overallPairInd += excessLoci;
		saveValues(LDmatChunk, output);
	}
	if (remainingPairs > 0) {
		std::vector<IndexedPairSimilarity> LDmatChunk(remainingPairs);
		const size_t nRemainPairsPerThread = std::max( remainingPairs / nThreads_, static_cast<size_t>(1) );
		threadCounts.count                 = nThreads_;
		threadCounts.size                  = nRemainPairsPerThread;
		threadRanges                       = makeThreadRanges(threadCounts);
		threadRanges.back().second         = remainingPairs;
		overallPairInd                     = hashJacThreaded_(threadRanges, overallPairInd, LDmatChunk);
		saveValues(LDmatChunk, output);
	}
	output.close();
}

void GenoTableHash::allHashLD(const InOutFileNames &bimAndLDnames) const {
	if (nLoci_ > maxPairs_) {
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	logMessages_ += "Getting locus names from the " + bimAndLDnames.inputFileName + " .bim file\n";
	std::vector<std::string> locusNames{getLocusNames(bimAndLDnames.inputFileName)};
	assert( (locusNames.size() == nLoci_) // NOLINT
			&& "ERROR: number of loci in the .bim file not the same as nLoci_");

	const size_t maxInRAM             = getAvailableRAM() / ( 2UL * sizeof(IndexedPairSimilarity) );      // use half to leave resources for other operations
	const size_t nPairs               = nLoci_ * (nLoci_ - 1) / 2;
	const size_t chunkSize            = std::min(nPairs, maxInRAM);
	const size_t nChunks              = std::max( nPairs / maxInRAM, static_cast<size_t>(1) );
	const size_t remainingPairs       = nPairs % nChunks;
	const size_t nLocusPairsPerThread = std::max( chunkSize / nThreads_, static_cast<size_t>(1) );

	logMessages_ += "Calculating all pairwise LD\n";
	logMessages_ += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "\n";
	logMessages_ += "calculating in " + std::to_string(nChunks) + " chunk(s)\n";

	size_t overallPairInd{0};
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = nLocusPairsPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	const size_t excessLoci    = chunkSize - threadRanges.back().second;
	threadRanges.back().second = chunkSize;

	std::fstream output;
	output.open(bimAndLDnames.outputFileName, std::ios::trunc | std::ios::out);
	output << "groupID\tlocus1\tlocus2\tjaccLD\n";
	for (size_t iChunk = 0; iChunk < nChunks; ++iChunk) {
		std::vector<IndexedPairSimilarity> LDmatChunk(chunkSize);
		overallPairInd  = hashJacThreaded_(threadRanges, overallPairInd, LDmatChunk);
		overallPairInd += excessLoci;
		saveValues(LDmatChunk, locusNames, output);
	}
	if (remainingPairs > 0) {
		std::vector<IndexedPairSimilarity> LDmatChunk(remainingPairs);
		const size_t nRemainPairsPerThread = std::max( remainingPairs / nThreads_, static_cast<size_t>(1) );
		threadCounts.count                 = nThreads_;
		threadCounts.size                  = nRemainPairsPerThread;
		threadRanges                       = makeThreadRanges(threadCounts);
		threadRanges.back().second         = remainingPairs;
		overallPairInd                     = hashJacThreaded_(threadRanges, overallPairInd, LDmatChunk);
		saveValues(LDmatChunk, locusNames, output);
	}
	output.close();
}

std::vector< std::vector<uint32_t> > GenoTableHash::makeLDgroups(const size_t &nRowsPerBand) const {
	assert( (nRowsPerBand != 0) // NOLINT
			&& "ERROR: nRowsPerBand must not be 0 in makeLDgroups()" );
	assert( (nRowsPerBand < kSketches_) // NOLINT
			&& "ERROR: nRowsPerBand must be less than kSketches_ in makeLDgroups()" );
	const size_t nBands = kSketches_ / nRowsPerBand;                                                          // only using full-size bands because smaller ones permit inclusion of low-similarity pairs
	assert( ( nBands >= std::numeric_limits<uint16_t>::max() ) // NOLINT
			&& "ERROR: number of bands cannot exceed uint16_t max in makeLDgroups()" );

	logMessages_ += "Grouping loci\n";
	logMessages_ += "Number of rows per band: " + std::to_string(nRowsPerBand) + "\n";
	logMessages_ += "Number of bands: "         + std::to_string(nBands) + "\n";

	const auto sketchSeed = static_cast<uint32_t>( rng_.ranInt() );
	std::unordered_map< uint32_t, std::vector<uint32_t> > ldGroups;                                           // the hash table

	for (uint32_t iLocus = 0; iLocus < nLoci_; ++iLocus) {
		size_t iSketch = 0;
		for (size_t iBand = 0; iBand < nBands; ++iBand) {
			std::vector<uint16_t> bandVec{static_cast<uint16_t>(iBand)};                                      // add the band index to the hash, so that only corresponding bands are compared
			const size_t firstSketchIdx = iSketch + iLocus * kSketches_;                                      // iSketch tracks band IDs
			for (size_t iInBand = firstSketchIdx; iInBand < firstSketchIdx + nRowsPerBand; ++iInBand) {
				bandVec.push_back(sketches_[iInBand]);
			}
			LocationWithLength bandVecWindow{0, 0};
			bandVecWindow.start  = 0;
			bandVecWindow.length = bandVec.size();
			const uint32_t hash  = murMurHash(bandVec, bandVecWindow, sketchSeed);
			ldGroups[hash].push_back(iLocus);
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
	return groups;
}

void GenoTableHash::makeLDgroups(const size_t &nRowsPerBand, const std::string &outFileName) const {
	std::vector< std::vector<uint32_t> > ldGroups{this->makeLDgroups(nRowsPerBand)};
	logMessages_ += "Saving group IDs only\n";
	std::fstream out;
	out.open(outFileName, std::ios::out | std::ios::trunc);
	out << "groupID\tlocusIdx\n";
	uint32_t groupID{1};
	for (const auto &eachGroup : ldGroups) {
		for (const auto &locusIdx : eachGroup) {
			out << "G" << groupID << "\t" << locusIdx + 1 << "\n";
		}
		++groupID;
	}
	out.close();
}

void GenoTableHash::makeLDgroups(const size_t &nRowsPerBand, const InOutFileNames &bimAndGroupNames) const {
	std::vector< std::vector<uint32_t> > ldGroups{this->makeLDgroups(nRowsPerBand)};
	logMessages_ += "Saving group IDs only\n";

	logMessages_ += "Getting locus names from the " + bimAndGroupNames.inputFileName + " .bim file\n";
	std::vector<std::string> locusNames{getLocusNames(bimAndGroupNames.inputFileName)};
	assert( (locusNames.size() == nLoci_) // NOLINT
			&& "ERROR: number of loci in the .bim file not the same as nLoci_");

	std::fstream out;
	out.open(bimAndGroupNames.outputFileName, std::ios::out | std::ios::trunc);
	out << "groupID\tlocusIdx\n";
	uint32_t groupID{1};
	for (const auto &eachGroup : ldGroups) {
		for (const auto &locusIdx : eachGroup) {
			out << "G" << groupID << "\t" << locusNames[locusIdx] << "\n";
		}
		++groupID;
	}
	out.close();
}

void GenoTableHash::ldInGroups(const size_t &nRowsPerBand, const std::string &outFileName) const {
	std::vector< std::vector<uint32_t> > ldGroups{this->makeLDgroups(nRowsPerBand)};
	
	logMessages_         += "Estimating LD in groups\n";
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * sizeof(IndexedPairSimilarity) );                          // use half to leave resources for other operations
	logMessages_         += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "; ";
	std::vector<size_t> groupSizes;                                                                               // number of locus pair in each group
	size_t totalPairNumber{0};                                                                                    // total number of pairs
	for (const auto &eachGrp : ldGroups) {
		assert( (eachGrp.size() > 1) // NOLINT
				&& "ERROR: groups cannot be empty in ldInGroups" );
		const size_t nPairs = eachGrp.size() * (eachGrp.size() - 1) / 2;
		groupSizes.push_back(nPairs);
		totalPairNumber += nPairs;
	}
	logMessages_ += "number of pairs in the hash table: " + std::to_string(totalPairNumber) + "\n";
	std::fstream output;
	output.open(outFileName, std::ios::trunc | std::ios::out);
	output << "groupID\tlocus1\tlocus2\tjaccLD\n";
	auto groupIt   = ldGroups.cbegin();
	auto grpSizeIt = groupSizes.cbegin();
	uint16_t iChunk{1};
	while ( groupIt != ldGroups.cend() ) {
		size_t nPairs{0};
		auto blockEndIt = groupIt;
		while ( (nPairs <= maxInRAM) && ( blockEndIt != ldGroups.cend() ) ) {
			nPairs += *grpSizeIt;
			++grpSizeIt;
			++blockEndIt;
		}
		std::vector<IndexedPairSimilarity> hashJacGroups{
			vectorizeGroups( static_cast<uint32_t>( std::distance(groupIt, blockEndIt) ), groupIt, blockEndIt )
		};
		logMessages_ += "\tChunk " + std::to_string(iChunk) + ":\n";
		logMessages_ += "\tNumber of locus pairs before removing duplicates: " + std::to_string( hashJacGroups.size() ) + "\n";
		std::sort(hashJacGroups.begin(), hashJacGroups.end(),
					[](const IndexedPairSimilarity &first, const IndexedPairSimilarity &second) {
						return (first.element1ind == second.element1ind ? first.element2ind < second.element2ind : first.element1ind < second.element1ind);
					}
				);
		// some locus pairs are present in multiple groups, so we eliminate duplicates
		auto lastUniqueIt = std::unique(hashJacGroups.begin(), hashJacGroups.end(),
					[](const IndexedPairSimilarity &first, const IndexedPairSimilarity &second) {
						return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
					}
				);
		hashJacGroups.erase( lastUniqueIt, hashJacGroups.end() );
		hashJacGroups.shrink_to_fit();
		logMessages_ += "Number of locus pairs after removing duplicates: " + std::to_string( hashJacGroups.size() ) + "\n";
		// estimate locus pair similarities
		const size_t nPairsPerThread = hashJacGroups.size() / nThreads_;
		if (nPairsPerThread > 0) {
			CountAndSize threadCounts{0, 0};
			threadCounts.count = nThreads_;
			threadCounts.size  = nPairsPerThread;
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
			// assuming the number of threads << number of groups, this should not lead to noticeable thread imbalance
			threadRanges.back().second = hashJacGroups.size();
			hashJacThreaded_(threadRanges, hashJacGroups);
		} else {
			hashJacBlock_( hashJacGroups.begin(), hashJacGroups.end() );
		}

		saveValues(hashJacGroups, output);
		groupIt = blockEndIt;
		++iChunk;
	}
	output.close();
}

void GenoTableHash::ldInGroups(const size_t &nRowsPerBand, const InOutFileNames &bimAndLDnames) const {
	std::vector< std::vector<uint32_t> > ldGroups{this->makeLDgroups(nRowsPerBand)};
	
	logMessages_         += "Estimating LD in groups\n";
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * sizeof(IndexedPairSimilarity) );                          // use half to leave resources for other operations
	logMessages_         += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "; ";
	std::vector<size_t> groupSizes;                                                                               // number of locus pair in each group
	size_t totalPairNumber{0};                                                                                    // total number of pairs
	for (const auto &eachGrp : ldGroups) {
		assert( (eachGrp.size() > 1) // NOLINT
				&& "ERROR: groups cannot be empty in ldInGroups" );
		const size_t nPairs = eachGrp.size() * (eachGrp.size() - 1) / 2;
		groupSizes.push_back(nPairs);
		totalPairNumber += nPairs;
	}
	logMessages_ += "number of pairs in the hash table: " + std::to_string(totalPairNumber) + "\n";
	logMessages_ += "Getting locus names from the " + bimAndLDnames.inputFileName + " .bim file\n";
	std::vector<std::string> locusNames{getLocusNames(bimAndLDnames.inputFileName)};
	assert( (locusNames.size() == nLoci_) // NOLINT
			&& "ERROR: number of loci in the .bim file not the same as nLoci_" );
	logMessages_ += "calculating in chunks\n";
	// too many loci to fit the LD matrix in RAM
	// will work on chunks and save as we go
	std::fstream output;
	output.open(bimAndLDnames.outputFileName, std::ios::trunc | std::ios::out);
	output << "groupID\tlocus1\tlocus2\tjaccLD\n";
	auto groupIt   = ldGroups.cbegin();
	auto grpSizeIt = groupSizes.cbegin();
	uint16_t iChunk{1};
	while ( groupIt != ldGroups.cend() ) {
		size_t nPairs{0};
		auto blockEndIt = groupIt;
		while ( (nPairs <= maxInRAM) && ( blockEndIt != ldGroups.cend() ) ) {
			nPairs += *grpSizeIt;
			++grpSizeIt;
			++blockEndIt;
		}
		std::vector<IndexedPairSimilarity> hashJacGroups{
			vectorizeGroups( static_cast<uint32_t>( std::distance(groupIt, blockEndIt) ), groupIt, blockEndIt )
		};
		logMessages_ += "\tChunk " + std::to_string(iChunk) + ":\n";
		logMessages_ += "\tNumber of locus pairs before removing duplicates: " + std::to_string( hashJacGroups.size() ) + "\n";
		std::sort(hashJacGroups.begin(), hashJacGroups.end(),
					[](const IndexedPairSimilarity &first, const IndexedPairSimilarity &second) {
						return (first.element1ind == second.element1ind ? first.element2ind < second.element2ind : first.element1ind < second.element1ind);
					}
				);
		auto lastUniqueIt = std::unique(hashJacGroups.begin(), hashJacGroups.end(),
					[](const IndexedPairSimilarity &first, const IndexedPairSimilarity &second) {
						return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
					}
				);
		hashJacGroups.erase( lastUniqueIt, hashJacGroups.end() );
		hashJacGroups.shrink_to_fit();
		logMessages_ += "Number of locus pairs after removing duplicates: " + std::to_string( hashJacGroups.size() ) + "\n";
		// estimate locus pair similarities
		const size_t nPairsPerThread = hashJacGroups.size() / nThreads_;
		if (nPairsPerThread > 0) {
			CountAndSize threadCounts{0, 0};
			threadCounts.count = nThreads_;
			threadCounts.size  = nPairsPerThread;
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
			// assuming the number of threads << number of groups, this should not lead to noticeable thread imbalance
			threadRanges.back().second = hashJacGroups.size();
			hashJacThreaded_(threadRanges, hashJacGroups);
		} else {
			hashJacBlock_( hashJacGroups.begin(), hashJacGroups.end() );
		}

		saveValues(hashJacGroups, locusNames, output);
		groupIt = blockEndIt;
		++iChunk;
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

void GenoTableHash::locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds, std::vector<uint8_t> &binLocus) {
	// Start with a permutation to make OPH
	permuteBits_(permutation, binLocus);
	// Now make the sketches
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
		filledIndexes.push_back( rng_.sampleInt(kSketches_) );
	}
	size_t emptyCount = kSketches_ - filledIndexes.size();
	while (emptyCount > 0) {
		for (const auto eachFI : filledIndexes) {
			std::array<uint32_t, SIZE_OF_SIZET> key{};
			memcpy(key.data(), &eachFI, sketchSize_);
			auto newIdx = static_cast<uint32_t>(murMurHash(key, seeds[iSeed]) % kSketches_ + sketchBeg);
			// should be safe: each thread accesses different vector elements
			if (sketches_[newIdx] == emptyBinToken_) {
				sketches_[newIdx] = sketches_[eachFI + sketchBeg];
				--emptyCount;
				break;
			}
		}
		++iSeed;
		std::lock_guard<std::mutex> lock(mtx_);      // lock before measuring to ensure that the size is valid
		if ( iSeed == seeds.size() ) {
			seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
		}
	}
}

void GenoTableHash::bed2ophBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t>&bedLocusIndRange, const LocationWithLength &bedLocusSpan,
									const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint64_t locSeed{0};
	{
		std::lock_guard<std::mutex> lock(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
	size_t iLocus{bedLocusSpan.start};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus) {
		std::vector<uint8_t> binLocus(locusSize_, 0);
		LocationWithLength bedWindow{0, 0};
		bedWindow.start  = iBedLocus * bedLocusSpan.length;
		bedWindow.length = bedLocusSpan.length;
		LocationWithLength binWindow{0, 0};
		binWindow.length = locusSize_;
		binarizeBedLocus(bedWindow, bedData, nIndividuals_, locPRNG, binWindow, binLocus);
		// pad the locus to have an whole number of sketches 
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
		locusOPH_(iLocus, permutation, seeds, binLocus);
		++iLocus;
	}
}

size_t GenoTableHash::bed2ophThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const LocationWithLength &bedLocusSpan,
											const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds) {
	size_t locusInd = bedLocusSpan.start;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges) {
		LocationWithLength currentBedLocusSpan{0, 0};
		currentBedLocusSpan.start  = locusInd;
		currentBedLocusSpan.length = bedLocusSpan.length;
		tasks.emplace_back(
			std::async([this, &bedData, &eachTR, &currentBedLocusSpan, &permutation, &padIndiv, &seeds]{
				bed2ophBlk_(bedData, eachTR, currentBedLocusSpan, permutation, padIndiv, seeds);
			})
		);
		locusInd += eachTR.second - eachTR.first;
	}
	for (const auto &eachTask : tasks) {
		eachTask.wait();
	}
	return locusInd;
}

size_t GenoTableHash::bed2oph_(const BedDataStats &locusGroupStats, std::fstream &bedStream, const std::vector<size_t> &permutation,
								const std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds) {
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
		locusInd            = bed2ophThreaded_(bedChunkToRead, threadRanges, bedLocusSpan, permutation, padIndiv, seeds);
		locusInd           += excessLoci;
	}
	return locusInd;
}

void GenoTableHash::mac2ophBlk_(const std::vector<int> &macData, const LocationWithLength &locusBlock, const size_t &randVecLen, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint64_t locSeed{0};
	{
		std::lock_guard<std::mutex> lock(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
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
			randValue = locPRNG.ranInt();
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
		locusOPH_(iLocus, permutation, seeds, binLocus);
	}
}

void GenoTableHash::hashJacBlock_(const std::pair<size_t, size_t> &blockRange, const size_t &blockStartAll, std::vector<IndexedPairSimilarity> &hashJacVec) const {
	const size_t nnLoci = nLoci_ * (nLoci_ - 1) / 2 - 1; // overflow checked in the calling function; do this here for encapsulation, may move to the calling function later
	size_t curJacMatInd = blockStartAll;
	const auto kSkDst   = static_cast<std::vector<uint16_t>::difference_type>(kSketches_);
	for (size_t iVecInd = blockRange.first; iVecInd < blockRange.second; ++iVecInd) {
		// compute row and column indexes from the vectorized by column lower triangle index
		// got these expressions by combining various web sources and verifying
		const size_t kpIdx = nnLoci - curJacMatInd;
		const size_t pIdx  = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kpIdx) ) ) - 1) / 2;
		const auto row     = static_cast<uint32_t>(nLoci_ - 2 - pIdx);
		const auto col     = static_cast<uint32_t>(nLoci_ - (kpIdx - pIdx * (pIdx + 1) / 2) - 1);
		const auto rowSk   = static_cast<std::vector<uint16_t>::difference_type>(row * kSketches_);
		const auto colSk   = static_cast<std::vector<uint16_t>::difference_type>(col * kSketches_);
		auto start         = sketches_.begin() + rowSk;
		// count equal elements using the inner_product idiom
		int simVal = std::inner_product( start, start + kSkDst, sketches_.begin() + colSk, 0, std::plus<>(), std::equal_to<>() );
		// should be safe: each thread accesses different vector elements
		hashJacVec[iVecInd].groupID         = 0;
		hashJacVec[iVecInd].element1ind     = row;
		hashJacVec[iVecInd].element2ind     = col;
		hashJacVec[iVecInd].similarityValue = static_cast<float>(simVal) / fSketches_;
		++curJacMatInd;
	}
}

void GenoTableHash::hashJacBlock_(const LocationWithLength &blockRange, std::vector<IndexedPairSimilarity> &indexedJacc) const {
	const auto kSkDst = static_cast<std::vector<uint16_t>::difference_type>(kSketches_);
	const size_t blockEnd{blockRange.start + blockRange.length};
	for (size_t iBlock = blockRange.start; iBlock < blockEnd; ++iBlock) {
		const auto start1 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(indexedJacc[iBlock].element1ind * kSketches_);
		const auto start2 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(indexedJacc[iBlock].element2ind * kSketches_);
		// count equal elements using the inner_product idiom
		int simVal = std::inner_product( start1, start1 + kSkDst, start2, 0, std::plus<>(), std::equal_to<>() );
		// should be safe: each thread accesses different vector elements
		indexedJacc[iBlock].similarityValue = static_cast<float>(simVal) / fSketches_;
	}
}

void GenoTableHash::hashJacBlock_(const std::vector<IndexedPairSimilarity>::iterator blockStart, const std::vector<IndexedPairSimilarity>::iterator blockEnd) const {
	const auto kSkDst = static_cast<std::vector<uint16_t>::difference_type>(kSketches_);
	for (auto ipsIt = blockStart; ipsIt != blockEnd; ++ipsIt) {
		const auto start1 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(ipsIt->element1ind * kSketches_);
		const auto start2 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(ipsIt->element2ind * kSketches_);
		// count equal elements using the inner_product idiom
		int simVal = std::inner_product( start1, start1 + kSkDst, start2, 0, std::plus<>(), std::equal_to<>() );
		// should be safe: each thread accesses different vector elements
		ipsIt->similarityValue = static_cast<float>(simVal) / fSketches_;
	}
}

size_t GenoTableHash::hashJacThreaded_(const std::vector< std::pair<size_t, size_t> > &threadRanges, const size_t &blockStartAll, std::vector<IndexedPairSimilarity> &hashJacVec) const {
	size_t overallPairInd = blockStartAll;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges) {
		tasks.emplace_back(
			std::async([this, &hashJacVec, &eachTR, overallPairInd]{
				hashJacBlock_(eachTR, overallPairInd, hashJacVec);
			})
		);
		overallPairInd += eachTR.second - eachTR.first;
	}
	for (const auto &eachTask : tasks) {
		eachTask.wait();
	}
	return overallPairInd;
}

void GenoTableHash::hashJacThreaded_(const std::vector< std::pair<size_t, size_t> > &threadRanges, std::vector<IndexedPairSimilarity> &hashJacVec) const {
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges) {
		const auto hjvStartIt = hashJacVec.begin() + static_cast<std::vector<IndexedPairSimilarity>::difference_type>(eachTR.first);
		const auto hjvEndIt   = hashJacVec.begin() + static_cast<std::vector<IndexedPairSimilarity>::difference_type>(eachTR.second);
		tasks.emplace_back(
			std::async([this, hjvStartIt, hjvEndIt]{
				hashJacBlock_(hjvStartIt, hjvEndIt);
			})
		);
	}
	for (const auto &eachTask : tasks) {
		eachTask.wait();
	}
}
