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
	if (nIndividuals <= 1){
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nIndividuals > std::numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too big to make a square relationship matrix; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: the number of individuals (") + std::to_string(nIndividuals) + std::string( ") is too big to make a square relationship matrix in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ){
		nThreads_ = std::thread::hardware_concurrency();
	}
	const size_t nBedBytesPerLocus = nIndividuals_ / bedGenoPerByte_ + static_cast<size_t>( (nIndividuals_ % bedGenoPerByte_) > 0);
	std::fstream inStr;
	// Start by measuring file size
	inStr.open(inputFileName, std::ios::in | std::ios::binary | std::ios::ate);
	if ( inStr.fail() ){
		logMessages_ += "ERROR: failed to open file " + inputFileName + "\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const uint64_t endPosition = static_cast<uint64_t>( inStr.tellg() );
	if (endPosition < nMagicBytes_){
		logMessages_ += "ERROR: no genotype records in file " + inputFileName + "\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t nTotalBedBytes{static_cast<uint64_t>(endPosition) - nMagicBytes_};
	inStr.close();
	nLoci_        = nTotalBedBytes / nBedBytesPerLocus;
	logMessages_ += "Number of individuals: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: " + std::to_string(nLoci_) + "\n";
	logMessages_ += "Number of threads: " + std::to_string(nThreads_) + "\n";

	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStr.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	binLocusSize_ = nIndividuals_ / byteSize_ + static_cast<size_t>( (nIndividuals_ % byteSize_) > 0 );
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);
	const size_t ramSize        = getAvailableRAM() / 2UL;                               // measuring here, after all the major allocations; use half to leave resources for other operations
	size_t nBedLociToRead       = ramSize / nBedBytesPerLocus;                           // number of .bed loci to read at a time
	nBedLociToRead              = (nBedLociToRead < nLoci_ ? nBedLociToRead : nLoci_);
	const size_t remainingLoci  = nLoci_ % nBedLociToRead;
	const size_t remainingBytes = remainingLoci * nBedBytesPerLocus;
	const size_t nChunks        = nLoci_ / nBedLociToRead;
	size_t nBedBytesToRead      = nBedLociToRead * nBedBytesPerLocus;
	nBedBytesToRead             = (nBedBytesToRead > std::numeric_limits<std::streamsize>::max() ? std::numeric_limits<std::streamsize>::max() : nBedBytesToRead);
	size_t nLociPerThread       = nBedLociToRead / nThreads_;
	logMessages_               += "RAM available for reading the .bed file: " + std::to_string(ramSize) + " bytes\n";
	logMessages_               += ".bed file will be read in " + std::to_string(nChunks) + " chunk(s)\n";
	std::vector<char> bedChunkToRead(nBedBytesToRead, 0);
	assert( ( remainingBytes < std::numeric_limits<std::streamsize>::max() )
			&& "ERROR: remainingBytes larger than maximum streamsize in GenoTableBin constructor");

	size_t locusInd = 0;
	if (nLociPerThread > 0){
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLociPerThread)};
		const size_t excessLoci = nBedLociToRead - threadRanges.back().second;
		threadRanges.back().second = nBedLociToRead;
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			assert( ( nBedBytesToRead < std::streamsize::max() ) && "ERROR: nBedBytesToRead exceeds maximum streamsize in GenoTableBin constructor");
			inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(nBedBytesToRead) );
			locusInd  = bed2binThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus);
			locusInd += excessLoci;
		}
	} else {
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nBedLociToRead, 1)};
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			assert( ( nBedBytesToRead < std::streamsize::max() ) && "ERROR: nBedBytesToRead exceeds maximum streamsize in GenoTableBin constructor");
			inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(nBedBytesToRead) );
			locusInd = bed2binThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus);
		}
	}
	if (remainingLoci > 0) {
		bedChunkToRead.resize(remainingBytes);
		assert( ( remainingBytes < std::streamsize::max() ) && "ERROR: remainingBytes exceeds maximum streamsize in GenoTableBin constructor");
		inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(remainingBytes) );
		nLociPerThread = remainingLoci / nThreads_;
		if (nLociPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLociPerThread)};
			const size_t excessLoci = remainingLoci - threadRanges.back().second;
			threadRanges.back().second = remainingLoci;
			locusInd  = bed2binThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus);
			locusInd += excessLoci;
		} else {
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(remainingLoci, 1)};
			locusInd = bed2binThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus);
		}
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
	if (nIndividuals <= 1){
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( (maCounts.size() % nIndividuals) > 0 ){
		logMessages_ += "ERROR: length of allele count vector (" + std::to_string( maCounts.size() ) + " is not divisible by the provided number of individuals (" +
			std::to_string(nIndividuals) + "\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(nIndividuals) + std::string( ") in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( maCounts.empty() ){
		logMessages_ += "ERROR: empty vector of minor allele counts\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ) {
		nThreads_ = std::thread::hardware_concurrency();
	}
	logMessages_ += "Number of individuals: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: " + std::to_string(nLoci_) + "\n";
	logMessages_ += "Number of threads: " + std::to_string(nThreads_) + "\n";
	binLocusSize_ = nIndividuals_ / byteSize_ + static_cast<size_t>( (nIndividuals_ % byteSize_) > 0 );
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);
	const size_t ranVecSize     = nIndividuals_ / llWordSize_ + static_cast<size_t>( (nIndividuals_ % llWordSize_) > 0 );
	const size_t nLociPerThread = nLoci_ / nThreads_;
	if (nLociPerThread > 0){
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLociPerThread)};
		threadRanges.back().second = nLoci_;
		std::vector< std::future<void> > tasks;
		tasks.reserve(nThreads_);
		for (const auto &eachTR : threadRanges){
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
	if (this != &toMove){
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
	assert( ( binGenotypes_.size() < std::streamsize::max() ) && "ERROR: binGenotypes_ size exceeds maximum streamsize in GenoTableBin.saveGenoBinary()");
	out.open(outFileName, std::ios::out | std::ios::binary | std::ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), static_cast<std::streamsize>( binGenotypes_.size() ) ); // OK because we are casting to const char*
	out.close();
}

void GenoTableBin::allJaccardLD(const std::string &ldFileName) const {
	if (nLoci_ > maxNlocusPairs_){
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxNlocusPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * sizeof(IndexedPairLD) );      // use half to leave resources for other operations
	const size_t nPairs   = nLoci_ * (nLoci_ - 1) / 2;
	logMessages_         += "Calculating all pairwise LD using full Jaccard similarity estimates\n";
	logMessages_         += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "; ";
	if (nPairs > maxInRAM){
		// too many loci to fit the LD matrix in RAM
		// will work on chunks and save as we go
		logMessages_ += "calculating in chunks\n";
		const size_t nChunks        = nPairs / maxInRAM;
		const size_t remainingPairs = nPairs % maxInRAM;
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		output << "locus1\tlocus2\tjaccard\trSq\n";
		const size_t nLocusPairsPerThread = maxInRAM / nThreads_;
		size_t overallPairInd = 0;
		if (nLocusPairsPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLocusPairsPerThread)};
			const size_t excessLoci    = maxInRAM - threadRanges.back().second;
			threadRanges.back().second = maxInRAM;
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<IndexedPairLD> LDmatChunk(maxInRAM);
				overallPairInd  = jaccardThreaded_(threadRanges, overallPairInd, LDmatChunk);
				overallPairInd += excessLoci;
				saveValues(LDmatChunk, output);
			}
		} else {
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<IndexedPairLD> LDmatChunk(maxInRAM);
				const std::pair<size_t, size_t> ramRange{0, maxInRAM};
				jaccardBlock_(ramRange, overallPairInd, LDmatChunk);
				overallPairInd += maxInRAM;
				saveValues(LDmatChunk, output);
			}
		}
		overallPairInd = maxInRAM * nChunks;
		if (remainingPairs > 0){
			std::vector<IndexedPairLD> LDmatChunk(remainingPairs);
			const size_t nRemainPairsPerThread = remainingPairs / nThreads_;
			if (nRemainPairsPerThread > 0){
				std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nRemainPairsPerThread)};
				threadRanges.back().second = remainingPairs;
				overallPairInd = jaccardThreaded_(threadRanges, overallPairInd, LDmatChunk);
			} else {
				const std::pair<size_t, size_t> remainingRange{0, remainingPairs};
				jaccardBlock_(remainingRange, overallPairInd, LDmatChunk);
			}
			saveValues(LDmatChunk, output);
		}
		output.close();
	} else {
		logMessages_ += "calculating in one go\n";
		std::vector<IndexedPairLD> LDmat(nPairs);
		const size_t nLocusPairsPerThread = LDmat.size() / nThreads_;
		if (nLocusPairsPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLocusPairsPerThread)};
			threadRanges.back().second = LDmat.size();
			jaccardThreaded_(threadRanges, 0, LDmat);
		} else {
			const std::pair<size_t, size_t> fullRange{0, LDmat.size()};
			jaccardBlock_(fullRange, 0, LDmat);
		}
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		output << "locus1\tlocus2\tjaccard\trSq\n";
		saveValues(LDmat, output);
		output.close();
	}
}

void GenoTableBin::saveLogFile() const {
	std::fstream outLog;
	outLog.open(logFileName_, std::ios::out | std::ios::trunc);
	outLog << logMessages_;
	outLog.close();
}

void GenoTableBin::bed2binBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const size_t &firstLocusInd, const size_t &bedLocusLength) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint64_t locSeed{0};
	{
		std::lock_guard<std::mutex> lockGuard(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
	size_t begByte{firstLocusInd * binLocusSize_};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus){
		const size_t begInd{iBedLocus * bedLocusLength};
		binarizeBedLocus(begInd, bedLocusLength, bedData, nIndividuals_, locPRNG, begByte, binLocusSize_, binGenotypes_);
		begByte += binLocusSize_;
	}
}

size_t GenoTableBin::bed2binThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const size_t &firstLocusInd, const size_t &bedLocusLength){
	size_t locusInd = firstLocusInd;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges){
		tasks.emplace_back(
			std::async([this, &bedData, &eachTR, locusInd, bedLocusLength]{
				bed2binBlk_(bedData, eachTR, locusInd, bedLocusLength);
			})
		);
		locusInd += eachTR.second - eachTR.first;
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
	for (size_t iLocus = locusIndRange.first; iLocus < locusIndRange.second; ++iLocus){
		// Fill the random byte vector
		for (auto &randValue : rand){
			randValue = locPRNG.ranInt();
		}
		size_t iIndiv         = 0;
		size_t iRB            = 0;                                                                 // randBytes index
		const size_t begIndiv = iLocus * nIndividuals_;
		const size_t begByte  = iLocus * binLocusSize_;
		std::vector<uint8_t> missMasks(binLocusSize_, 0);
		size_t i0Byte = 0;                                                                         // to index the missMasks vector
		for (size_t iByte = begByte; iByte < begByte + binLocusSize_ - 1; ++iByte){                // treat the last byte separately
			uint8_t binByte = 0;
			for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte){
				auto curIndiv            = static_cast<uint8_t>(macData[begIndiv + iIndiv]);       // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= middleMask;                                            // mask everything in the middle
				const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
				missMasks[i0Byte]        |= (missingMask << iInByte);
				curIndiv                 &= endTwoBitMask;
				const uint8_t randMask    = (randBytes[iRB] >> iInByte) & oneBit_;                 // 0b00000000 or 0b00000001 with equal chance
				uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
				curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
				binByte                  |= curBitMask << iInByte;
				++iIndiv;
			}
			// should be safe: each thread accesses different vector elements
			binGenotypes_[iByte] = binByte;
			++i0Byte;
			++iRB;
		}
		// now deal with the last byte in the individual
		uint8_t lastBinByte = 0;
		for (uint8_t iRem = 0; iRem < remainderInd; ++iRem){
			auto curIndiv             = static_cast<uint8_t>(macData[begIndiv + iIndiv]);          // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= middleMask;                                                // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
			missMasks.back()         |= (missingMask << iRem);
			curIndiv                 &= endTwoBitMask;                                                
			const uint8_t randMask    = (randBytes[binLocusSize_ - 1] >> iRem) & oneBit_;          // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                   // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                              // zero it out if missing value is set

			lastBinByte              |= curBitMask << iRem;
			++iIndiv;
		}
		// should be safe: each thread accesses different vector elements
		binGenotypes_[begByte + binLocusSize_ - 1] = lastBinByte;
		float maf = static_cast<float>( countSetBits(binGenotypes_, begByte, binLocusSize_) ) / static_cast<float>(nIndividuals_);
		if (maf > 0.5){ // always want the alternative to be the minor allele
			i0Byte = 0;
			for (size_t i = begByte; i < begByte + binLocusSize_; ++i){
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
	for (size_t iVecInd = blockVecRange.first; iVecInd < blockVecRange.second; ++iVecInd){
		// compute row and column indexes from the vectorized lower triangle index
		// got these expressions by combining various web sources and verifying
		const size_t kpIdx  = nnLoci - curJacMatInd;
		const size_t pIdx   = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kpIdx) ) ) - 1) / 2;
		const size_t row    = nLoci_ - 2 - pIdx;
		const size_t col    = nLoci_ - (kpIdx - pIdx * (pIdx + 1) / 2) - 1;
		const size_t rowBin = row * binLocusSize_;
		const size_t colBin = col * binLocusSize_;
		for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc){
			locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] & binGenotypes_[colBin + iBinLoc];
		}
		const uint64_t isect = countSetBits(locus);
		for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc){
			locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] | binGenotypes_[colBin + iBinLoc];
		}
		const uint64_t uni     = countSetBits(locus);
		const uint64_t locus1n = countSetBits(binGenotypes_, rowBin, binLocusSize_);
		const uint64_t locus2n = countSetBits(binGenotypes_, colBin, binLocusSize_);
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
		ldVec[iVecInd].rSq     = dValue * dValue / ( pApB * (1.0F - locus1p) * (1.0F - locus2p) );  // nan if either of the loci monomorphic; could be if all hets are called as '0'
		++curJacMatInd;
	}
}

size_t GenoTableBin::jaccardThreaded_(const std::vector< std::pair<size_t, size_t> > &pairIndRanges, const size_t &blockStartAll, std::vector<IndexedPairLD> &ldVec) const {
	size_t overallPairInd = blockStartAll;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : pairIndRanges){
		tasks.emplace_back(
			std::async([this, &ldVec, &eachTR, overallPairInd]{
				jaccardBlock_(eachTR, overallPairInd, ldVec);
			})
		);
		overallPairInd += eachTR.second - eachTR.first;
	}
	for (const auto &eachThread : tasks){
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
GenoTableHash::GenoTableHash(const std::string &inputFileName, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, std::string logFileName)
								: kSketches_{kSketches}, fSketches_{static_cast<float>(kSketches)}, nLoci_{0}, nThreads_{nThreads}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime{std::time(nullptr)};
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype hashing from a .bed file started on " + logStream.str() + "\n";
	logStream.clear();
	if (nIndividuals <= 1){
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ < 3){
		logMessages_ += "ERROR: number of sketches (" + std::to_string(kSketches_) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: sketch number must be at least three in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	// Round up the number of individuals to nearest divisible by kSketches_
	sketchSize_   = nIndividuals / kSketches_ + static_cast<size_t>( (nIndividuals % kSketches_) > 0 );
	nIndividuals_ = sketchSize_ * kSketches_;
	if (sketchSize_ >= emptyBinToken_){
		logMessages_ += "ERROR: sketch size (" + std::to_string(sketchSize_) + ") is too big; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: Number of sketches (") + std::to_string(kSketches_) + std::string(") implies sketch size (") +
			std::to_string(sketchSize_) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t nBedBytes = nIndividuals / bedGenoPerByte_ + static_cast<size_t>( (nIndividuals % bedGenoPerByte_) > 0 );
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ){
		nThreads_ = std::thread::hardware_concurrency();
	}
	logMessages_ += "Number of threads used: " + std::to_string(nThreads_) + "\n";
	std::fstream inStr;
	// Start by measuring file size
	inStr.open(inputFileName, std::ios::in | std::ios::binary | std::ios::ate);
	if ( inStr.fail() ){
		logMessages_ += "ERROR: failed to open file " + inputFileName + "; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const std::streamoff endPosition = inStr.tellg();
	if ( endPosition < static_cast<std::streamoff>(nMagicBytes_) ){
		logMessages_ += "ERROR: no loci in the input .bed file " + inputFileName + "; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	inStr.close();
	const size_t fileSize  = static_cast<uint64_t>(endPosition) - nMagicBytes_;
	nLoci_          = fileSize / nBedBytes;
	logMessages_   += "Number of individuals: " + std::to_string(nIndividuals) + "\n";
	logMessages_   += "Number of individuals to hash: " + std::to_string(nIndividuals_) + "\n";
	logMessages_   += "Number of loci: " + std::to_string(nLoci_) + "\n";
	logMessages_   += "Hash size: " + std::to_string(kSketches_) + "\n";
	locusSize_      = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;                    // round up to the nearest multiple of 8
	nFullWordBytes_ = (nIndividuals_ - 1) / byteSize_;
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStr.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	const size_t nBedBytesPerLocus = nIndividuals / bedGenoPerByte_ + static_cast<size_t>(nIndividuals % bedGenoPerByte_ > 0);
	const size_t ramSize           = getAvailableRAM() / 2UL;                               // measuring here, after all the major allocations; use half to leave resources for other operations
	size_t nBedLociToRead          = ramSize / nBedBytesPerLocus;                           // number of .bed loci to read at a time
	nBedLociToRead                 = (nBedLociToRead < nLoci_ ? nBedLociToRead : nLoci_);
	const size_t remainingLoci     = nLoci_ % nBedLociToRead;
	const size_t remainingBytes    = remainingLoci * nBedBytesPerLocus;
	const size_t nChunks           = nLoci_ / nBedLociToRead;
	size_t nBedBytesToRead         = nBedLociToRead * nBedBytesPerLocus;
	nBedBytesToRead                = (nBedBytesToRead > std::numeric_limits<std::streamsize>::max() ? std::numeric_limits<std::streamsize>::max() : nBedBytesToRead);
	size_t nLociPerThread          = nBedLociToRead / nThreads_;
	std::vector<char> bedChunkToRead(nBedBytesToRead, 0);
	logMessages_ += "RAM available for reading the .bed file: " + std::to_string(ramSize) + " bytes\n";
	logMessages_ += ".bed file will be read in " + std::to_string(nChunks) + " chunk(s)\n";
	// Sample with replacement additional individuals to pad out the total
	std::vector< std::pair<size_t, size_t> > addIndv;
	for (size_t iAddIndiv = nIndividuals; iAddIndiv < nIndividuals_; ++iAddIndiv){
		addIndv.emplace_back( iAddIndiv, rng_.sampleInt(nIndividuals) );
	}
	if ( !addIndv.empty() ){
		std::string addIndexes;
		for (const auto &eachIdx : addIndv){
			addIndexes += std::to_string(eachIdx.second) + " ";
		}
		logMessages_ += "Re-sampled individuals: " + addIndexes + "\n";
	}
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{rng_.fyIndexesUp(nIndividuals_)};
	std::vector<uint32_t> seeds{static_cast<uint32_t>( rng_.ranInt() )};

	size_t locusInd = 0;
	if (nLociPerThread > 0){
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLociPerThread)};
		const size_t excessLoci = nBedLociToRead - threadRanges.back().second;
		threadRanges.back().second = nBedLociToRead;
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			assert( ( nBedBytesToRead < std::numeric_limits<std::streamsize>::max() )
					&& "ERROR: amount to read larger than maximum streamsize in GenoTableHash constructor");
			inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(nBedBytesToRead) );
			locusInd = bed2ophThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus, ranInts, addIndv, seeds);
			locusInd += excessLoci;
		}
	} else {
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			assert( ( nBedBytesToRead < std::numeric_limits<std::streamsize>::max() )
					&& "ERROR: amount to read larger than maximum streamsize in GenoTableHash constructor");
			inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(nBedBytesToRead) );
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, 1)};
			locusInd = bed2ophThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus, ranInts, addIndv, seeds);
		}
	}
	if (remainingLoci > 0) {
		bedChunkToRead.resize(remainingBytes);
		assert( ( remainingBytes < std::numeric_limits<std::streamsize>::max() )
				&& "ERROR: amount left to read larger than maximum streamsize in GenoTableHash constructor");
		inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(remainingBytes) );
		nLociPerThread = remainingLoci / nThreads_;
		if (nLociPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLociPerThread)};
			const size_t excessLoci = remainingLoci - threadRanges.back().second;
			threadRanges.back().second = remainingLoci;
			locusInd = bed2ophThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus, ranInts, addIndv, seeds);
			locusInd += excessLoci;
		} else {
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, 1)};
			locusInd = bed2ophThreaded_(bedChunkToRead, threadRanges, locusInd, nBedBytesPerLocus, ranInts, addIndv, seeds);
		}
	}
	inStr.close();
}

GenoTableHash::GenoTableHash(const std::vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, std::string logFileName) 
		: nIndividuals_{nIndividuals}, kSketches_{kSketches}, fSketches_{static_cast<float>(kSketches)}, nLoci_{maCounts.size() / nIndividuals}, nThreads_{nThreads}, logFileName_{std::move(logFileName)} {
	std::stringstream logStream;
	const time_t startTime = std::time(nullptr);
	struct tm buf{};
	logStream << std::put_time(localtime_r(&startTime, &buf), "%b %e %Y %H:%M %Z");
	logMessages_ = "Genotype hashing from a minor allele count vector started on " + logStream.str() + "\n";
	logStream.clear();
	if (nIndividuals <= 1){
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( (maCounts.size() % nIndividuals) > 0){
		logMessages_ += "ERROR: minor allele vector size (" + std::to_string( maCounts.size() ) + ") is not evenly divisible by the number of individuals (" + std::to_string(nIndividuals_) + "); aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(nIndividuals) + std::string( ") in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( maCounts.empty() ){
		logMessages_ += "ERROR: minor allele count vector is empty; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ < 3){
		logMessages_ += "ERROR: sketch size (" + std::to_string(kSketches_) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: sketch size must be at least three in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ){
		nThreads_ = std::thread::hardware_concurrency();
	}
	sketchSize_   = nIndividuals / kSketches_ + static_cast<size_t>( (nIndividuals % kSketches_) > 0 );
	nIndividuals_ = sketchSize_ * kSketches_;
	if (sketchSize_ >= emptyBinToken_){
		logMessages_ += "ERROR: sketch size (" + std::to_string(sketchSize_) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: Number of sketches (") + std::to_string(kSketches_) + std::string(") implies sketch size (") +
			std::to_string(sketchSize_) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
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
	logMessages_ += "Number of individuals: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: " + std::to_string(nLoci_) + "\n";
	logMessages_ += "Hash size: " + std::to_string(kSketches_) + "\n";

	const size_t nLociPerThread = nLoci_ / nThreads_;
	if (nLociPerThread > 0){
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLociPerThread)};
		threadRanges.back().second = nLoci_;
		std::vector< std::future<void> > tasks;
		tasks.reserve(nThreads_);
		for (const auto &eachTR : threadRanges){
			tasks.emplace_back(
				std::async([this, &maCounts, &eachTR, ranVecSize, &ranInts, &seeds]{
					mac2ophBlk_(maCounts, eachTR.first, eachTR.second, ranVecSize, ranInts, seeds);
				})
			);
		}
		for (const auto &eachTask : tasks){
			eachTask.wait();
		}
	} else {
		mac2ophBlk_(maCounts, 0, nLoci_, ranVecSize, ranInts, seeds);
	}
}

GenoTableHash::GenoTableHash(GenoTableHash &&toMove) noexcept : nIndividuals_{0}, kSketches_{0}, fSketches_{0.0}, sketchSize_{0}, nLoci_{0}, locusSize_{0}, nFullWordBytes_{0}, nThreads_{0} {
	*this = std::move(toMove);
}

GenoTableHash& GenoTableHash::operator=(GenoTableHash &&toMove) noexcept {
	if (this != &toMove){
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

void GenoTableHash::allHashLD(const std::string &ldFileName) const {
	if (nLoci_ > maxPairs_){
		logMessages_ += "ERROR: too many loci (" + std::to_string(nLoci_) + ") to do all pairwise LD; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * sizeof(IndexedPairSimilarity) );      // use half to leave resources for other operations
	const size_t nPairs   = nLoci_ * (nLoci_ - 1) / 2;
	logMessages_         += "Calculating all pairwise LD\n";
	logMessages_         += "Maximum number of locus pairs that fit in RAM: " + std::to_string(maxInRAM) + "; ";
	if (nPairs > maxInRAM){
		logMessages_ += "calculating in chunks\n";
		// too many loci to fit the LD matrix in RAM
		// will work on chunks and save as we go
		const size_t nChunks        = nPairs / maxInRAM;
		const size_t remainingPairs = nPairs % maxInRAM;
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		output << "groupID\tlocus1\tlocus2\tjaccLD\n";
		const size_t nLocusPairsPerThread = maxInRAM / nThreads_;
		size_t overallPairInd = 0;
		if (nLocusPairsPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLocusPairsPerThread)};
			const size_t excessLoci    = maxInRAM - threadRanges.back().second;
			threadRanges.back().second = maxInRAM;
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<IndexedPairSimilarity> LDmatChunk(maxInRAM);
				overallPairInd  = hashJacThreaded_(threadRanges, overallPairInd, LDmatChunk);
				overallPairInd += excessLoci;
				saveValues(LDmatChunk, output);
			}
		} else {
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<IndexedPairSimilarity> LDmatChunk(maxInRAM);
				const std::pair<size_t, size_t> chunkRange(0, maxInRAM);
				hashJacBlock_(chunkRange, overallPairInd, LDmatChunk);
				overallPairInd += maxInRAM;
				saveValues(LDmatChunk, output);
			}
		}
		overallPairInd = maxInRAM * nChunks;
		if (remainingPairs > 0){
			std::vector<IndexedPairSimilarity> LDmatChunk(remainingPairs);
			const size_t nRemainPairsPerThread = remainingPairs / nThreads_;
			if (nRemainPairsPerThread > 0){
				std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nRemainPairsPerThread)};
				threadRanges.back().second = remainingPairs;
				overallPairInd             = hashJacThreaded_(threadRanges, overallPairInd, LDmatChunk);
			} else {
				const std::pair<size_t, size_t> remainRange(0, remainingPairs);
				hashJacBlock_(remainRange, overallPairInd, LDmatChunk);
			}
			saveValues(LDmatChunk, output);
		}
		output.close();
	} else {
		logMessages_ += "calculating in one go\n";
		std::vector<IndexedPairSimilarity> LDmat(nPairs);
		const size_t nLocusPairsPerThread = LDmat.size() / nThreads_;
		if (nLocusPairsPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nLocusPairsPerThread)};
			threadRanges.back().second = LDmat.size();
			hashJacThreaded_(threadRanges, 0, LDmat);
		} else {
			std::pair<size_t, size_t> wholeRange( 0, LDmat.size() );
			hashJacBlock_(wholeRange, 0, LDmat);
		}
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		output << "groupID\tlocus1\tlocus2\tjaccLD\n";
		saveValues(LDmat, output);
		output.close();
	}
}

std::unordered_map< uint32_t, std::vector<size_t> > GenoTableHash::makeLDgroups(const size_t &nRowsPerBand) const {
	assert( (nRowsPerBand != 0) && "ERROR: nRowsPerBand must not be 0 in makeLDgroups()" );
	assert( (nRowsPerBand < kSketches_) && "ERROR: nRowsPerBand must be less than kSketches_ in makeLDgroups()" );
	const size_t nBands = kSketches_ / nRowsPerBand;                                                               // only using full-size bands because smaller ones permit inclusion of low-similarity pairs
	assert( ( nBands >= std::numeric_limits<uint16_t>::max() ) &&
			"ERROR: number of bands cannot exceed uint16_t max in makeLDgroups()" );

	logMessages_ += "Grouping loci\n";
	logMessages_ += "Number of rows per band: " + std::to_string(nRowsPerBand) + "\n";
	logMessages_ += "Number of bands: " + std::to_string(nBands) + "\n";

	const auto sketchSeed = static_cast<uint32_t>( rng_.ranInt() );
	std::unordered_map< uint32_t, std::vector<size_t> > ldGroup;                                                   // the hash table

	for (size_t iLocus = 0; iLocus < nLoci_; ++iLocus){
		size_t iSketch = 0;
		for (size_t iBand = 0; iBand < nBands; ++iBand){
			std::vector<uint16_t> bandVec{static_cast<uint16_t>(iBand)};                                           // add the band index to the hash, so that only corresponding bands are compared
			const size_t firstSketchIdx = iSketch + iLocus * kSketches_;                                           // iSketch tracks band IDs
			for (size_t iInBand = firstSketchIdx; iInBand < firstSketchIdx + nRowsPerBand; ++iInBand){
				bandVec.push_back(sketches_[iInBand]);
			}
			const uint32_t hash = murMurHash(0, bandVec.size(), bandVec, sketchSeed);
			ldGroup[hash].push_back(iLocus);
			iSketch += nRowsPerBand;
		}
	}
	// remove groups with one locus
	auto ldgIt = ldGroup.begin();
	while ( ldgIt != ldGroup.end() ) {
		if (ldgIt->second.size() < 2) {
			ldgIt = ldGroup.erase(ldgIt);
		} else {
			++ldgIt;
		}
	}
	// Then maybe iterate through bands inside a locus, keep searching for the hash until found, if nothing found insert under the first band
	return ldGroup;
}

void GenoTableHash::makeLDgroups(const size_t &nRowsPerBand, const std::string &outFileName) const {
	std::unordered_map< uint32_t, std::vector<size_t> > ldGroups{this->makeLDgroups(nRowsPerBand)};
	logMessages_ += "Saving group IDs only\n";
	std::fstream out;
	out.open(outFileName, std::ios::out | std::ios::trunc);
	out << "groupID\tlocusIdx\n";
	for (const auto &eachGroup : ldGroups){
		for (const auto &locusIdx : eachGroup.second){
			out << "G" << eachGroup.first + 1 << "\t" << locusIdx + 1 << "\n";
		}
	}
	out.close();
}

void GenoTableHash::ldInGroups(const size_t &nRowsPerBand, const std::string &outFileName) const {
	std::unordered_map< uint32_t, std::vector<size_t> > ldGroups{this->makeLDgroups(nRowsPerBand)};
	
	// there is a possibility that all retained pairs will not fit in RAM
	// dealing with this is non-trivial (need the whole thing in RAM to eliminate duplicate pairs), so for now will let the system deal with memory issues
	logMessages_ += "Estimating LD in groups\n";
	
	// identify index pairs that are retained in groups
	// some will be repeatedly identified by makeLDgroups(), so we will need to eliminate duplicates
	std::vector<IndexedPairSimilarity> hashJacGroups;
	uint32_t groupInd = 0;
	for (auto &ldg : ldGroups){
		for (size_t iLocus = 0; iLocus < ldg.second.size() - 1; ++iLocus){
			for (size_t jLocus = iLocus + 1; jLocus < ldg.second.size(); ++jLocus){
				hashJacGroups.emplace_back(IndexedPairSimilarity{0.0, ldg.second[iLocus], ldg.second[jLocus], groupInd});
			}
		}
		++groupInd;
		ldg.second.clear();
	}
	logMessages_ += "Number of locus pairs before removing duplicates: " + std::to_string( hashJacGroups.size() ) + "\n";
	std::sort(hashJacGroups.begin(), hashJacGroups.end(),
				[](const IndexedPairSimilarity &first, const IndexedPairSimilarity &second){
					return (first.element1ind == second.element1ind ? first.element2ind < second.element2ind : first.element1ind < second.element1ind);
				}
			);
	auto lastUniqueIt = std::unique(hashJacGroups.begin(), hashJacGroups.end(),
				[](const IndexedPairSimilarity &first, const IndexedPairSimilarity &second){
					return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
				}
			);
	hashJacGroups.erase( lastUniqueIt, hashJacGroups.end() );
	hashJacGroups.shrink_to_fit();
	logMessages_ += "Number of locus pairs after removing duplicates: " + std::to_string( hashJacGroups.size() ) + "\n";
	// estimate locus pair similarities
	const size_t nPairsPerThread = hashJacGroups.size() / nThreads_;
	if (nPairsPerThread > 0){
		std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(nThreads_, nPairsPerThread)};
		threadRanges.back().second = hashJacGroups.size(); // assuming the number of threads << number of groups, this should not lead to noticeable thread imbalance
		std::vector< std::future<void> > tasks;
		tasks.reserve(nThreads_);
		for (const auto &eachTR : threadRanges){
			tasks.emplace_back(
				std::async([this, &eachTR, &hashJacGroups]{
					hashJacBlock_(eachTR.first, eachTR.second, hashJacGroups);
				})
			);
		}
		for (const auto &eachTask : tasks){
			eachTask.wait();
		}
	} else {
		std::vector< std::future<void> > tasks;
		tasks.reserve( hashJacGroups.size() );
		for (size_t iPair = 0; iPair < hashJacGroups.size(); ++iPair){
			tasks.emplace_back(
				std::async([this, iPair, &hashJacGroups]{
					hashJacBlock_(iPair, iPair + 1, hashJacGroups);
				})
			);
		}
		for (const auto &eachTask : tasks){
			eachTask.wait();
		}
	}

	std::fstream out;
	out.open(outFileName, std::ios::out | std::ios::trunc);
	out << "groupID\tlocus1\tlocus2\tjaccLD\n";
	saveValues(hashJacGroups, out);
	out.close();
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
	while(iByte < nFullWordBytes_){
		for (uint8_t iInLocusByte = 0; iInLocusByte < byteSize_; ++iInLocusByte){
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
			const uint16_t temp1     = ( bytePair ^ (bytePair >> shiftDistance) ) & mask;
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
	while (iIndiv < nIndividuals_ - 1){
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
		const uint16_t temp1     = ( bytePair ^ (bytePair >> shiftDistance) ) & mask;
		const auto temp2         = static_cast<uint16_t>(temp1 << shiftDistance);
		bytePair                ^= temp1 ^ temp2;
		// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus)
		// Must modify the current byte in each loop iteration because permutation indexes may fall into it
		binLocus[iByte]       ^= static_cast<uint8_t>( ( binLocus[iByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask) );
		binLocus[permByteInd] ^= static_cast<uint8_t>( ( binLocus[permByteInd] ^ static_cast<uint8_t>(bytePair >> byteSize_) ) & perMask );
		++iInLocusByte;
	}
}

void GenoTableHash::locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds, std::vector<uint8_t> &binLocus){
	// Start with a permutation to make OPH
	permuteBits_(permutation, binLocus);
	// Now make the sketches
	std::vector<size_t> filledIndexes;                                                       // indexes of the non-empty sketches
	size_t iByte{0};
	size_t sketchBeg{locusInd * kSketches_};
	size_t iSketch{0};
	uint64_t sketchTail{0};                                                                  // left over buts from beyond the last full byte of the previous sketch
	size_t locusChunkSize = (llWordSize_ > binLocus.size() ? binLocus.size() : llWordSize_);
	while ( iByte < binLocus.size() ){
		uint64_t nWordUnsetBits{wordSizeInBits_};
		uint64_t nSumUnsetBits{0};
		while ( (nWordUnsetBits == wordSizeInBits_) && ( iByte < binLocus.size() ) ){
			uint64_t locusChunk{allBitsSet_};
			assert( ( iByte < binLocus.size() ) && "ERROR: iByte must be less than locus size in bytes in locusOPH_()" );
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
		if (iSketch >= kSketches_){
			break;
		}
		filledIndexes.push_back(iSketch);
		sketches_[sketchBeg + iSketch] = static_cast<uint16_t>(nSumUnsetBits % sketchSize_);
		++iSketch;
		const uint64_t bitsDone{iSketch * sketchSize_};
		iByte      = bitsDone / byteSize_;
		sketchTail = bitsDone % byteSize_;
	}
	size_t iSeed = 0;                                             // index into the seed vector
	if (filledIndexes.size() == 1){
		for (size_t iSk = 0; iSk < kSketches_; ++iSk){ // this will overwrite the one assigned sketch, but the wasted operation should be swamped by the rest
			// should be safe: each thread accesses different vector elements
			sketches_[sketchBeg + iSk] = sketches_[filledIndexes[0] + sketchBeg];
		}
	} else if (filledIndexes.size() != kSketches_){
		if ( filledIndexes.empty() ){ // in the case where the whole locus is monomorphic, pick a random index as filled
			filledIndexes.push_back( rng_.sampleInt(kSketches_) );
		}
		size_t emptyCount = kSketches_ - filledIndexes.size();
		while (emptyCount > 0){
			for (const auto &eachFI : filledIndexes){
				auto newIdx = static_cast<uint32_t>(murMurHash(eachFI, seeds[iSeed]) % kSketches_ + sketchBeg);
				// should be safe: each thread accesses different vector elements
				if (sketches_[newIdx] == emptyBinToken_){
					sketches_[newIdx] = sketches_[eachFI + sketchBeg];
					--emptyCount;
					break;
				}
			}
			++iSeed;
			std::lock_guard<std::mutex> lock(mtx_);      // lock before measuring to ensure that the size is valid
			if ( iSeed == seeds.size() ){
				seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
			}
		}
	}
}

void GenoTableHash::bed2ophBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t>&bedLocusIndRange, const size_t &firstLocusInd,
					const size_t &bedLocusLength, const std::vector<size_t> &permutation, std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint64_t locSeed{0};
	{
		std::lock_guard<std::mutex> lock(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
	size_t iLocus{firstLocusInd};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus){
		std::vector<uint8_t> binLocus(locusSize_, 0);
		const size_t begInd{iBedLocus * bedLocusLength};
		binarizeBedLocus(begInd, bedLocusLength, bedData, nIndividuals_, locPRNG, 0, locusSize_, binLocus);
		// pad the locus to have an even number of sketches 
		for (const auto &addI : padIndiv){
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
			const uint16_t temp1     = ( bytePair ^ (bytePair >> shiftDistance) ) & mask;
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

size_t GenoTableHash::bed2ophThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const size_t &firstLocusInd,
							const size_t &bedLocusLength, const std::vector<size_t> &permutation, std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds){
	size_t locusInd = firstLocusInd;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges){
		tasks.emplace_back(
			std::async([this, &bedData, &eachTR, locusInd, bedLocusLength, &permutation, &padIndiv, &seeds]{
				bed2ophBlk_(bedData, eachTR, locusInd, bedLocusLength, permutation, padIndiv, seeds);
			})
		);
		locusInd += eachTR.second - eachTR.first;
	}
	for (const auto &eachTask : tasks){
		eachTask.wait();
	}
	return locusInd;
}

void GenoTableHash::mac2ophBlk_(const std::vector<int> &macData, const size_t &startLocusInd, const size_t &endLocusInd, const size_t &randVecLen, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds){
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
	auto *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	for (size_t iLocus = startLocusInd; iLocus < endLocusInd; ++iLocus){
		// Fill the random byte vector
		for (auto &randValue : rand){
			randValue = locPRNG.ranInt();
		}
		size_t iIndiv         = 0;
		const size_t begIndiv = iLocus * nIndividuals_;
		std::vector<uint8_t> missMasks(locusSize_, 0);
		std::vector<uint8_t> binLocus(locusSize_, 0);
		size_t i0Byte = 0;                                                                         // to index the missMasks vector
		for (size_t iByte = 0; iByte < locusSize_ - 1; ++iByte){                                   // treat the last byte separately
			for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte){
				auto curIndiv             = static_cast<uint8_t>(macData[begIndiv + iIndiv]);      // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= middleMask;                                            // mask everything in the middle
				const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
				missMasks[i0Byte]        |= (missingMask << iInByte);
				curIndiv                 &= endTwoBitMask;
				const uint8_t randMask    = (randBytes[iByte] >> iInByte) & oneBit_;               // 0b00000000 or 0b00000001 with equal chance
				uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
				curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
				binLocus[iByte]          |= curBitMask << iInByte;
				++iIndiv;
			}
			++i0Byte;
		}
		// now deal with the last byte in the individual
		for (uint8_t iRem = 0; iRem < remainderInd; ++iRem){
			auto curIndiv             = static_cast<uint8_t>(macData[begIndiv + iIndiv]);          // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= middleMask;                                                // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
			missMasks.back()         |= (missingMask << iRem);
			curIndiv                 &= endTwoBitMask;                                                
			const uint8_t randMask    = (randBytes[locusSize_ - 1] >> iRem) & oneBit_;             // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                   // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                              // zero it out if missing value is set

			binLocus.back()          |= curBitMask << iRem;
			++iIndiv;
		}
		float maf = static_cast<float>( countSetBits(binLocus) ) / static_cast<float>(nIndividuals_);
		if (maf > 0.5){ // always want the alternative to be the minor allele
			for (size_t iBL = 0; iBL < locusSize_; ++iBL){
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
	for (size_t iVecInd = blockRange.first; iVecInd < blockRange.second; ++iVecInd){
		// compute row and column indexes from the vectorized by column lower triangle index
		// got these expressions by combining various web sources and verifying
		const size_t kpIdx = nnLoci - curJacMatInd;
		const size_t pIdx  = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kpIdx) ) ) - 1) / 2;
		const size_t row   = nLoci_ - 2 - pIdx;
		const size_t col   = nLoci_ - (kpIdx - pIdx * (pIdx + 1) / 2) - 1;
		const auto rowSk   = static_cast<std::vector<uint16_t>::difference_type>(row * kSketches_);
		const auto colSk   = static_cast<std::vector<uint16_t>::difference_type>(col * kSketches_);
		auto start         = sketches_.begin() + rowSk;
		// count equal elements using the inner_product idiom
		int simVal = std::inner_product( start, start + static_cast<std::vector<uint16_t>::difference_type>(kSketches_),
				sketches_.begin() + colSk, 0, std::plus<>(), std::equal_to<>() );
		// should be safe: each thread accesses different vector elements
		hashJacVec[iVecInd].groupID         = 0;
		hashJacVec[iVecInd].element1ind     = row;
		hashJacVec[iVecInd].element2ind     = col;
		hashJacVec[iVecInd].similarityValue = static_cast<float>(simVal) / fSketches_;
		++curJacMatInd;
	}
}

void GenoTableHash::hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, std::vector<IndexedPairSimilarity> &indexedJacc) const {
	for (size_t iBlock = blockStartVec; iBlock < blockEndVec; ++iBlock){
		const auto start1 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(indexedJacc[iBlock].element1ind * kSketches_);
		const auto start2 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(indexedJacc[iBlock].element2ind * kSketches_);
		// count equal elements using the inner_product idiom
		int simVal = std::inner_product( start1, start1 + static_cast<std::vector<uint16_t>::difference_type>(kSketches_),
				start2, 0, std::plus<>(), std::equal_to<>() );
		// should be safe: each thread accesses different vector elements
		indexedJacc[iBlock].similarityValue = static_cast<float>(simVal) / fSketches_;
	}
}

size_t GenoTableHash::hashJacThreaded_(const std::vector< std::pair<size_t, size_t> > &threadRanges, const size_t &blockStartAll, std::vector<IndexedPairSimilarity> &hashJacVec) const {
	size_t overallPairInd = blockStartAll;
	std::vector< std::future<void> > tasks;
	tasks.reserve(nThreads_);
	for (const auto &eachTR : threadRanges){
		tasks.emplace_back(
			std::async([this, &hashJacVec, &eachTR, overallPairInd]{
				hashJacBlock_(eachTR, overallPairInd, hashJacVec);
			})
		);
		overallPairInd += eachTR.second - eachTR.first;
	}
	for (const auto &eachTask : tasks){
		eachTask.wait();
	}
	return overallPairInd;
}



