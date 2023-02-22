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
#include <ctime>
#include <iomanip>
#include <immintrin.h>

#include "gvarHash.hpp"
#include "random.hpp"

using namespace BayesicSpace;
// External functions
uint16_t BayesicSpace::countSetBits(uint16_t inVal) {
	uint16_t totSet = 0;
	for (; inVal > 0; ++totSet) {
		inVal &= inVal - 1;
	}
	return totSet;
}

uint64_t BayesicSpace::countSetBits(const std::vector<uint8_t> &inVec) {
	constexpr size_t wordSize{8};
	constexpr uint64_t roundMask{0xfffffffffffffff8};
	uint32_t totSet{0};
	const size_t nWholeWords{inVec.size() & roundMask};
	size_t iByte{0};
	while (iByte < nWholeWords){
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, wordSize);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
		iByte += wordSize;
	}
	if ( nWholeWords < inVec.size() ){
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, inVec.size() - nWholeWords);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
	}
	return totSet;
}

uint64_t BayesicSpace::countSetBits(const std::vector<uint8_t> &inVec, const size_t &start, const size_t &length) {
	constexpr size_t wordSize{8};
	constexpr uint64_t roundMask{0xfffffffffffffff8};
	uint32_t totSet{0};
	const size_t roundLength{length & roundMask};
	const size_t nWholeWords{start + roundLength};
	size_t iByte{start};
	while (iByte < nWholeWords){
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, wordSize);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
		iByte += wordSize;
	}
	if (roundLength < length){
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, length - roundLength);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
	}
	return totSet;
}

size_t BayesicSpace::getAvailableRAM() {
	constexpr size_t defaultSize{2147483648};
	if ( std::ifstream("/proc/meminfo").good() ){
		constexpr size_t bytesInKb{1024};
		constexpr size_t memTokenLength{13};
		std::string memLine;
		std::fstream memInfoStream;
		memInfoStream.open("/proc/meminfo", std::ios::in);
		while ( getline(memInfoStream, memLine) ){
			if (memLine.compare(0, memTokenLength, "MemAvailable:") == 0) {
				break;
			}
		}
		memInfoStream.close();
		std::stringstream memLineStream(memLine);
		std::string freeMemStr;
		memLineStream >> freeMemStr;
		memLineStream >> freeMemStr;
		return static_cast<size_t>( stoi(freeMemStr) ) * bytesInKb; // memory is in kB in the file
	}
	return defaultSize;
}

void testBedMagicBytes(std::array<char, 3> &bytesToTest) {
	constexpr std::array<char, 3> magicBytes{0x6c, 0x1b, 0x01};   // Leading bytes for .bed files
	if (bytesToTest[0] != magicBytes[0]){
		throw std::string("ERROR: first magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (bytesToTest[1] != magicBytes[1]){
		throw std::string("ERROR: second magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (bytesToTest[2] != magicBytes[2]){
		throw std::string("ERROR: third magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
}

std::vector< std::pair<size_t, size_t> > makeThreadRanges(const size_t &nThreads, const size_t &nElementsPerThread){
	std::vector< std::pair<size_t, size_t> > threadRanges;
	size_t bedInd = 0;
	for (size_t iThread = 0; iThread < nThreads; ++iThread){
		threadRanges.emplace_back(bedInd, bedInd + nElementsPerThread);
		bedInd += nElementsPerThread;
	}
	return threadRanges;
}

// GenoTableBin methods
constexpr size_t   GenoTableBin::nMagicBytes_    = 3;                // number of leading bytes for .bed files
constexpr uint8_t  GenoTableBin::oneBit_         = 0b00000001;       // One set bit for masking
constexpr uint8_t  GenoTableBin::byteSize_       = 8;                // Size of one byte in bits
constexpr uint8_t  GenoTableBin::bedGenoPerByte_ = 4;                // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableBin::llWordSize_     = 8;                // 64 bit word size in bytes
constexpr size_t   GenoTableBin::maxNlocusPairs_ = 6074000999UL;     // approximate number of loci that does not overflow with n*(n-1)/2

// Constructors
GenoTableBin::GenoTableBin(const std::string &inputFileName, const size_t &nIndividuals, const size_t &nThreads) : nIndividuals_{nIndividuals}, nThreads_{nThreads} {
	if (nIndividuals <= 1){
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nIndividuals > std::numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
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
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const uint64_t endPosition = static_cast<uint64_t>( inStr.tellg() );
	if (endPosition < nMagicBytes_){
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t nTotalBedBytes{static_cast<uint64_t>(endPosition) - nMagicBytes_};
	inStr.close();
	nLoci_ = nTotalBedBytes / nBedBytesPerLocus;

	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStr.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	binLocusSize_           = nIndividuals_ / byteSize_ + static_cast<size_t>( (nIndividuals_ % byteSize_) > 0 );
	const size_t ranVecSize = nBedBytesPerLocus / llWordSize_ + static_cast<size_t>( (nBedBytesPerLocus % llWordSize_) > 0 );
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
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &eachTR : threadRanges){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, &eachTR, locusInd, nBedBytesPerLocus, ranVecSize]{
						bed2binBlk_(bedChunkToRead, eachTR, locusInd, nBedBytesPerLocus, ranVecSize);
					})
				);
				locusInd += nLociPerThread;
			}
			locusInd += excessLoci;
		}
	} else {
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			assert( ( nBedBytesToRead < std::streamsize::max() ) && "ERROR: nBedBytesToRead exceeds maximum streamsize in GenoTableBin constructor");
			inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(nBedBytesToRead) );
			std::vector< std::future<void> > tasks;
			tasks.reserve(nBedLociToRead);
			for (size_t iBedLocus = 0; iBedLocus < nBedLociToRead; ++iBedLocus){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, iBedLocus, locusInd, nBedBytesPerLocus, ranVecSize]{
						const std::pair<size_t, size_t> oneLocusRange{iBedLocus, iBedLocus + 1};
						bed2binBlk_(bedChunkToRead, oneLocusRange, locusInd, nBedBytesPerLocus, ranVecSize);})
				);
				++locusInd;
			}
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
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &eachTR : threadRanges){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, &eachTR, locusInd, nBedBytesPerLocus, ranVecSize]{
						bed2binBlk_(bedChunkToRead, eachTR, locusInd, nBedBytesPerLocus, ranVecSize);
					})
				);
				locusInd += nLociPerThread;
			}
			locusInd += excessLoci;
		} else {
			std::vector< std::future<void> > tasks;
			tasks.reserve(remainingLoci);
			for (size_t iBedLocus = 0; iBedLocus < remainingLoci; ++iBedLocus){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, iBedLocus, locusInd, nBedBytesPerLocus, ranVecSize]{
						const std::pair<size_t, size_t> oneLocusRange{iBedLocus, iBedLocus + 1};
						bed2binBlk_(bedChunkToRead, oneLocusRange, locusInd, nBedBytesPerLocus, ranVecSize);})
				);
				++locusInd;
			}
		}
	}
	inStr.close();
}

GenoTableBin::GenoTableBin(const std::vector<int> &maCounts, const size_t &nIndividuals, const size_t &nThreads) : nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals}, nThreads_{nThreads} {
	if (nIndividuals <= 1){
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( (maCounts.size() % nIndividuals) > 0 ){
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(nIndividuals) + std::string( ") in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( maCounts.empty() ){
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ) {
		nThreads_ = std::thread::hardware_concurrency();
	}
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
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxNlocusPairs_) + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * sizeof(float) );      // use half to leave resources for other operations
	const size_t nPairs   = nLoci_ * (nLoci_ - 1) / 2;
	if (nPairs > maxInRAM){
		// too many loci to fit the LD matrix in RAM
		// will work on chunks and save as we go
		const size_t nChunks        = nPairs / maxInRAM;
		const size_t remainingPairs = nPairs % maxInRAM;
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		const size_t nLocusPairsPerThread = maxInRAM / nThreads_;
		size_t overallPairInd = 0;
		if (nLocusPairsPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges;
			size_t locusPairInd = 0;
			for (size_t iThread = 0; iThread < nThreads_; ++iThread){
				threadRanges.emplace_back(locusPairInd, locusPairInd + nLocusPairsPerThread);
				locusPairInd += nLocusPairsPerThread;
			}
			const size_t excessLoci    = maxInRAM - threadRanges.back().second;
			threadRanges.back().second = maxInRAM;
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<float> LDmatChunk(maxInRAM, 0.0);
				std::vector< std::future<void> > tasks;
				tasks.reserve(nThreads_);
				for (const auto &eachTR : threadRanges){
					tasks.emplace_back(
						std::async([this, &LDmatChunk, &eachTR, overallPairInd]{
							jaccardBlock_(eachTR, overallPairInd, LDmatChunk);
						})
					);
					overallPairInd += nLocusPairsPerThread;
				}
				overallPairInd += excessLoci;
				for (const auto &eachThread : tasks){
					eachThread.wait();
				}
				for (const auto &ldValue : LDmatChunk){
					output << ldValue << " ";
				}
			}
		} else {
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<float> LDmatChunk(maxInRAM, 0.0);
				const std::pair<size_t, size_t> ramRange{0, maxInRAM};
				jaccardBlock_(ramRange, overallPairInd, LDmatChunk);
				overallPairInd += maxInRAM;
				for (const auto &ldValue : LDmatChunk){
					output << ldValue << " ";
				}
			}
		}
		overallPairInd = maxInRAM * nChunks;
		if (remainingPairs > 0){
			std::vector<float> LDmatChunk(remainingPairs, 0.0);
			const size_t nRemainPairsPerThread = remainingPairs / nThreads_;
			if (nRemainPairsPerThread > 0){
				std::vector< std::pair<size_t, size_t> > threadRanges;
				size_t locusPairInd = 0;
				for (size_t iThread = 0; iThread < nThreads_; ++iThread){
					threadRanges.emplace_back(locusPairInd, locusPairInd + nRemainPairsPerThread);
					locusPairInd += nRemainPairsPerThread;
				}
				threadRanges.back().second = remainingPairs;
				std::vector< std::future<void> > tasks;
				tasks.reserve(nThreads_);
				for (const auto &eachTR : threadRanges){
					tasks.emplace_back(
						std::async([this, &LDmatChunk, &eachTR, overallPairInd]{
							jaccardBlock_(eachTR, overallPairInd, LDmatChunk);
						})
					);
					overallPairInd += nRemainPairsPerThread;
				}
			} else {
				const std::pair<size_t, size_t> remainingRange{0, remainingPairs};
				jaccardBlock_(remainingRange, overallPairInd, LDmatChunk);
			}
			for (const auto &ldValue : LDmatChunk){
				output << ldValue << " ";
			}
		}
		output.close();
	} else {
		std::vector<float> LDmat(nPairs, 0.0);
		const size_t nLocusPairsPerThread = LDmat.size() / nThreads_;
		if (nLocusPairsPerThread > 0){
			std::vector< std::pair<size_t, size_t> > threadRanges;
			size_t locusPairInd = 0;
			for (size_t iThread = 0; iThread < nThreads_; ++iThread){
				threadRanges.emplace_back(locusPairInd, locusPairInd + nLocusPairsPerThread);
				locusPairInd += nLocusPairsPerThread;
			}
			threadRanges.back().second = LDmat.size();
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &eachTR : threadRanges){
				tasks.emplace_back(
					std::async([this, &LDmat, &eachTR]{
						jaccardBlock_(eachTR, eachTR.first, LDmat);
					})
				);
			}
		} else {
			const std::pair<size_t, size_t> fullRange{0, LDmat.size()};
			jaccardBlock_(fullRange, 0, LDmat);
		}
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		for (const auto &ldValue : LDmat){
			output << ldValue << " ";
		}
		output.close();
	}
}

void GenoTableBin::bed2binBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const size_t &firstLocusInd, const size_t &bedLocusLength, const size_t &randVecLen) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint64_t locSeed{0};
	{
		std::lock_guard<std::mutex> lockGuard(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
	const size_t addIndv       = nIndividuals_ % 4UL;
	const size_t addBL         = bedLocusLength - static_cast<size_t>(addIndv > 0);
	size_t begByte             = firstLocusInd * binLocusSize_;
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> static_cast<uint8_t>(binLocusSize_ * byteSize_ - nIndividuals_);
	std::vector<uint64_t> rand(randVecLen);
	auto *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus){
		const size_t begInd      = iBedLocus * bedLocusLength;
		const size_t endWholeBed = begInd + addBL;
		// Fill the random byte vector
		for (auto &randValue : rand){
			randValue = locPRNG.ranInt();
		}
		std::vector<uint8_t> missMasks(binLocusSize_, 0);
		size_t iBinGeno = begByte;                       // binLocus vector index
		size_t iMissMsk = 0;
		size_t iRB      = 0;                             // randBytes index
		uint8_t bedByte = 0;
		// Two bytes of .bed code go into one byte of my binary representation
		// Therefore, work on two consecutive bytes of .bed code in the loop
		for (size_t iBed = begInd; iBed < endWholeBed ; iBed += 2){                // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
			uint8_t binByte{0};
			bedByte = static_cast<uint8_t>(~bedData[iBed]);                                  // flip so that homozygous second allele (usually minor) is set to 11
			uint8_t offsetToBin{0};                                               // move the .bed mask by this much to align with the binarized byte
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks[iMissMsk]      |= curMissMask >> offsetToBin;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask |= randBytes[iRB] & (firstBitMask << 1);
				firstBitMask  &= secondBitMask >> 1;
				binByte       |= firstBitMask >> offsetToBin;
				++offsetToBin;
			}
			const size_t nextIbed{iBed + 1};
			bedByte = static_cast<uint8_t>(~bedData[nextIbed]);
			++iRB;
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks[iMissMsk]      |= curMissMask << offsetToBin;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask |= randBytes[iRB] & (firstBitMask << 1);
				firstBitMask  &= secondBitMask >> 1;
				binByte       |= firstBitMask << offsetToBin; // keep adding to the current binarized byte, so switch the direction of shift
				--offsetToBin;
			}
			// should be safe: each thread accesses different vector elements
			binGenotypes_[iBinGeno] = binByte;
			++iBinGeno;
			++iMissMsk;
			++iRB;
		}
		uint8_t inBedByteOffset = 0;
		uint8_t binByte         = 0;
		if (addIndv > 0) {
			const auto lastBedByte = static_cast<uint8_t>(~bedData[endWholeBed]);
			for (size_t iInd = 0; iInd < addIndv; ++iInd){
				uint8_t firstBitMask    = lastBedByte & (oneBit_ << inBedByteOffset);
				const uint8_t secondBBO = inBedByteOffset + 1;
				uint8_t secondBitMask   = lastBedByte & (oneBit_ << secondBBO);
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks.back()         |= (curMissMask >> inBedByteOffset) << iInd;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask   |= randBytes[iRB] & (firstBitMask << 1);                                       // iRB incremented at the end of the previous loop
				firstBitMask    &= secondBitMask >> 1;
				firstBitMask     = firstBitMask >> inBedByteOffset;
				firstBitMask     = static_cast<uint8_t>(firstBitMask << iInd);
				binByte         |= firstBitMask;
				inBedByteOffset += 2;
				inBedByteOffset  = inBedByteOffset % byteSize_;
			}
			// should be safe: each thread accesses different vector elements
			binGenotypes_[begByte + binLocusSize_ - 1] = binByte;
		}
		float aaCount = static_cast<float>( countSetBits(binGenotypes_, begByte, binLocusSize_) ) / static_cast<float>(nIndividuals_);
		if (aaCount > 0.5){ // always want the alternative to be the minor allele
			for (size_t iBL = 0; iBL < binLocusSize_; ++iBL){
				const size_t biBL = begByte + iBL;
				// should be safe: each thread accesses different vector elements
				binGenotypes_[biBL] = (~binGenotypes_[biBL]) & (~missMasks[iBL]);
			}
			// should be safe: each thread accesses different vector elements
			binGenotypes_[begByte + binLocusSize_ - 1] &= lastByteMask; // unset the remainder bits
		}
		begByte += binLocusSize_;
	}
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
		{
			// should be safe: each thread accesses different vector elements
			binGenotypes_[begByte + binLocusSize_ - 1] = lastBinByte;
		}
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

void GenoTableBin::jaccardBlock_(const std::pair<size_t, size_t> &blockVecRange, const size_t &blockStartAll, std::vector<float> &jaccardVec) const {
	const size_t nnLoci = nLoci_ * (nLoci_ - 1) / 2 - 1; // overflow checked in the calling function; do this here for encapsulation, may move to the calling function later
	std::vector<uint8_t> locus(binLocusSize_);
	size_t curJacMatInd = blockStartAll;
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
		const uint64_t uni = countSetBits(locus);
		for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc){
			locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] | binGenotypes_[colBin + iBinLoc];
		}
		const uint64_t isect = countSetBits(locus);
		// should be safe: each thread accesses different vector elements
		jaccardVec[iVecInd]  = static_cast<float>(uni) / static_cast<float>(isect);
		++curJacMatInd;
	}
}

// GenoTableHash methods
constexpr std::array<char, 3> GenoTableHash::magicBytes_{0x6c, 0x1b, 0x01};               // Leading bytes for .bed files 
constexpr uint8_t  GenoTableHash::oneBit_         = 0b00000001;                           // One set bit for masking 
constexpr uint8_t  GenoTableHash::byteSize_       = 8;                                    // Size of one byte in bits 
constexpr uint8_t  GenoTableHash::bedGenoPerByte_ = 4;                                    // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableHash::llWordSize_     = 8;                                    // 64 bit word size in bytes 
constexpr size_t   GenoTableHash::maxPairs_       = 6074000999UL;                         // approximate maximum number that does not overflow with n*(n-1)/2
constexpr uint64_t GenoTableHash::roundMask_      = 0xfffffffffffffff8;                   // mask for rounding down to nearest whole-byte value
constexpr uint64_t GenoTableHash::allBitsSet_     = std::numeric_limits<uint64_t>::max(); // 64-bit word with all bits set
constexpr size_t   GenoTableHash::wordSizeInBits_ = 64;                                   // 64-bit word size
constexpr size_t   GenoTableHash::nblocks32_      = sizeof(size_t) / sizeof(uint32_t);    // MurMurHash number of 32 bit blocks in size_t
constexpr uint32_t GenoTableHash::mmhKeyLen_      = sizeof(size_t);                       // MurMurHash key length 
constexpr uint16_t GenoTableHash::emptyBinToken_  = std::numeric_limits<uint16_t>::max(); // Value corresponding to an empty token 
constexpr uint32_t GenoTableHash::c1_             = 0xcc9e2d51;                           // MurMurHash c1 constant 
constexpr uint32_t GenoTableHash::c2_             = 0x1b873593;                           // MurMurHash c2 constant 

// Constructors
GenoTableHash::GenoTableHash(const std::string &inputFileName, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, const std::string &logFileName) : kSketches_{kSketches}, nLoci_{0}, nThreads_{nThreads}, logFileName_{logFileName} {
	std::stringstream logStream;
	const time_t t = time(nullptr);
	logStream << std::put_time(localtime(&t), "%b %e %Y %H:%M %Z");
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
	const size_t nBedBytes = nIndividuals / bedGenoPerByte_ + static_cast<bool>(nIndividuals % bedGenoPerByte_);
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
	if ( endPosition < static_cast<std::streamoff>( magicBytes_.size() ) ){
		logMessages_ += "ERROR: no loci in the input .bed file " + inputFileName + "; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	inStr.close();
	const size_t N  = static_cast<uint64_t>(endPosition) - magicBytes_.size();
	nLoci_          = N / nBedBytes;
	logMessages_   += "Number of individuals: " + std::to_string(nIndividuals) + "\n";
	logMessages_   += "Number of individuals to hash: " + std::to_string(nIndividuals_) + "\n";
	logMessages_   += "Number of loci: " + std::to_string(nLoci_) + "\n";
	logMessages_   += "Hash size: " + std::to_string(kSketches_) + "\n";
	locusSize_      = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;                    // round up to the nearest multiple of 8
	nFullWordBytes_ = (nIndividuals_ - 1) / byteSize_;
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	char magicBuf[magicBytes_.size()]{};
	inStr.read( magicBuf, magicBytes_.size() );
	if (magicBuf[0] != magicBytes_[0]){
		logMessages_ += "ERROR: file " + inputFileName + " does not appear to be in .bed format; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: first magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	} else if (magicBuf[1] != magicBytes_[1]){
		logMessages_ += "ERROR: file " + inputFileName + " does not appear to be in .bed format; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: second magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	} else if (magicBuf[2] != magicBytes_[2]){
		logMessages_ += "ERROR: file " + inputFileName + " does not appear to be in .bed format; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: third magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	// Generate the binary genotype table while reading the .bed file
	const size_t nBedBytesPerLocus = nIndividuals / bedGenoPerByte_ + static_cast<size_t>(nIndividuals % bedGenoPerByte_ > 0);
	const size_t ranVecSize        = nBedBytes / llWordSize_ + static_cast<size_t>(nBedBytes % llWordSize_ > 0);
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
		addIndv.emplace_back( std::pair<size_t, size_t>{iAddIndiv, rng_.sampleInt(nIndividuals)} );
	}
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{rng_.fyIndexesUp(nIndividuals_)};
	std::vector<uint32_t> seeds{static_cast<uint32_t>( rng_.ranInt() )};

	size_t locusInd = 0;
	if (nLociPerThread){
		std::vector< std::pair<size_t, size_t> > threadRanges;
		size_t bedInd = 0;
		for (size_t iThread = 0; iThread < nThreads_; ++iThread){
			threadRanges.emplace_back(std::pair<size_t, size_t>{bedInd, bedInd + nLociPerThread});
			bedInd += nLociPerThread;
		}
		const size_t excessLoci = nBedLociToRead - threadRanges.back().second;
		threadRanges.back().second = nBedLociToRead;
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			// TODO: add cassert() for nBedBytesToRead < streamsize::max()
			inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(nBedBytesToRead) );
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &tr : threadRanges){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, &tr, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &addIndv, &seeds]{
						bed2ophBlk_(bedChunkToRead, tr.first, tr.second, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, addIndv, seeds);
					})
				);
				locusInd += nLociPerThread;
			}
			for (const auto &th : tasks){
				th.wait();
			}
			locusInd += excessLoci;
		}
	} else {
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			// TODO: add cassert() for nBedBytesToRead < streamsize::max()
			inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(nBedBytesToRead) );
			std::vector< std::future<void> > tasks;
			tasks.reserve(nBedLociToRead);
			for (size_t iBedLocus = 0; iBedLocus < nBedLociToRead; ++iBedLocus){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, iBedLocus, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &addIndv, &seeds]{
							bed2ophBlk_(bedChunkToRead, iBedLocus, iBedLocus + 1, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, addIndv, seeds);
						})
				);
				++locusInd;
			}
			for (const auto &th : tasks){
				th.wait();
			}
		}
	}
	if (remainingLoci) {
		bedChunkToRead.resize(remainingBytes);
		// TODO: add cassert() for remainingBytes < streamsize::max()
		inStr.read( bedChunkToRead.data(), static_cast<std::streamsize>(remainingBytes) );
		nLociPerThread = remainingLoci / nThreads_;
		if (nLociPerThread){
			std::vector< std::pair<size_t, size_t> > threadRanges;
			size_t bedInd = 0;
			for (size_t iThread = 0; iThread < nThreads_; ++iThread){
				threadRanges.emplace_back(std::pair<size_t, size_t>{bedInd, bedInd + nLociPerThread});
				bedInd += nLociPerThread;
			}
			const size_t excessLoci = remainingLoci - threadRanges.back().second;
			threadRanges.back().second = remainingLoci;
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &tr : threadRanges){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, &tr, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &addIndv, &seeds]{
						bed2ophBlk_(bedChunkToRead, tr.first, tr.second, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, addIndv, seeds);
					})
				);
				locusInd += nLociPerThread;
			}
			for (const auto &th : tasks){
				th.wait();
			}
			locusInd += excessLoci;
		} else {
			std::vector< std::future<void> > tasks;
			tasks.reserve(remainingLoci);
			for (size_t iBedLocus = 0; iBedLocus < remainingLoci; ++iBedLocus){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, iBedLocus, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &addIndv, &seeds]{
							bed2ophBlk_(bedChunkToRead, iBedLocus, iBedLocus + 1, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, addIndv, seeds);
						})
				);
				++locusInd;
			}
			for (const auto &th : tasks){
				th.wait();
			}
		}
	}
	inStr.close();
}

GenoTableHash::GenoTableHash(const std::vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, const std::string &logFileName) : nIndividuals_{nIndividuals}, kSketches_{kSketches}, nLoci_{maCounts.size() / nIndividuals}, nThreads_{nThreads}, logFileName_{logFileName} {
	std::stringstream logStream;
	const time_t t = time(nullptr);
	logStream << std::put_time(localtime(&t), "%b %e %H:%M %Z");
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
	if (maCounts.size() % nIndividuals){
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
	sketchSize_ = nIndividuals_ / kSketches_;    // want only full sketches (all from the same number of individuals)
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
	locusSize_              = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	const size_t ranVecSize = locusSize_ / llWordSize_ + static_cast<bool>(locusSize_ % llWordSize_);
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

	std::vector< std::pair<size_t, size_t> > threadRanges;
	const size_t nLociPerThread = nLoci_ / nThreads_;
	if (nLociPerThread){
		size_t locusInd = 0;
		for (size_t iThread = 0; iThread < nThreads_; ++iThread){
			threadRanges.emplace_back(std::pair<size_t, size_t>{locusInd, locusInd + nLociPerThread});
			locusInd += nLociPerThread;
		}
		threadRanges.back().second = nLoci_;
		std::vector< std::future<void> > tasks;
		tasks.reserve(nThreads_);
		for (const auto &tr : threadRanges){
			tasks.emplace_back(
				std::async([this, &maCounts, &tr, ranVecSize, &ranInts, &seeds]{
					mac2ophBlk_(maCounts, tr.first, tr.second, ranVecSize, ranInts, seeds);
				})
			);
		}
		for (const auto &th : tasks){
			th.wait();
		}
	} else {
		mac2ophBlk_(maCounts, 0, nLoci_, ranVecSize, ranInts, seeds);
	}
}

GenoTableHash::GenoTableHash(GenoTableHash &&in) noexcept {
	*this = std::move(in);
}

GenoTableHash& GenoTableHash::operator=(GenoTableHash &&in) noexcept {
	// TODO: add the new variables
	if (this != &in){
		sketches_     = std::move(in.sketches_);
		nIndividuals_ = in.nIndividuals_;
		kSketches_    = in.kSketches_;
		nLoci_        = in.nLoci_;
		sketchSize_   = in.sketchSize_;
		nThreads_     = in.nThreads_;
		locusSize_    = in.locusSize_;
		logFileName_  = std::move(in.logFileName_);
		logMessages_  = std::move(in.logMessages_);

		in.nIndividuals_ = 0;
		in.kSketches_    = 0;
		in.nLoci_        = 0;
		in.sketchSize_   = 0;
		in.nThreads_     = 0;
		in.locusSize_    = 0;
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
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * sizeof(float) );      // use half to leave resources for other operations
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
		const size_t nLocusPairsPerThread = maxInRAM / nThreads_;
		size_t overallPairInd = 0;
		if (nLocusPairsPerThread){
			std::vector< std::pair<size_t, size_t> > threadRanges;
			size_t locusPairInd = 0;
			for (size_t iThread = 0; iThread < nThreads_; ++iThread){
				threadRanges.emplace_back(std::pair<size_t, size_t>{locusPairInd, locusPairInd + nLocusPairsPerThread});
				locusPairInd += nLocusPairsPerThread;
			}
			const size_t excessLoci    = maxInRAM - threadRanges.back().second;
			threadRanges.back().second = maxInRAM;
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<float> LDmatChunk(maxInRAM, 0.0);
				std::vector< std::future<void> > tasks;
				tasks.reserve(nThreads_);
				for (const auto &tr : threadRanges){
					tasks.emplace_back(
						std::async([this, &LDmatChunk, &tr, overallPairInd]{
							hashJacBlock_(tr.first, tr.second, overallPairInd, LDmatChunk);
						})
					);
					overallPairInd += nLocusPairsPerThread;
				}
				overallPairInd += excessLoci;
				for (const auto &th : tasks){
					th.wait();
				}
				for (const auto &ld : LDmatChunk){
					output << ld << " ";
				}
			}
		} else {
			for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
				std::vector<float> LDmatChunk(maxInRAM, 0.0);
				hashJacBlock_(0, maxInRAM, overallPairInd, LDmatChunk);
				overallPairInd += maxInRAM;
				for (const auto &ld : LDmatChunk){
					output << ld << " ";
				}
			}
		}
		overallPairInd = maxInRAM * nChunks;
		if (remainingPairs){
			std::vector<float> LDmatChunk(remainingPairs, 0.0);
			const size_t nRemainPairsPerThread = remainingPairs / nThreads_;
			if (nRemainPairsPerThread){
				std::vector< std::pair<size_t, size_t> > threadRanges;
				size_t locusPairInd = 0;
				for (size_t iThread = 0; iThread < nThreads_; ++iThread){
					threadRanges.emplace_back(std::pair<size_t, size_t>{locusPairInd, locusPairInd + nRemainPairsPerThread});
					locusPairInd += nRemainPairsPerThread;
				}
				threadRanges.back().second = remainingPairs;
				std::vector< std::future<void> > tasks;
				tasks.reserve(nThreads_);
				for (const auto &tr : threadRanges){
					tasks.emplace_back(
						std::async([this, &LDmatChunk, &tr, overallPairInd]{
							hashJacBlock_(tr.first, tr.second, overallPairInd, LDmatChunk);
						})
					);
					overallPairInd += nRemainPairsPerThread;
				}
			} else {
				hashJacBlock_(0, remainingPairs, overallPairInd, LDmatChunk);
			}
			for (const auto &ld : LDmatChunk){
				output << ld << " ";
			}
		}
		output.close();
	} else {
		logMessages_ += "calculating in one go\n";
		std::vector<float> LDmat(nPairs, 0.0);
		const size_t nLocusPairsPerThread = LDmat.size() / nThreads_;
		if (nLocusPairsPerThread){
			std::vector< std::pair<size_t, size_t> > threadRanges;
			size_t locusPairInd = 0;
			for (size_t iThread = 0; iThread < nThreads_; ++iThread){
				threadRanges.emplace_back(std::pair<size_t, size_t>{locusPairInd, locusPairInd + nLocusPairsPerThread});
				locusPairInd += nLocusPairsPerThread;
			}
			threadRanges.back().second = LDmat.size();
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &tr : threadRanges){
				tasks.emplace_back(
					std::async([this, &LDmat, &tr]{
						hashJacBlock_(tr.first, tr.second, tr.first, LDmat);
					})
				);
			}
		} else {
			hashJacBlock_(0, LDmat.size(), 0, LDmat);
		}
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		for (const auto &ld : LDmat){
			output << ld << " ";
		}
		output.close();
	}
}

std::unordered_map< uint32_t, std::vector<size_t> > GenoTableHash::makeLDgroups(const size_t &nRowsPerBand) const {
	// TODO: add an assert() for nRowsPerBand != 0 and nRowsPerBand < kSketches_
	const size_t nBands = kSketches_ / nRowsPerBand; // only using full-size bands because smaller ones permit inclusion of low-similarity pairs

	logMessages_ += "Grouping loci\n";
	logMessages_ += "Number of rows per band: " + std::to_string(nRowsPerBand) + "\n";
	logMessages_ += "Number of bands: " + std::to_string(nBands) + "\n";

	const uint32_t sketchSeed = static_cast<uint32_t>( rng_.ranInt() );
	const uint32_t indexSeed  = static_cast<uint32_t>( rng_.ranInt() );
	std::unordered_map< uint32_t, std::vector<size_t> > ldGroup;       // the hash table; indexed by hashing the index vector

	size_t iSketch = 0;
	for (size_t iBand = 0; iBand < nBands; ++iBand){
		std::unordered_map< uint32_t, std::vector<size_t> > localLDG;  // hash table local to each band, indexed by sketch hashes
		for (size_t iLocus = 0; iLocus < nLoci_; ++iLocus){
		 	const uint32_t hash = murMurHash_(iSketch + iLocus * kSketches_, nRowsPerBand, sketchSeed);
			localLDG[hash].push_back(iLocus);
		}
		for (auto el : localLDG){
			const uint32_t elementHash = murMurHash_(el.second, indexSeed);
			if ( ldGroup[elementHash].empty() ){
				ldGroup[elementHash] = std::move(el.second);
			}
		}
		iSketch += nRowsPerBand;
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
	return ldGroup;
}

void GenoTableHash::ldInGroups(const size_t &nRowsPerBand, const std::string &outFileName) const {
	std::unordered_map< uint32_t, std::vector<size_t> > ldGroup = this->makeLDgroups(nRowsPerBand);
	
	// there is a possibility that all retained pairs will not fit in RAM
	// dealing with this is non-trivial (need the whole thing in RAM to eliminate duplicate pairs), so for now will let the system deal with memory issues
	logMessages_ += "Estimating LD in groups\n";
	
	// identify index pairs that are retained in groups
	// some will be repeatedly identified by makeLDgroups(), so we will need to eliminate duplicates
	std::vector<IndexedPairSimilarity> hashJacGroups;
	uint32_t groupInd = 0;
	for (auto &ldg : ldGroup){
		for (size_t iLocus = 0; iLocus < ldg.second.size() - 1; ++iLocus){
			for (size_t jLocus = iLocus + 1; jLocus < ldg.second.size(); ++jLocus){
				hashJacGroups.emplace_back(IndexedPairSimilarity{0.0, ldg.second[iLocus], ldg.second[jLocus], groupInd});
			}
		}
		++groupInd;
		ldg.second.clear();
	}
	logMessages_ += "Number of locus pairs before removing duplicates: " + std::to_string( hashJacGroups.size() )+ "\n";
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
	if (nPairsPerThread){
		std::vector< std::pair<size_t, size_t> > threadRanges;
		threadRanges.reserve(nThreads_);
		size_t pairInd = 0;
		for (size_t iThread = 0; iThread < nThreads_; ++iThread){
			threadRanges.emplace_back(std::pair<size_t, size_t>{pairInd, pairInd + nPairsPerThread});
			pairInd += nPairsPerThread;
		}
		threadRanges.back().second = hashJacGroups.size(); // assuming the number of threads << number of groups, this should not lead to noticeable thread imbalance
		std::vector< std::future<void> > tasks;
		tasks.reserve(nThreads_);
		for (const auto &tr : threadRanges){
			tasks.emplace_back(
				std::async([this, &tr, &hashJacGroups]{
					hashJacBlock_(tr.first, tr.second, hashJacGroups);
				})
			);
		}
		for (const auto &th : tasks){
			th.wait();
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
		for (const auto &th : tasks){
			th.wait();
		}
	}

	std::fstream out;
	out.open(outFileName, std::ios::out | std::ios::trunc);
	out << "groupID\tlocus1\tlocus2\tjaccLD\n";
	for (const auto idxRes : hashJacGroups){
		out << "G" << idxRes.groupID + 1 << "\t" << idxRes.element1ind + 1 << "\t" << idxRes.element2ind + 1 << "\t" << idxRes.similarityValue << "\n";
	}
	out.close();
}

void GenoTableHash::saveLogFile() const {
	std::fstream outLog;
	outLog.open(logFileName_, std::ios::out | std::ios::trunc);
	outLog << logMessages_;
	outLog.close();
}

void GenoTableHash::locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds, std::vector<uint8_t> &binLocus){
	// Start with a permutation to make OPH
	size_t iIndiv = 0;
	size_t iByte  = 0;
	while(iByte < nFullWordBytes_){
		for (uint8_t iInLocusByte = 0; iInLocusByte < byteSize_; ++iInLocusByte){
			auto bytePair            = static_cast<uint16_t>(binLocus[iByte]);
			const size_t perIndiv    = permutation[iIndiv++];                                   // post-increment to use current value for index first
			const size_t permByteInd = perIndiv / byteSize_;
			const auto permInByteInd = static_cast<uint8_t>( perIndiv - (perIndiv & roundMask_) );
			// Pair the current locus byte with the byte containing the value to be swapped
			// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
			bytePair                    |= static_cast<uint16_t>(binLocus[permByteInd]) << byteSize_;
			const auto mask              = static_cast<uint16_t>(oneBit_ << iInLocusByte);
			const auto perMask           = static_cast<uint8_t>(oneBit_ << permInByteInd);
			const uint16_t shiftDistance = (byteSize_ - iInLocusByte) + permInByteInd;           // subtraction is safe b/c byteSize is the loop terminator
			const uint16_t temp1         = ( bytePair ^ (bytePair >> shiftDistance) ) & mask;
			const auto temp2             = static_cast<uint16_t>(temp1 << shiftDistance);
			bytePair                    ^= temp1 ^ temp2;
			// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus)
			// Must modify the current byte in each loop iteration because permutation indexes may fall into it
			binLocus[iByte]             ^= ( binLocus[iByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask);
			binLocus[permByteInd]       ^= ( binLocus[permByteInd] ^ static_cast<uint8_t>(bytePair >> byteSize_) ) & perMask;
 		}
		++iByte;
	}
	// Finish the individuals in the remaining partial byte, if any
	uint8_t iInLocusByte = 0;
	while (iIndiv < nIndividuals_ - 1){
		auto bytePair            = static_cast<uint16_t>(binLocus[iByte]);
		const size_t perIndiv    = permutation[iIndiv++];                                   // post-increment to use current value for index first
		const size_t permByteInd = perIndiv / byteSize_;
		const auto permInByteInd = static_cast<uint8_t>( perIndiv - (perIndiv & roundMask_) );
		// Pair the current locus byte with the byte containing the value to be swapped
		// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
		bytePair                    |= static_cast<uint16_t>(binLocus[permByteInd]) << byteSize_;
		const auto mask              = static_cast<uint16_t>(oneBit_ << iInLocusByte);
		const auto perMask           = static_cast<uint8_t>(oneBit_ << permInByteInd);
		const uint16_t shiftDistance = (byteSize_ - iInLocusByte) + permInByteInd;           // subtraction is safe b/c byteSize is the loop terminator
		const uint16_t temp1         = ( bytePair ^ (bytePair >> shiftDistance) ) & mask;
		const auto temp2             = static_cast<uint16_t>(temp1 << shiftDistance);
		bytePair                    ^= temp1 ^ temp2;
		// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus)
		// Must modify the current byte in each loop iteration because permutation indexes may fall into it
		binLocus[iByte]       ^= ( binLocus[iByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask);
		binLocus[permByteInd] ^= ( binLocus[permByteInd] ^ static_cast<uint8_t>(bytePair >> byteSize_) ) & perMask;
		++iInLocusByte;
	}
	// Now make the sketches
	std::vector<size_t> filledIndexes;                                                       // indexes of the non-empty sketches
	iByte = 0;
	size_t sketchBeg{locusInd * kSketches_};
	size_t iSketch{0};
	uint64_t sketchTail{0};                                                                  // left over buts from beyond the last full byte of the previous sketch
	size_t locusChunkSize = (llWordSize_ > binLocus.size() ? binLocus.size() : llWordSize_);
	while ( iByte < binLocus.size() ){
		uint64_t nWordUnsetBits{wordSizeInBits_};
		uint64_t nSumUnsetBits{0};
		while ( (nWordUnsetBits == wordSizeInBits_) && ( iByte < binLocus.size() ) ){
			uint64_t locusChunk{allBitsSet_};
			// TODO: add cassert() for iByte < binLocus.size(); must be true since this is the loop conditions
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
		while (emptyCount){
			for (const auto &f : filledIndexes){
				auto newIdx = static_cast<uint32_t>(murMurHash_(f, seeds[iSeed]) % kSketches_ + sketchBeg);
				// should be safe: each thread accesses different vector elements
				if (sketches_[newIdx] == emptyBinToken_){
					sketches_[newIdx] = sketches_[f + sketchBeg];
					--emptyCount;
					break;
				}
			}
			++iSeed;
			std::lock_guard<std::mutex> lk(mtx_);      // lock before measuring to ensure that the size is valid
			if ( iSeed == seeds.size() ){
				seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
			}
		}
	}
}

void GenoTableHash::bed2ophBlk_(const std::vector<char> &bedData, const size_t &firstBedLocusInd, const size_t &lastBedLocusInd, const size_t &firstLocusInd,
					const size_t &bedLocusLength, const size_t &randVecLen, const std::vector<size_t> &permutation, std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint64_t locSeed;
	{
		std::lock_guard<std::mutex> lk(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
	const size_t addIndv    = nIndividuals_ % 4UL;
	const size_t addBL      = bedLocusLength - static_cast<size_t>(addIndv > 0);
	const auto lastByteMask = static_cast<uint8_t>( std::numeric_limits<uint8_t>::max() << static_cast<uint8_t>(nIndividuals_ % byteSize_) );
	// Fill the random byte vector
	std::vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	size_t iLocus      = firstLocusInd;
	for (size_t iBedLocus = firstBedLocusInd; iBedLocus < lastBedLocusInd; ++iBedLocus){
		const size_t begInd      = iBedLocus * bedLocusLength;
		const size_t endWholeBed = begInd + addBL;
		for (auto &rv : rand){
			rv = locPRNG.ranInt();
		}
		std::vector<uint8_t> binLocus(locusSize_, 0);
		std::vector<uint8_t> missMasks(locusSize_, 0);
		size_t iBinGeno = 0;                       // binLocus vector index
		uint8_t bedByte = 0;
		size_t iRB      = 0;                       // random byte index
		// Two bytes of .bed code go into one byte of my binary representation
		// Therefore, work on two consecutive bytes of .bed code in the loop
		for (size_t iBed = begInd; iBed < endWholeBed; iBed += 2){                 // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
			bedByte = static_cast<uint8_t>(~bedData[iBed]);                        // flip so that homozygous second allele (usually minor) is set to 11
			uint8_t offsetToBin{0};                                                // move the .bed mask by this much to align with the binarized byte
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks[iBinGeno]      |= curMissMask >> offsetToBin;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask      |= randBytes[iRB] & (firstBitMask << 1);
				firstBitMask       &= secondBitMask >> 1;
				binLocus[iBinGeno] |= firstBitMask >> offsetToBin;
				++offsetToBin;
			}
			const size_t nextIbed = iBed + 1;
			bedByte = static_cast<uint8_t>(~bedData[nextIbed]);
			++iRB;
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks[iBinGeno]      |= curMissMask << offsetToBin;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask      |= randBytes[iRB] & (firstBitMask << 1);
				firstBitMask       &= secondBitMask >> 1;
				binLocus[iBinGeno] |= firstBitMask << offsetToBin; // keep adding to the current binarized byte, so switch the direction of shift
				--offsetToBin;
			}
			++iBinGeno;
			++iRB;
		}
		uint8_t inBedByteOffset = 0;
		if (addIndv) {
			const uint8_t lastBedByte = static_cast<uint8_t>(~bedData[endWholeBed]);
			for (size_t iInd = 0; iInd < addIndv; ++iInd){
				uint8_t firstBitMask    = lastBedByte & (oneBit_ << inBedByteOffset);
				const uint8_t secondBBO = inBedByteOffset + 1;
				uint8_t secondBitMask   = lastBedByte & (oneBit_ << secondBBO);
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks.back()         |= (curMissMask >> inBedByteOffset) << iInd;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask   |= randBytes[iRB] & (firstBitMask << 1);                                       // iRB incremented at the end of the previous loop
				firstBitMask    &= secondBitMask >> 1;
				firstBitMask     = firstBitMask >> inBedByteOffset;
				firstBitMask     = static_cast<uint8_t>(firstBitMask << iInd);
				binLocus.back() |= firstBitMask;
				inBedByteOffset += 2;
				inBedByteOffset  = inBedByteOffset % 8;
			}
		}
		for (const auto &addI : padIndiv){
			const size_t iLocByte    = addI.first / byteSize_;
			const auto iInLocByte    = static_cast<uint8_t>(addI.first % byteSize_);
			auto bytePair            = static_cast<uint16_t>(binLocus[iLocByte]);
			const size_t permByteInd = addI.second / byteSize_;
			const auto permInByteInd = static_cast<uint8_t>(addI.second % byteSize_);
			// Pair the current locus byte with the byte containing the value to be swapped
			// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
			bytePair                    |= static_cast<uint16_t>(binLocus[permByteInd]) << byteSize_;
			const auto mask              = static_cast<uint16_t>(1 << iInLocByte);
			const uint16_t shiftDistance = (byteSize_ - iInLocByte) + permInByteInd;                        // subtraction is safe b/c iInLocByte is modulo byteSize
			const uint16_t temp1         = ( bytePair ^ (bytePair >> shiftDistance) ) & mask;
			const auto temp2             = static_cast<uint16_t>(temp1 << shiftDistance);
			bytePair                    ^= temp1 ^ temp2;
			// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus1)
			// Must modify the current byte in each loop iteration because permutation indexes may fall into it
			binLocus[iLocByte]         ^= ( binLocus[iLocByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask);
		}
		float aaCount = static_cast<float>( countSetBits(binLocus) ) / static_cast<float>(nIndividuals_);
		if (aaCount > 0.5){ // always want the alternative to be the minor allele
			for (size_t iBL = 0; iBL < locusSize_; ++iBL){
				binLocus[iBL] = (~binLocus[iBL]) & (~missMasks[iBL]);
			}
		}
		if ( lastByteMask < std::numeric_limits<uint8_t>::max() ){
			binLocus.back() |= lastByteMask;
		}
		locusOPH_(iLocus, permutation, seeds, binLocus);
		++iLocus;
	}
}

void GenoTableHash::mac2ophBlk_(const std::vector<int> &macData, const size_t &startLocusInd, const size_t &endLocusInd, const size_t &randVecLen, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint64_t locSeed;
	{
		std::lock_guard<std::mutex> lk(mtx_);
		locSeed = rng_.ranInt();
	}
	RanDraw locPRNG(locSeed);
	uint8_t remainderInd       = static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
	remainderInd               = byteSize_ - remainderInd;
	// Create a vector to store random bytes for stochastic heterozygote resolution
	std::vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	for (size_t iLocus = startLocusInd; iLocus < endLocusInd; ++iLocus){
		// Fill the random byte vector
		for (auto &rv : rand){
			rv = locPRNG.ranInt();
		}
		size_t iIndiv         = 0;
		const size_t begIndiv = iLocus * nIndividuals_;
		std::vector<uint8_t> missMasks(locusSize_, 0);
		std::vector<uint8_t> binLocus(locusSize_, 0);
		size_t i0Byte = 0;                                                                         // to index the missMasks vector
		for (size_t iByte = 0; iByte < locusSize_ - 1; ++iByte){                                   // treat the last byte separately
			for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte){
				uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);      // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= 0b10000011;                                            // mask everything in the middle
				const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
				missMasks[i0Byte]        |= (missingMask << iInByte);
				curIndiv                 &= 0b00000011;
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
			uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);          // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= 0b10000011;                                                // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
			missMasks.back()         |= (missingMask << iRem);
			curIndiv                 &= 0b00000011;                                                
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

uint32_t GenoTableHash::murMurHash_(const size_t &key, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	auto blocks = reinterpret_cast<const uint32_t *>(&key);

	for (size_t i = 0; i < nblocks32_; ++i){
		uint32_t k1 = blocks[i];

		k1 *= c1_;
		k1  = (k1 << 15) | (k1 >> 17);
		k1 *= c2_;

		hash ^= k1;
		hash  = (hash << 13) | (hash >> 19);
		hash  = hash * 5 + 0xe6546b64;
	}

	// there is no tail since the input is fixed (at eight bytes typically)
	// finalize
	hash ^= mmhKeyLen_;
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;

	return hash;
}

uint32_t GenoTableHash::murMurHash_(const std::vector<size_t> &key, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	auto blocks = reinterpret_cast<const uint32_t *>( key.data() );
	const size_t nBlocks = nblocks32_ * key.size();

	for (size_t i = 0; i < nBlocks; ++i){
		uint32_t k1 = blocks[i];

		k1 *= c1_;
		k1  = (k1 << 15) | (k1 >> 17);
		k1 *= c2_;

		hash ^= k1;
		hash  = (hash << 13) | (hash >> 19);
		hash  = hash * 5 + 0xe6546b64;
	}

	// there is no tail since the input is fixed (at eight bytes typically)
	// finalize
	hash ^= mmhKeyLen_;
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;

	return hash;
}

uint32_t GenoTableHash::murMurHash_(const size_t &startInd, const size_t &nElements, const uint32_t &seed) const {
	// TODO: add an assert() on nElements != 0 and startInd < sketches_.size()
	uint32_t hash  = seed;
	auto blocks    = reinterpret_cast<const uint32_t *>(sketches_.data() + startInd);
	size_t nBlocks = nElements / 2; // each sketch is 16 bits

	// body
	for (size_t iBlock = 0; iBlock < nBlocks; ++iBlock){
		uint32_t k1 = blocks[iBlock];

		k1 *= c1_;
		k1  = (k1 << 15) | (k1 >> 17);
		k1 *= c2_;

		hash ^= k1;
		hash  = (hash << 13) | (hash >> 19);
		hash  = hash * 5 + 0xe6546b64;
	}

	// tail, if exists
	if (nElements % 2){
		uint32_t k1 = static_cast<uint32_t>(sketches_[startInd + nElements - 1]);

		k1 *= c1_;
		k1  = (k1 << 15) | (k1 >> 17);
		k1 *= c2_;

		hash ^= k1;
		hash  = (hash << 13) | (hash >> 19);
		hash  = hash * 5 + 0xe6546b64;
	}

	// finalize
	hash ^= mmhKeyLen_;
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;

	return hash;
}

void GenoTableHash::hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, const size_t &blockStartAll, std::vector<float> &hashJacVec) const {
	const size_t nnLoci = nLoci_ * (nLoci_ - 1) / 2 - 1; // overflow checked in the calling function; do this here for encapsulation, may move to the calling function later
	size_t curJacMatInd = blockStartAll;
	for (size_t iVecInd = blockStartVec; iVecInd < blockEndVec; ++iVecInd){
		// compute row and column indexes from the vectorized by column lower triangle index
		// got these expressions by combining various web sources and verifying
		const size_t kp  = nnLoci - curJacMatInd;
		const size_t p   = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kp) ) ) - 1) / 2;
		const size_t row = nLoci_ - 2 - p;
		const size_t col = nLoci_ - (kp - p * (p + 1) / 2) - 1;
		const auto rowSk = static_cast<std::vector<uint16_t>::difference_type>(row * kSketches_);
		const auto colSk = static_cast<std::vector<uint16_t>::difference_type>(col * kSketches_);
		auto start       = sketches_.begin() + rowSk;
		// count equal elements using the inner_product idiom
		int simVal = std::inner_product( start, start + static_cast<std::vector<uint16_t>::difference_type>(kSketches_),
				sketches_.begin() + colSk, 0, std::plus<>(), std::equal_to<>() );
		// should be safe: each thread accesses different vector elements
		hashJacVec[iVecInd] = static_cast<float>(simVal) / static_cast<float>(kSketches_);
		++curJacMatInd;
	}
}

void GenoTableHash::hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, std::vector<IndexedPairSimilarity> &hashJacVec) const {
	for (size_t iBlock = blockStartVec; iBlock < blockEndVec; ++iBlock){
		const auto start1 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(hashJacVec[iBlock].element1ind * kSketches_);
		const auto start2 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(hashJacVec[iBlock].element2ind * kSketches_);
		// count equal elements using the inner_product idiom
		int simVal = std::inner_product( start1, start1 + static_cast<std::vector<uint16_t>::difference_type>(kSketches_),
				start2, 0, std::plus<>(), std::equal_to<>() );
		// should be safe: each thread accesses different vector elements
		hashJacVec[iBlock].similarityValue = static_cast<float>(simVal) / static_cast<float>(kSketches_);
	}
}

