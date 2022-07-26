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
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <utility>  // for std::pair
#include <iterator>
#include <algorithm>
#include <limits>
#include <future>
#include <thread>
#include <ctime>
#include <iomanip>

#include "gvarHash.hpp"

using namespace BayesicSpace;
// External functions
uint16_t BayesicSpace::countSetBits(uint16_t inVal) {
	uint16_t totSet = 0;
	for (; inVal; ++totSet) {
		inVal &= inVal - 1;
	}
	return totSet;
}

uint32_t BayesicSpace::countSetBits(const std::vector<uint8_t> &inVec) {
	uint32_t totSet = 0;
	for (const auto &in : inVec){
		uint8_t v = in;
		for (; v; ++totSet) {
			v &= v - 1;
		}
	}
	return totSet;
}

uint32_t BayesicSpace::countSetBits(const std::vector<uint8_t> &inVec, const size_t &start, const size_t &length) {
	uint32_t totSet = 0;
	for (size_t i = start; i < start + length; ++i){
		uint8_t v = inVec[i];
		for (; v; ++totSet) {
			v &= v - 1;
		}
	}
	return totSet;
}

size_t BayesicSpace::getAvailableRAM() {
	if ( std::ifstream("/proc/meminfo").good() ) {
		std::string memLine;
		std::fstream memInfoStream;
		memInfoStream.open("/proc/meminfo", std::ios::in);
		while ( getline(memInfoStream, memLine) ){
			if (memLine.compare(0, 13, "MemAvailable:") == 0) {
				break;
			}
		}
		memInfoStream.close();
		std::stringstream memLineStream(memLine);
		std::string freeMemStr;
		memLineStream >> freeMemStr;
		memLineStream >> freeMemStr;
		return static_cast<size_t>( stoi(freeMemStr) ) * 1024UL; // memory is in kB in the file
	} else {
		return 2147483648UL;
	}
}

// GenoTableBin methods
constexpr std::array<char, 3> GenoTableBin::magicBytes_ = {0x6c, 0x1b, 0x01};   // Leading bytes for .bed files
constexpr uint8_t  GenoTableBin::oneBit_                = 0b00000001;           // One set bit for masking
constexpr uint8_t  GenoTableBin::byteSize_              = 8;                    // Size of one byte in bits
constexpr uint8_t  GenoTableBin::llWordSize_            = 8;                    // 64 bit word size in bytes
constexpr size_t   GenoTableBin::maxNlocusPairs_        = 6074000999UL;         // approximate number of loci that does not overflow with n*(n-1)/2

// Constructors
GenoTableBin::GenoTableBin(const std::string &inputFileName, const size_t &nIndividuals, const size_t &nThreads) : nIndividuals_{nIndividuals}, nThreads_{nThreads} {
	if (nIndividuals <= 1){
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string(__FUNCTION__);
	} else if (nIndividuals > std::numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
		throw std::string("ERROR: the number of individuals (") + std::to_string(nIndividuals) + std::string( ") is too big to make a square relationship matrix in ") + std::string(__FUNCTION__);
	}
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ){
		nThreads_ = std::thread::hardware_concurrency();
	}
	const size_t nBedBytesPerLocus = nIndividuals_ / 4 + static_cast<bool>(nIndividuals_ % 4);
	std::fstream inStr;
	// Start by measuring file size
	inStr.open(inputFileName, std::ios::in | std::ios::binary | std::ios::ate);
	if ( inStr.fail() ){
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string(__FUNCTION__);
	}
	const int32_t endPosition = inStr.tellg();
	if ( endPosition < magicBytes_.size() ){
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string(__FUNCTION__);
	}
	const size_t N = static_cast<uint64_t>(endPosition) - magicBytes_.size();
	inStr.close();
	nLoci_ = N / nBedBytesPerLocus;

	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	char magicBuf[magicBytes_.size()]{};
	inStr.read( magicBuf, magicBytes_.size() );
	if (magicBuf[0] != magicBytes_[0]){
		throw std::string("ERROR: first magic byte in input .bed file is not the expected value in ") + std::string(__FUNCTION__);
	} else if (magicBuf[1] != magicBytes_[1]){
		throw std::string("ERROR: second magic byte in input .bed file is not the expected value in ") + std::string(__FUNCTION__);
	} else if (magicBuf[2] != magicBytes_[2]){
		throw std::string("ERROR: third magic byte in input .bed file is not the expected value in ") + std::string(__FUNCTION__);
	}
	// Generate the binary genotype table while reading the .bed file
	binLocusSize_           = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	const size_t ranVecSize = nBedBytesPerLocus / llWordSize_ + static_cast<bool>(nBedBytesPerLocus % llWordSize_);
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);
	const size_t ramSize         = getAvailableRAM() / 2UL;                               // measuring here, after all the major allocations; use half to leave resources for other operations
	size_t nBedLociToRead        = ramSize / nBedBytesPerLocus;                           // number of .bed loci to read at a time
	nBedLociToRead               = (nBedLociToRead < nLoci_ ? nBedLociToRead : nLoci_);
	const size_t remainingLoci   = nLoci_ % nBedLociToRead;
	const size_t remainingBytes  = remainingLoci * nBedBytesPerLocus;
	const size_t nChunks         = nLoci_ / nBedLociToRead;
	const size_t nBedBytesToRead = nBedLociToRead * nBedBytesPerLocus;
	size_t nLociPerThread        = nBedLociToRead / nThreads_;
	std::vector<char> bedChunkToRead(nBedBytesToRead, 0);

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
			inStr.read(bedChunkToRead.data(), nBedBytesToRead);
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &tr : threadRanges){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, &tr, locusInd, nBedBytesPerLocus, ranVecSize]{
						bed2binBlk_(bedChunkToRead, tr.first, tr.second, locusInd, nBedBytesPerLocus, ranVecSize);
					})
				);
				locusInd += nLociPerThread;
			}
			locusInd += excessLoci;
		}
	} else {
		for (size_t iChunk = 0; iChunk < nChunks; ++iChunk){
			inStr.read(bedChunkToRead.data(), nBedBytesToRead);
			std::vector< std::future<void> > tasks;
			tasks.reserve(nBedLociToRead);
			for (size_t iBedLocus = 0; iBedLocus < nBedLociToRead; ++iBedLocus){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, iBedLocus, locusInd, nBedBytesPerLocus, ranVecSize]{
						bed2binBlk_(bedChunkToRead, iBedLocus, iBedLocus + 1, locusInd, nBedBytesPerLocus, ranVecSize);})
				);
				++locusInd;
			}
		}
	}
	if (remainingLoci) {
		bedChunkToRead.resize(remainingBytes);
		inStr.read(bedChunkToRead.data(), remainingBytes);
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
					std::async([this, &bedChunkToRead, &tr, locusInd, nBedBytesPerLocus, ranVecSize]{
						bed2binBlk_(bedChunkToRead, tr.first, tr.second, locusInd, nBedBytesPerLocus, ranVecSize);
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
						bed2binBlk_(bedChunkToRead, iBedLocus, iBedLocus + 1, locusInd, nBedBytesPerLocus, ranVecSize);})
				);
				++locusInd;
			}
		}
	}
	inStr.close();
}

GenoTableBin::GenoTableBin(const std::vector<int> &maCounts, const size_t &nIndividuals, const size_t &nThreads) : nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals}, nThreads_{nThreads} {
	if (nIndividuals <= 1){
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string(__FUNCTION__);
	}
	if (maCounts.size() % nIndividuals){
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(nIndividuals) + std::string( ") in ") + std::string(__FUNCTION__);
	}
	if ( maCounts.empty() ){
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string(__FUNCTION__);
	}
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ) {
		nThreads_ = std::thread::hardware_concurrency();
	}
	binLocusSize_ = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);
	const size_t ranVecSize     = nIndividuals_ / llWordSize_ + static_cast<bool>(nIndividuals_ % llWordSize_);
	const size_t nLociPerThread = nLoci_ / nThreads_;
	if (nLociPerThread){
		std::vector< std::pair<size_t, size_t> > threadRanges;
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
				std::async([this, &maCounts, &tr, ranVecSize]{
					mac2binBlk_(maCounts, tr.first, tr.second, ranVecSize);
				})
			);
		}
	} else {
		mac2binBlk_(maCounts, 0, nLoci_, ranVecSize);
	}
}

GenoTableBin::GenoTableBin(GenoTableBin &&in) noexcept {
	*this = std::move(in);
}

GenoTableBin& GenoTableBin::operator=(GenoTableBin &&in) noexcept {
	if (this != &in){
		binGenotypes_ = move(in.binGenotypes_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;
		binLocusSize_ = in.binLocusSize_;
		nThreads_     = in.nThreads_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;
		in.binLocusSize_ = 0;
		in.nThreads_     = 0;
	}
	return *this;
}

void GenoTableBin::saveGenoBinary(const std::string &outFileName) const {
	std::fstream out;
	out.open(outFileName, std::ios::out | std::ios::binary | std::ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), binGenotypes_.size() );
	out.close();
}

void GenoTableBin::allJaccardLD(const std::string &ldFileName) const {
	if (nLoci_ > maxNlocusPairs_){
		throw std::string("ERROR: Too many loci (") + std::to_string(nLoci_) + std::string(") to do all pairwise LD. Maximum supported is ") +
			 std::to_string(maxNlocusPairs_) + std::string(" in ") + std::string(__FUNCTION__);
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
							jaccardBlock_(tr.first, tr.second, overallPairInd, LDmatChunk);
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
				jaccardBlock_(0, maxInRAM, overallPairInd, LDmatChunk);
				overallPairInd += maxInRAM;
				for (const auto &ld : LDmatChunk){
					output << ld << " ";
				}
			}
		}
		overallPairInd = maxInRAM * nChunks;
		if (remainingPairs){
			std::vector<float> LDmatChunk(remainingPairs, 0.0);
			const size_t nLocusPairsPerThread = remainingPairs / nThreads_;
			if (nLocusPairsPerThread){
				std::vector< std::pair<size_t, size_t> > threadRanges;
				size_t locusPairInd = 0;
				for (size_t iThread = 0; iThread < nThreads_; ++iThread){
					threadRanges.emplace_back(std::pair<size_t, size_t>{locusPairInd, locusPairInd + nLocusPairsPerThread});
					locusPairInd += nLocusPairsPerThread;
				}
				threadRanges.back().second = remainingPairs;
				std::vector< std::future<void> > tasks;
				tasks.reserve(nThreads_);
				for (const auto &tr : threadRanges){
					tasks.emplace_back(
						std::async([this, &LDmatChunk, &tr, overallPairInd]{
							jaccardBlock_(tr.first, tr.second, overallPairInd, LDmatChunk);
						})
					);
					overallPairInd += nLocusPairsPerThread;
				}
			} else {
				jaccardBlock_(0, remainingPairs, overallPairInd, LDmatChunk);
			}
			for (const auto &ld : LDmatChunk){
				output << ld << " ";
			}
		}
		output.close();
	} else {
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
						jaccardBlock_(tr.first, tr.second, tr.first, LDmat);
					})
				);
			}
		} else {
			jaccardBlock_(0, LDmat.size(), 0, LDmat);
		}
		std::fstream output;
		output.open(ldFileName, std::ios::trunc | std::ios::out);
		for (const auto &ld : LDmat){
			output << ld << " ";
		}
		output.close();
	}
}

void GenoTableBin::bed2binBlk_(const std::vector<char> &bedData, const size_t &firstBedLocusInd, const size_t &lastBedLocusInd, const size_t &firstLocusInd, const size_t &bedLocusLength, const size_t &randVecLen) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	const size_t addIndv       = nIndividuals_ % 4UL;
	const size_t addBL         = bedLocusLength - static_cast<size_t>(addIndv > 0);
	size_t begByte             = firstLocusInd * binLocusSize_;
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> static_cast<uint8_t>(binLocusSize_ * byteSize_ - nIndividuals_);
	std::vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	size_t iLocus      = firstLocusInd;
	for (size_t iBedLocus = firstBedLocusInd; iBedLocus < lastBedLocusInd; ++iBedLocus){
		const size_t begInd      = iBedLocus * bedLocusLength;
		const size_t endWholeBed = begInd + addBL;
		// Fill the random byte vector
		for (auto &rv : rand){
			rv = rng_.ranInt();
		}
		std::vector<uint8_t> missMasks(binLocusSize_, 0);
		size_t iBinGeno = begByte;                       // binLocus vector index
		size_t iMissMsk = 0;
		size_t iRB      = 0;                             // randBytes index
		uint8_t bedByte = 0;
		// Two bytes of .bed code go into one byte of my binary representation
		// Therefore, work on two consecutive bytes of .bed code in the loop
		for (size_t iBed = begInd; iBed < endWholeBed ; iBed += 2){                // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
			uint8_t binByte     = 0;
			bedByte             = ~bedData[iBed];                                  // flip so that homozygous second allele (usually minor) is set to 11
			uint8_t offsetToBin = 0;                                               // move the .bed mask by this much to align with the binarized byte
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
			const size_t nextIbed = iBed + 1;
			bedByte               = ~bedData[nextIbed];
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
			// mutex may have a performance hit sometimes
			// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
			std::lock_guard<std::mutex> lk(mtx_);
			binGenotypes_[iBinGeno] = binByte;
			++iBinGeno;
			++iMissMsk;
			++iRB;
		}
		uint8_t inBedByteOffset = 0;
		uint8_t binByte         = 0;
		if (addIndv) {
			const uint8_t lastBedByte = ~bedData[endWholeBed];
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
				firstBitMask     = firstBitMask << iInd;
				binByte         |= firstBitMask;
				inBedByteOffset += 2;
				inBedByteOffset  = inBedByteOffset % 8;
			}
			// mutex may have a performance hit sometimes
			// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
			std::lock_guard<std::mutex> lk(mtx_);
			binGenotypes_[begByte + binLocusSize_ - 1] = binByte;
		}
		float aaCount = static_cast<float>( countSetBits(binGenotypes_, begByte, binLocusSize_) ) / static_cast<float>(nIndividuals_);
		if (aaCount > 0.5){ // always want the alternative to be the minor allele
			for (size_t iBL = 0; iBL < binLocusSize_; ++iBL){
				const size_t biBL = begByte + iBL;
				// mutex may have a performance hit sometimes
				// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
				std::lock_guard<std::mutex> lk(mtx_);
				binGenotypes_[biBL] = (~binGenotypes_[biBL]) & (~missMasks[iBL]);
			}
			// mutex may have a performance hit sometimes
			// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
			std::lock_guard<std::mutex> lk(mtx_);
			binGenotypes_[begByte + binLocusSize_ - 1] &= lastByteMask; // unset the remainder bits
		}
		begByte     += binLocusSize_;
		++iLocus;
	}
}

void GenoTableBin::mac2binBlk_(const std::vector<int> &macData, const size_t &startLocusInd, const size_t &endLocusInd, const size_t &randVecLen) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint8_t remainderInd       = static_cast<uint8_t>(binLocusSize_ * byteSize_ - nIndividuals_);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
	remainderInd               = byteSize_ - remainderInd;
	// Create a vector to store random bytes for stochastic heterozygote resolution
	std::vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	for (size_t iLocus = startLocusInd; iLocus < endLocusInd; ++iLocus){
		// Fill the random byte vector
		for (auto &rv : rand){
			rv = rng_.ranInt();
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
				uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);      // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= 0b10000011;                                            // mask everything in the middle
				const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
				missMasks[i0Byte]        |= (missingMask << iInByte);
				curIndiv                 &= 0b00000011;
				const uint8_t randMask    = (randBytes[iRB] >> iInByte) & oneBit_;                 // 0b00000000 or 0b00000001 with equal chance
				uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
				curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
				binByte                  |= curBitMask << iInByte;
				++iIndiv;
			}
			// mutex may have a performance hit sometimes
			// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
			std::lock_guard<std::mutex> lk(mtx_);
			binGenotypes_[iByte] = binByte;
			++i0Byte;
			++iRB;
		}
		// now deal with the last byte in the individual
		uint8_t lastBinByte = 0;
		for (uint8_t iRem = 0; iRem < remainderInd; ++iRem){
			uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);          // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= 0b10000011;                                                // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
			missMasks.back()         |= (missingMask << iRem);
			curIndiv                 &= 0b00000011;                                                
			const uint8_t randMask    = (randBytes[binLocusSize_ - 1] >> iRem) & oneBit_;          // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                   // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                              // zero it out if missing value is set

			lastBinByte              |= curBitMask << iRem;
			++iIndiv;
		}
		{
			// mutex may have a performance hit sometimes
			// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
			std::lock_guard<std::mutex> lk(mtx_);
			binGenotypes_[begByte + binLocusSize_ - 1] = lastBinByte;
		}
		float maf = static_cast<float>( countSetBits(binGenotypes_, begByte, binLocusSize_) ) / static_cast<float>(nIndividuals_);
		if (maf > 0.5){ // always want the alternative to be the minor allele
			i0Byte = 0;
			for (size_t i = begByte; i < begByte + binLocusSize_; ++i){
				// mutex may have a performance hit sometimes
				// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
				std::lock_guard<std::mutex> lk(mtx_);
				binGenotypes_[i] = (~binGenotypes_[i]) & (~missMasks[i0Byte]);
				++i0Byte;
			}
			// mutex may have a performance hit sometimes
			// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
			std::lock_guard<std::mutex> lk(mtx_);
			binGenotypes_[begByte + binLocusSize_ - 1] &= lastByteMask; // unset the remainder bits
		}
	}
}

void GenoTableBin::jaccardBlock_(const size_t &blockStartVec, const size_t &blockEndVec, const size_t &blockStartAll, std::vector<float> &jaccardVec) const {
	const size_t nnLoci = nLoci_ * (nLoci_ - 1) / 2 - 1; // overflow checked in the calling function; do this here for encapsulation, may move to the calling function later
	std::vector<uint8_t> locus(binLocusSize_);
	size_t curJacMatInd = blockStartAll;
	for (size_t iVecInd = blockStartVec; iVecInd < blockEndVec; ++iVecInd){
		// compute row and column indexes from the vectorized lower triangle index
		// got these expressions by combining various web sources and verifying
		const size_t kp     = nnLoci - curJacMatInd;
		const size_t p      = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kp) ) ) - 1) / 2;
		const size_t row    = nLoci_ - 2 - p;
		const size_t col    = nLoci_ - (kp - p * (p + 1) / 2) - 1;
		const size_t rowBin = row * binLocusSize_;
		const size_t colBin = col * binLocusSize_;
		for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc){
			locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] & binGenotypes_[colBin + iBinLoc];
		}
		const uint32_t uni = countSetBits(locus);
		for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc){
			locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] | binGenotypes_[colBin + iBinLoc];
		}
		const uint32_t isect = countSetBits(locus);
		// mutex may have a performance hit sometimes
		// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
		std::lock_guard<std::mutex> lk(mtx_);
		jaccardVec[iVecInd]  = static_cast<float>(uni) / static_cast<float>(isect);
		++curJacMatInd;
	}
}

// GenoTableHash methods
constexpr std::array<char, 3> GenoTableHash::magicBytes_ = {0x6c, 0x1b, 0x01};                   // Leading bytes for .bed files 
constexpr uint8_t  GenoTableHash::oneBit_                = 0b00000001;                           // One set bit for masking 
constexpr uint8_t  GenoTableHash::byteSize_              = 8;                                    // Size of one byte in bits 
constexpr uint8_t  GenoTableHash::llWordSize_            = 8;                                    // 64 bit word size in bytes 
constexpr size_t   GenoTableHash::maxPairs_              = 6074000999UL;                         // approximate maximum number that does not overflow with n*(n-1)/2
constexpr size_t   GenoTableHash::nblocks_               = sizeof(size_t) / 4;                   // MurMurHash number of blocks 
constexpr uint32_t GenoTableHash::mmhKeyLen_             = sizeof(size_t);                       // MurMurHash key length 
constexpr uint16_t GenoTableHash::emptyBinToken_         = std::numeric_limits<uint16_t>::max(); // Value corresponding to an empty token 
constexpr uint32_t GenoTableHash::c1_                    = 0xcc9e2d51;                           // MurMurHash c1 constant 
constexpr uint32_t GenoTableHash::c2_                    = 0x1b873593;                           // MurMurHash c2 constant 

// Constructors
GenoTableHash::GenoTableHash(const std::string &inputFileName, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, const std::string &logFileName) : nIndividuals_{nIndividuals}, kSketches_{kSketches}, nLoci_{0}, nThreads_{nThreads}, logFileName_{logFileName} {
	std::stringstream logStream;
	const time_t t = time(nullptr);
	logStream << std::put_time(localtime(&t), "%b %e %H:%M %Z");
	logMessages_ = "Genotype hashing from a .bed file started on " + logStream.str() + "\n";
	logStream.clear();
	if (nIndividuals <= 1){
		logMessages_ += "ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string(__FUNCTION__);
	}
	if (kSketches_ < 3){
		logMessages_ += "ERROR: sketch size (" + std::to_string(kSketches_) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: sketch size must be at least three in ") + std::string(__FUNCTION__);
	}
	sketchSize_ = nIndividuals_ / kSketches + static_cast<bool>(nIndividuals_ % kSketches);
	if (sketchSize_ >= emptyBinToken_){
		logMessages_ += "ERROR: sketch size (" + std::to_string(sketchSize_) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: Number of sketches (") + std::to_string(kSketches_) + std::string(") implies sketch size (") +
			std::to_string(sketchSize_) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string(__FUNCTION__);
	}
	const size_t nBedBytes = nIndividuals_ / 4 + static_cast<bool>(nIndividuals_ % 4);
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
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string(__FUNCTION__);
	}
	const int32_t endPosition = inStr.tellg();
	if ( endPosition < magicBytes_.size() ){
		logMessages_ += "ERROR: no loci in the input .bed file " + inputFileName + "; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string(__FUNCTION__);
	}
	const size_t N = static_cast<uint64_t>(endPosition) - magicBytes_.size();
	inStr.close();
	nLoci_ = N / nBedBytes;
	logMessages_ += "Number of individuals: " + std::to_string(nIndividuals_) + "\n";
	logMessages_ += "Number of loci: " + std::to_string(nLoci_) + "\n";
	logMessages_ += "Hash size: " + std::to_string(kSketches_) + "\n";
	// Calculate the actual sketch number based on the realized sketch size
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	locusSize_ = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	char magicBuf[magicBytes_.size()]{};
	inStr.read( magicBuf, magicBytes_.size() );
	if (magicBuf[0] != magicBytes_[0]){
		logMessages_ += "ERROR: file " + inputFileName + " does not appear to be in .bed format; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: first magic byte in input .bed file is not the expected value in ") + std::string(__FUNCTION__);
	} else if (magicBuf[1] != magicBytes_[1]){
		logMessages_ += "ERROR: file " + inputFileName + " does not appear to be in .bed format; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: second magic byte in input .bed file is not the expected value in ") + std::string(__FUNCTION__);
	} else if (magicBuf[2] != magicBytes_[2]){
		logMessages_ += "ERROR: file " + inputFileName + " does not appear to be in .bed format; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: third magic byte in input .bed file is not the expected value in ") + std::string(__FUNCTION__);
	}
	// Generate the binary genotype table while reading the .bed file
	const size_t nBedBytesPerLocus = nIndividuals_ / 4 + static_cast<bool>(nIndividuals_ % 4);
	const size_t ranVecSize        = nBedBytes / llWordSize_ + static_cast<bool>(nBedBytes % llWordSize_);
	const size_t ramSize           = getAvailableRAM() / 2UL;                               // measuring here, after all the major allocations; use half to leave resources for other operations
	size_t nBedLociToRead          = ramSize / nBedBytesPerLocus;                           // number of .bed loci to read at a time
	nBedLociToRead                 = (nBedLociToRead < nLoci_ ? nBedLociToRead : nLoci_);
	const size_t remainingLoci     = nLoci_ % nBedLociToRead;
	const size_t remainingBytes    = remainingLoci * nBedBytesPerLocus;
	const size_t nChunks           = nLoci_ / nBedLociToRead;
	const size_t nBedBytesToRead   = nBedLociToRead * nBedBytesPerLocus;
	size_t nLociPerThread          = nBedLociToRead / nThreads_;
	std::vector<char> bedChunkToRead(nBedBytesToRead, 0);
	logMessages_ += "RAM available for reading the .bed file: " + std::to_string(ramSize) + " bytes\n";
	logMessages_ += ".bed file will be read in " + std::to_string(nChunks) + " chunk(s)\n";
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{rng_.shuffleUint(nIndividuals_)};
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
			inStr.read(bedChunkToRead.data(), nBedBytesToRead);
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &tr : threadRanges){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, &tr, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &seeds]{
						bed2ophBlk_(bedChunkToRead, tr.first, tr.second, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, seeds);
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
			inStr.read(bedChunkToRead.data(), nBedBytesToRead);
			std::vector< std::future<void> > tasks;
			tasks.reserve(nBedLociToRead);
			for (size_t iBedLocus = 0; iBedLocus < nBedLociToRead; ++iBedLocus){
				tasks.emplace_back(
					std::async([this, &bedChunkToRead, iBedLocus, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &seeds]{
							bed2ophBlk_(bedChunkToRead, iBedLocus, iBedLocus + 1, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, seeds);
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
		inStr.read(bedChunkToRead.data(), remainingBytes);
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
					std::async([this, &bedChunkToRead, &tr, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &seeds]{
						bed2ophBlk_(bedChunkToRead, tr.first, tr.second, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, seeds);
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
					std::async([this, &bedChunkToRead, iBedLocus, locusInd, nBedBytesPerLocus, ranVecSize, &ranInts, &seeds]{
							bed2ophBlk_(bedChunkToRead, iBedLocus, iBedLocus + 1, locusInd, nBedBytesPerLocus, ranVecSize, ranInts, seeds);
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

GenoTableHash::GenoTableHash(const std::vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, const std::string &logFileName) : nIndividuals_{nIndividuals}, kSketches_{kSketches}, nLoci_{maCounts.size() / nIndividuals}, logFileName_{logFileName} {
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
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string(__FUNCTION__);
	}
	if (maCounts.size() % nIndividuals){
		logMessages_ += "ERROR: minor allele vector size (" + std::to_string( maCounts.size() ) + ") is not evenly divisible by the number of individuals (" + std::to_string(nIndividuals_) + "); aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(nIndividuals) + std::string( ") in ") + std::string(__FUNCTION__);
	}
	if ( maCounts.empty() ){
		logMessages_ += "ERROR: minor allele count vector is empty; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string(__FUNCTION__);
	}
	if (kSketches_ < 3){
		logMessages_ += "ERROR: sketch size (" + std::to_string(kSketches_) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: sketch size must be at least three in ") + std::string(__FUNCTION__);
	}
	if (nThreads_ == 0){
		nThreads_ = 1;
	} else if ( nThreads_ > std::thread::hardware_concurrency() ){
		nThreads_ = std::thread::hardware_concurrency();
	}
	sketchSize_ = nIndividuals_ / kSketches_ + static_cast<bool>(nIndividuals_ % kSketches_);
	if (sketchSize_ >= emptyBinToken_){
		logMessages_ += "ERROR: sketch size (" + std::to_string(sketchSize_) + ") is too small; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: Number of sketches (") + std::to_string(kSketches_) + std::string(") implies sketch size (") +
			std::to_string(sketchSize_) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string(__FUNCTION__);
	}
	locusSize_              = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	const size_t ranVecSize = locusSize_ / llWordSize_ + static_cast<bool>(locusSize_ % llWordSize_);
	// Calculate the actual sketch number based on the realized sketch size
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{rng_.shuffleUint(nIndividuals_)};
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
	if (this != &in){
		sketches_     = move(in.sketches_);
		nIndividuals_ = in.nIndividuals_;
		kSketches_    = in.kSketches_;
		nLoci_        = in.nLoci_;
		sketchSize_   = in.sketchSize_;
		nThreads_     = in.nThreads_;
		locusSize_    = in.locusSize_;
		logFileName_  = move(in.logFileName_);
		logMessages_  = move(in.logMessages_);

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
			 std::to_string(maxPairs_) + std::string(" in ") + std::string(__FUNCTION__);
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
			const size_t nLocusPairsPerThread = remainingPairs / nThreads_;
			if (nLocusPairsPerThread){
				std::vector< std::pair<size_t, size_t> > threadRanges;
				size_t locusPairInd = 0;
				for (size_t iThread = 0; iThread < nThreads_; ++iThread){
					threadRanges.emplace_back(std::pair<size_t, size_t>{locusPairInd, locusPairInd + nLocusPairsPerThread});
					locusPairInd += nLocusPairsPerThread;
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
					overallPairInd += nLocusPairsPerThread;
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
void GenoTableHash::testSimHash(const size_t &kSketchSubset, const std::string &outFileName) const {
	std::fstream output;
	output.open(outFileName, std::ios::trunc | std::ios::out);
	for (size_t iSketch = 0; iSketch < kSketches_; ++iSketch){
		for (size_t iLocus = 0; iLocus < nLoci_; ++iLocus){
			output << sketches_[iLocus * kSketches_ + iSketch] << " ";
		}
		output << "\n";
	}

	const uint32_t seed = static_cast<uint32_t>( rng_.ranInt() );
	for (size_t iLocus = 0; iLocus < nLoci_; ++iLocus){
		output << simHash_(iLocus * kSketches_, kSketches_, seed) << "\n";
	}

}
std::vector< std::vector<size_t> > GenoTableHash::makeLDgroups(const uint16_t &hammingCutoff, const size_t &kSketchSubset, const size_t &lookBackNumber) const {
	if (kSketchSubset == 0){
		logMessages_ += "ERROR: cannot have zero sketches to simHash; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: number of OPH sketches to simHash cannot be 0 in ") + std::string(__FUNCTION__);
	} else if (kSketchSubset > kSketches_){
		logMessages_ += "ERROR: subset of sketches (" + std::to_string(kSketchSubset) + ") larger than the total; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: number of OPH sketches to simHash cannot exceed the total number of sketches per locus in ") + std::string(__FUNCTION__);
	} else if (lookBackNumber == 0){
		logMessages_ += "ERROR: look back number must not be zero; aborting\n";
		std::fstream outLog;
		outLog.open(logFileName_, std::ios::out | std::ios::trunc);
		outLog << logMessages_;
		outLog.close();
		throw std::string("ERROR: zero look back number in ") + std::string(__FUNCTION__);
	}

	logMessages_ += "Grouping loci\n";
	logMessages_ += "Hamming distance cut-off: " + std::to_string(hammingCutoff) + "\n";
	logMessages_ += "OPH subset size: " + std::to_string(kSketchSubset) + "\n";
	logMessages_ += "Number of preceding groups to consider: " + std::to_string(lookBackNumber) + "\n";

	const uint32_t seed = static_cast<uint32_t>( rng_.ranInt() );
	std::vector< std::vector<size_t> > ldGroup;                         // each element is a vector of indexes into sketches_
	std::vector<uint16_t> activeHashes(lookBackNumber, 0);              // hashes that are under consideration in a simplified ring buffer

	// Start by filling the active hash buffer
	//
	// initialize with the first locus
	activeHashes.back() = simHash_(0, kSketchSubset, seed);
	auto latestHashIt   = activeHashes.rbegin();                   // points to the latest hash under consideration
	size_t locusInd     = 1;                                       // current locus index (starts from the second locus)
	ldGroup.emplace_back(std::vector<size_t>{0});
	// loop through loci until buffer is full or no more loci left
	while ( ( latestHashIt != activeHashes.rend() ) && (locusInd < nLoci_) ){
		const uint16_t curHash = simHash_(locusInd * kSketches_, kSketchSubset, seed);
		auto fwdHashIt         = activeHashes.cbegin() + distance( latestHashIt, activeHashes.rend() ) - 1;
		auto ldgIt             = ldGroup.rbegin();   // latest group is at the end
		bool noCollision       = true;
		while ( fwdHashIt != activeHashes.cend() ){
			if (hammingDistance_(*fwdHashIt, curHash) <= hammingCutoff){
				noCollision = false;
				break;
			}
			++fwdHashIt;
			++ldgIt;
		}
		if (noCollision){
			++latestHashIt;
			if ( latestHashIt == activeHashes.rend() ){
				ldGroup.emplace_back(std::vector<size_t>{locusInd});
				activeHashes[0] = curHash;
				++locusInd;
				break;
			} else {
				*latestHashIt = curHash;
				ldGroup.emplace_back(std::vector<size_t>{locusInd});
				++locusInd;
			}
		} else {
			ldgIt->push_back(locusInd);
			++locusInd;
		}
	}
	// If there are loci left after filling the buffer, continue.
	while (locusInd < nLoci_){
		const uint16_t curHash = simHash_(locusInd * kSketches_, kSketchSubset, seed);
		size_t latestHashInd   = 0;
		auto ldgIt             = ldGroup.rbegin();   // latest group is at the end
		bool noCollision       = true;
		for (size_t iTestHash = 0; iTestHash < lookBackNumber; ++iTestHash){
			if (hammingDistance_(activeHashes[(latestHashInd + iTestHash) % lookBackNumber], curHash) <= hammingCutoff){
				noCollision = false;
				break;
			}
			++ldgIt;
		}
		if (noCollision){
			++latestHashInd;
			latestHashInd %= lookBackNumber;
			ldGroup.emplace_back(std::vector<size_t>{locusInd});
			++locusInd;
		} else {
			ldgIt->push_back(locusInd);
			++locusInd;
		}
	}

	return ldGroup;
}

void GenoTableHash::ldInGroups(const uint16_t &hammingCutoff, const size_t &kSketchSubset, const size_t &lookBackNumber, const size_t &smallestGrpSize, const std::string &outFileName) const {
	std::vector< std::vector<size_t> > ldGroup = this->makeLDgroups(hammingCutoff, kSketchSubset, lookBackNumber);
	size_t totNpairs = 0;
	for (const auto &ldg : ldGroup){
		if (ldg.size() >= smallestGrpSize){
			size_t nLowTri = ldg.size() * ( ldg.size() - 1 ) / 2;
			if (totNpairs >= std::numeric_limits<size_t>::max() - nLowTri){
				logMessages_ += "ERROR: number of locus pairs in grouped LD exceeds the 64-bit unsigned integer limit; aborting\n";
				std::fstream outLog;
				outLog.open(logFileName_, std::ios::out | std::ios::trunc);
				outLog << logMessages_;
				outLog.close();
				throw std::string("ERROR: Number of locus pairs exceeds 64-bit unsigned integer limit in ") + std::string(__FUNCTION__);
			}
			totNpairs += nLowTri;
		}
	}
	const size_t maxPairsInRAM = getAvailableRAM() / ( 2UL * ( sizeof(float) + sizeof(std::pair<size_t, size_t>) + sizeof(size_t) ) );      // use half to leave resources for other operations

	logMessages_ += "Estimating LD in groups\n";
	logMessages_ += "Smallest group size: " + std::to_string(smallestGrpSize) + "\n";
	logMessages_ += "Number of pairs considered: " + std::to_string(totNpairs) + "\n";
	logMessages_ += "Maximum number fitting in RAM: " + std::to_string(maxPairsInRAM) + "; ";
	if (totNpairs > maxPairsInRAM){
		logMessages_                += "calculating in chunks\n";
		const size_t nRAMchunks      = totNpairs / maxPairsInRAM;
		const size_t nRemainingPairs = totNpairs % maxPairsInRAM;
		size_t iGroup                = 0;
		size_t iLocus                = 0;
		size_t jLocus                = 1;

		std::fstream out;
		out.open(outFileName, std::ios::out | std::ios::trunc);
		out << "groupID\tlocus1\tlocus2\tjaccLD\n";
		size_t iKeptGroup = 0;

		for (size_t iChunk = 0; iChunk < nRAMchunks; ++iChunk){
			std::vector<float> hashJacGroups(maxPairsInRAM, 0.0);
			std::vector< std::pair<size_t, size_t> > idxPairs;
			std::vector<size_t> groupID;
			idxPairs.reserve(maxPairsInRAM);
			size_t iPair = 0;
			for (; iGroup < ldGroup.size(); ++iGroup) {
				if (ldGroup[iGroup].size() >= smallestGrpSize){
					for (; iLocus < ldGroup[iGroup].size() - 1; ++iLocus){
						for (; jLocus < ldGroup[iGroup].size(); ++jLocus){
							idxPairs.emplace_back(std::pair<size_t, size_t>{ldGroup[iGroup][iLocus], ldGroup[iGroup][jLocus]});
							groupID.push_back(iKeptGroup);
							++iPair;
							if (iPair == maxPairsInRAM){
								++jLocus;
								if ( jLocus == ldGroup[iGroup].size() ){
									if (iLocus < ldGroup[iGroup].size() - 1){
										++iLocus;
										jLocus = iLocus + 1;
									} else {
										++iGroup;
										++iKeptGroup;
										iLocus = 0;
										jLocus = 1;
									}
								}
								goto stopAllLoops;
							}
						}
						jLocus = iLocus + 2;
					}
					++iKeptGroup;
					iLocus = 0;
					jLocus = 1;
				}
			}
		stopAllLoops:
			const size_t nPairsPerThread = maxPairsInRAM / nThreads_;
			if (nPairsPerThread){
				std::vector< std::pair<size_t, size_t> > threadRanges;
				threadRanges.reserve(nThreads_);
				size_t pairInd = 0;
				for (size_t iThread = 0; iThread < nThreads_; ++iThread){
					threadRanges.emplace_back(std::pair<size_t, size_t>{pairInd, pairInd + nPairsPerThread});
					pairInd += nPairsPerThread;
				}
				threadRanges.back().second = maxPairsInRAM;
				std::vector< std::future<void> > tasks;
				tasks.reserve(nThreads_);
				for (const auto &tr : threadRanges){
					tasks.emplace_back(
						std::async([this, &tr, &idxPairs, &hashJacGroups]{
							hashJacBlock_(tr.first, tr.second, idxPairs, hashJacGroups);
						})
					);
				}
				for (const auto &th : tasks){
					th.wait();
				}
			} else {
				std::vector< std::future<void> > tasks;
				tasks.reserve(maxPairsInRAM);
				for (size_t iPair = 0; iPair < maxPairsInRAM; ++iPair){
					tasks.emplace_back(
						std::async([this, iPair, &idxPairs, &hashJacGroups]{
							hashJacBlock_(iPair, iPair + 1, idxPairs, hashJacGroups);
						})
					);
				}
				for (const auto &th : tasks){
					th.wait();
				}
			}
			// Save the results
			for (size_t iPair = 0; iPair < idxPairs.size(); ++iPair){
				out << groupID[iPair] + 1 << "\t" << idxPairs[iPair].first + 1 << "\t" << idxPairs[iPair].second + 1 << "\t" << hashJacGroups[iPair] << "\n";
			}
		}
		if (nRemainingPairs){
			std::vector<float> hashJacGroups(nRemainingPairs, 0.0);
			std::vector< std::pair<size_t, size_t> > idxPairs;
			std::vector<size_t> groupID;
			idxPairs.reserve(nRemainingPairs);
			for (; iGroup < ldGroup.size(); ++iGroup) {
				if (ldGroup[iGroup].size() >= smallestGrpSize){
					for (; iLocus < ldGroup[iGroup].size() - 1; ++iLocus){
						for (; jLocus < ldGroup[iGroup].size(); ++jLocus){
							idxPairs.emplace_back(std::pair<size_t, size_t>{ldGroup[iGroup][iLocus], ldGroup[iGroup][jLocus]});
							groupID.push_back(iKeptGroup);
						}
						jLocus = iLocus + 2;
					}
					++iKeptGroup;
					iLocus = 0;
					jLocus = 1;
				}
			}
			const size_t nPairsPerThread = nRemainingPairs / nThreads_;
			if (nPairsPerThread){
				std::vector< std::pair<size_t, size_t> > threadRanges;
				threadRanges.reserve(nThreads_);
				size_t pairInd = 0;
				for (size_t iThread = 0; iThread < nThreads_; ++iThread){
					threadRanges.emplace_back(std::pair<size_t, size_t>{pairInd, pairInd + nPairsPerThread});
					pairInd += nPairsPerThread;
				}
				threadRanges.back().second = nRemainingPairs;
				std::vector< std::future<void> > tasks;
				tasks.reserve(nThreads_);
				for (const auto &tr : threadRanges){
					tasks.emplace_back(
						std::async([this, &tr, &idxPairs, &hashJacGroups]{
							hashJacBlock_(tr.first, tr.second, idxPairs, hashJacGroups);
						})
					);
				}
				for (const auto &th : tasks){
					th.wait();
				}
			} else {
				std::vector< std::future<void> > tasks;
				tasks.reserve(nRemainingPairs);
				for (size_t iPair = 0; iPair < nRemainingPairs; ++iPair){
					tasks.emplace_back(
						std::async([this, iPair, &idxPairs, &hashJacGroups]{
							hashJacBlock_(iPair, iPair + 1, idxPairs, hashJacGroups);
						})
					);
				}
				for (const auto &th : tasks){
					th.wait();
				}
			}
			// Save the results
			for (size_t iPair = 0; iPair < idxPairs.size(); ++iPair){
				out << groupID[iPair] + 1 << "\t" << idxPairs[iPair].first + 1 << "\t" << idxPairs[iPair].second + 1 << "\t" << hashJacGroups[iPair] << "\n";
			}
		}
		out.close();
	} else {
		logMessages_ += "calculating in one go\n";
		// estimate Jaccard similarities within groups
		std::vector<float> hashJacGroups(totNpairs, 0.0);
		std::vector< std::pair<size_t, size_t> > idxPairs;
		std::vector<size_t> groupID;
		idxPairs.reserve(totNpairs);
		size_t iKeptGroup = 0;
		for (auto &ldg : ldGroup){
			if (ldg.size() >= smallestGrpSize){
				for (size_t iLocus = 0; iLocus < ldg.size() - 1; ++iLocus){
					for (size_t jLocus = iLocus + 1; jLocus < ldg.size(); ++jLocus){
						idxPairs.emplace_back(std::pair<size_t, size_t>{ldg[iLocus], ldg[jLocus]});
						groupID.push_back(iKeptGroup);
					}
				}
				++iKeptGroup;
				ldg.clear();
			}
		}
		const size_t nPairsPerThread = totNpairs / nThreads_;
		if (nPairsPerThread){
			std::vector< std::pair<size_t, size_t> > threadRanges;
			threadRanges.reserve(nThreads_);
			size_t pairInd = 0;
			for (size_t iThread = 0; iThread < nThreads_; ++iThread){
				threadRanges.emplace_back(std::pair<size_t, size_t>{pairInd, pairInd + nPairsPerThread});
				pairInd += nPairsPerThread;
			}
			threadRanges.back().second = totNpairs;
			std::vector< std::future<void> > tasks;
			tasks.reserve(nThreads_);
			for (const auto &tr : threadRanges){
				tasks.emplace_back(
					std::async([this, &tr, &idxPairs, &hashJacGroups]{
						hashJacBlock_(tr.first, tr.second, idxPairs, hashJacGroups);
					})
				);
			}
			for (const auto &th : tasks){
				th.wait();
			}
		} else {
			std::vector< std::future<void> > tasks;
			tasks.reserve(totNpairs);
			for (size_t iPair = 0; iPair < totNpairs; ++iPair){
				tasks.emplace_back(
					std::async([this, iPair, &idxPairs, &hashJacGroups]{
						hashJacBlock_(iPair, iPair + 1, idxPairs, hashJacGroups);
					})
				);
			}
			for (const auto &th : tasks){
				th.wait();
			}
		}
		// Save the results
		std::fstream out;
		out.open(outFileName, std::ios::out | std::ios::trunc);
		out << "groupID\tlocus1\tlocus2\tjaccLD\n";
		for (size_t iPair = 0; iPair < idxPairs.size(); ++iPair){
			out << groupID[iPair] + 1 << "\t" << idxPairs[iPair].first + 1 << "\t" << idxPairs[iPair].second + 1 << "\t" << hashJacGroups[iPair] << "\n";
		}
		out.close();
	}
}

void GenoTableHash::saveLogFile() const {
	std::fstream outLog;
	outLog.open(logFileName_, std::ios::out | std::ios::trunc);
	outLog << logMessages_;
	outLog.close();
}

void GenoTableHash::locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds, std::vector<uint8_t> &binLocus){
	// Start with a permutation to make OPH
	size_t iIndiv = nIndividuals_ - 1UL; // safe b/c nIndividuals_ > 1 is checked at construction
	for (const auto &ri : permutation){
		uint16_t firstIdx  = iIndiv % byteSize_;
		size_t firstByte   = (iIndiv / byteSize_);
		uint16_t secondIdx = ri % byteSize_;
		size_t secondByte  = (ri / byteSize_);
		uint16_t diff      = byteSize_ * (firstByte != secondByte); // will be 0 if the same byte is being accessed; then need to swap bits within byte

		// swapping bits within a two-byte variable
		// using the method in https://graphics.stanford.edu/~seander/bithacks.html#SwappingBitsXOR
		// if the same byte is being accessed, secondIdx is not shifted to the second byte
		// This may be affected by endianness (did not test)
		uint16_t twoBytes  = (static_cast<uint16_t>(binLocus[secondByte]) << 8) | ( static_cast<uint16_t>(binLocus[firstByte]) );
		secondIdx         += diff;
		uint16_t x         = ( (twoBytes >> firstIdx) ^ (twoBytes >> secondIdx) ) & 1;
		twoBytes          ^= ( (x << firstIdx) | (x << secondIdx) );

		memcpy( binLocus.data() + firstByte, &twoBytes, sizeof(uint8_t) );
		twoBytes = twoBytes >> diff;
		memcpy( binLocus.data() + secondByte, &twoBytes, sizeof(uint8_t) );
		--iIndiv;
	}
	// Now make the sketches
	std::vector<size_t> filledIndexes;                     // indexes of the non-empty sketches
	size_t iByte     = 0;
	size_t colEnd    = iByte + locusSize_;
	size_t sketchBeg = locusInd * kSketches_;
	size_t iSeed     = 0;                             // index into the seed vector
	uint8_t iInByte  = 0;
	// A possible optimization is to test a whole byte for 0
	// Will test later
	for (size_t iSketch = 0; iSketch < kSketches_; ++iSketch){
		uint16_t firstSetBitPos = 0;
		while ( (iByte != colEnd) && ( ( (oneBit_ << iInByte) & binLocus[iByte] ) == 0 ) &&
				(firstSetBitPos < sketchSize_) ){
			++iInByte;
			// these are instead of an if statement
			iByte  += iInByte == byteSize_;
			iInByte = iInByte % byteSize_;
			++firstSetBitPos;
		}
		if ( (iByte < colEnd) && (firstSetBitPos < sketchSize_) ){
			filledIndexes.push_back(iSketch);
			{
				std::lock_guard<std::mutex> lk(mtx_);
				sketches_[sketchBeg + iSketch] = firstSetBitPos;
			}

			uint16_t remainder = sketchSize_ - firstSetBitPos;
			uint16_t inByteMod = remainder % byteSize_;
			uint16_t inByteSum = iInByte + inByteMod;

			iByte  += remainder / byteSize_ + inByteSum / byteSize_;
			iInByte = inByteSum % byteSize_;
		}
	}
	if (filledIndexes.size() == 1){
		for (size_t iSk = 0; iSk < kSketches_; ++iSk){ // this will overwrite the one assigned sketch, but the wasted operation should be swamped by the rest
			std::lock_guard<std::mutex> lk(mtx_);
			sketches_[sketchBeg + iSk] = sketches_[filledIndexes[0] + sketchBeg];
		}
	} else if (filledIndexes.size() != kSketches_){
		if ( filledIndexes.empty() ){ // in the case where the whole locus is monomorphic, pick a random index as filled
			filledIndexes.push_back( rng_.sampleInt(kSketches_) );
		}
		size_t emptyCount = kSketches_ - filledIndexes.size();
		while (emptyCount){
			for (const auto &f : filledIndexes){
				uint32_t newIdx = murMurHash_(f, seeds[iSeed]) % kSketches_ + sketchBeg;
				std::lock_guard<std::mutex> lk(mtx_);
				if ( sketches_[newIdx] == emptyBinToken_ ){
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
					const size_t &bedLocusLength, const size_t &randVecLen, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	const size_t addIndv       = nIndividuals_ % 4UL;
	const size_t addBL         = bedLocusLength - static_cast<size_t>(addIndv > 0);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	// Fill the random byte vector
	std::vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	size_t iLocus      = firstLocusInd;
	for (size_t iBedLocus = firstBedLocusInd; iBedLocus < lastBedLocusInd; ++iBedLocus){
		const size_t begInd      = iBedLocus * bedLocusLength;
		const size_t endWholeBed = begInd + addBL;
		for (auto &rv : rand){
			rv = rng_.ranInt();
		}
		std::vector<uint8_t> binLocus(locusSize_, 0);
		std::vector<uint8_t> missMasks(locusSize_, 0);
		size_t iBinGeno = 0;                       // binLocus vector index
		uint8_t bedByte = 0;
		size_t iRB      = 0;                       // random byte index
		// Two bytes of .bed code go into one byte of my binary representation
		// Therefore, work on two consecutive bytes of .bed code in the loop
		for (size_t iBed = begInd; iBed < endWholeBed; iBed += 2){                 // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
			bedByte             = ~bedData[iBed];                                  // flip so that homozygous second allele (usually minor) is set to 11
			uint8_t offsetToBin = 0;                                               // move the .bed mask by this much to align with the binarized byte
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
			bedByte               = ~bedData[nextIbed];
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
			const uint8_t lastBedByte = ~bedData[endWholeBed];
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
				firstBitMask     = firstBitMask << iInd;
				binLocus.back() |= firstBitMask;
				inBedByteOffset += 2;
				inBedByteOffset  = inBedByteOffset % 8;
			}
		}
		float aaCount = static_cast<float>( countSetBits(binLocus) ) / static_cast<float>(nIndividuals_);
		if (aaCount > 0.5){ // always want the alternative to be the minor allele
			for (size_t iBL = 0; iBL < locusSize_; ++iBL){
				binLocus[iBL] = (~binLocus[iBL]) & (~missMasks[iBL]);
			}
			binLocus.back() &= lastByteMask; // unset the remainder bits
		}
		locusOPH_(iLocus, permutation, seeds, binLocus);
		++iLocus;
	}
}

void GenoTableHash::mac2ophBlk_(const std::vector<int> &macData, const size_t &startLocusInd, const size_t &endLocusInd, const size_t &randVecLen, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	uint8_t remainderInd       = static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
	remainderInd               = byteSize_ - remainderInd;
	// Create a vector to store random bytes for stochastic heterozygote resolution
	std::vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	for (size_t iLocus = startLocusInd; iLocus < endLocusInd; ++iLocus){
		// Fill the random byte vector
		for (auto &rv : rand){
			rv = rng_.ranInt();
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

	for(size_t i = 0; i < nblocks_; ++i){
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

uint16_t GenoTableHash::murMurHash_(const uint16_t &key, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	uint32_t block = static_cast<uint32_t>(key);

	uint32_t k1 = block;

	k1 *= c1_;
	k1 = (k1 << 15) | (k1 >> 17);
	k1 *= c2_;

	hash ^= k1;
	hash  = (hash << 13) | (hash >> 19);
	hash  = hash * 5 + 0xe6546b64;

	// there is no tail since the input is fixed (at eight bytes typically)
	// finalize
	hash ^= mmhKeyLen_;
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;
	
	auto hashV = reinterpret_cast<const uint16_t *>(&hash);

	return hashV[0] ^ hashV[1];
}

uint32_t GenoTableHash::murMurHash_(const size_t &startInd, const size_t &nElements, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	auto blocks    = reinterpret_cast<const uint32_t *>(sketches_.data() + startInd);
	size_t nBlocks = nElements / 2; // each sketch is 16 bits; this implies that only even number of sketches is considered

	for(size_t i = 0; i < nBlocks; ++i){
		uint32_t k1 = blocks[i];

		k1 *= c1_;
		k1 = (k1 << 15) | (k1 >> 17);
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

uint32_t GenoTableHash::simHash_(const size_t &startInd, const size_t &kSketches, const uint32_t &seed) const {
	uint32_t hash        = 0;
	const uint32_t one   = 1;
	const size_t twoByte = 2 * byteSize_;
	std::array<int32_t, twoByte> v{};
	for (size_t iSketch = startInd; iSketch < startInd + kSketches; ++iSketch){
		const uint32_t skHash = murMurHash_(sketches_[iSketch], seed);
		for (uint32_t j = 0; j < twoByte; ++j){
			v[j] += -1 + 2 * static_cast<int32_t>( one & (skHash >> j) );
		}
	}
	for (uint32_t i = 0; i < twoByte; ++i){
		hash |= (static_cast<uint32_t>(v[i] > 0) << i);
	}
	return hash;
}

void GenoTableHash::hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, const size_t &blockStartAll, std::vector<float> &hashJacVec) const {
	const size_t nnLoci = nLoci_ * (nLoci_ - 1) / 2 - 1; // overflow checked in the calling function; do this here for encapsulation, may move to the calling function later
	size_t curJacMatInd = blockStartAll;
	for (size_t iVecInd = blockStartVec; iVecInd < blockEndVec; ++iVecInd){
		// compute row and column indexes from the vectorized by column lower triangle index
		// got these expressions by combining various web sources and verifying
		const size_t kp    = nnLoci - curJacMatInd;
		const size_t p     = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kp) ) ) - 1) / 2;
		const size_t row   = nLoci_ - 2 - p;
		const size_t col   = nLoci_ - (kp - p * (p + 1) / 2) - 1;
		const size_t rowSk = row * kSketches_;
		const size_t colSk = col * kSketches_;
		float simVal       = 0.0;
		for (size_t iSk = 0; iSk < kSketches_; ++iSk){
			simVal += static_cast<float>(sketches_[rowSk + iSk] == sketches_[colSk + iSk]);
		}
		// mutex may have a performance hit sometimes
		// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
		std::lock_guard<std::mutex> lk(mtx_);
		hashJacVec[iVecInd] = simVal / static_cast<float>(kSketches_);
		++curJacMatInd;
	}
}

void GenoTableHash::hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, const std::vector< std::pair<size_t, size_t> > &idxVector, std::vector<float> &hashJacVec) const {
	for (size_t iBlock = blockStartVec; iBlock < blockEndVec; ++iBlock){
		float simVal = 0.0;
		for (size_t iSk = 0; iSk < kSketches_; ++iSk){
			simVal += static_cast<float>(sketches_[idxVector[iBlock].first * kSketches_ + iSk] == sketches_[idxVector[iBlock].second * kSketches_ + iSk]);
		}
		// mutex may have a performance hit sometimes
		// ASSUMING everything works correctly, it is not necessary (each thread accesses different vector elements)
		std::lock_guard<std::mutex> lk(mtx_);
		hashJacVec[iBlock] = simVal / static_cast<float>(kSketches_);
	}
}

