/*
 * Copyright (c) 2023 Anthony J. Greenberg
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

/// Auxiliary functions for variant hashing
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023 Anthony J. Greenberg
 * \version 0.5
 *
 * Implementation of class-external functions needed by hashing classes.
 *
 */

#include <cstring>
#include <cassert>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <utility>  // for std::pair
#include <immintrin.h>

#include "vashFunctions.hpp"

using namespace BayesicSpace;

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
	uint64_t totSet{0};
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
	uint64_t totSet{0};
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

uint32_t BayesicSpace::murMurHashMixer(const size_t &key, const uint32_t &seed){
	constexpr uint32_t const1{0xcc9e2d51};
	constexpr uint32_t const2{0x1b873593};
	constexpr size_t   nBlocks32{sizeof(size_t) / sizeof(uint32_t)};    // number of 32 bit blocks in size_t
	constexpr uint32_t keyLen{sizeof(size_t)};                          // key length 
	constexpr uint32_t hashMultiplier{5};
	constexpr uint32_t hashAdder{0xe6546b64};

	constexpr std::array<uint32_t, 4> blockShifts{15, 17, 13, 19};

	uint32_t hash{seed};
	std::array<uint32_t, nBlocks32> blocks{0};
	memcpy(blocks.data(), &key, keyLen);

	// body
	for (auto &eachBlock : blocks){
		eachBlock *= const1;
		eachBlock  = (eachBlock << blockShifts[0]) | (eachBlock >> blockShifts[1]);
		eachBlock *= const2;

		hash ^= eachBlock;
		hash  = (hash << blockShifts[2]) | (hash >> blockShifts[3]);
		hash  = hash * hashMultiplier + hashAdder;
	}
	return hash;
}

uint32_t BayesicSpace::murMurHashFinalizer(const uint32_t &inputHash){
	constexpr uint32_t keyLen{sizeof(size_t)};                          // key length 
	constexpr std::array<uint32_t, 2> finalizeShifts{16, 13};
	constexpr std::array<uint32_t, 2> finalizeMult{0x85ebca6b, 0xc2b2ae35};

	uint32_t hash = inputHash;
	hash ^= keyLen;
	hash ^= hash >> finalizeShifts[0];
	hash *= finalizeMult[0];
	hash ^= hash >> finalizeShifts[1];
	hash *= finalizeMult[1];
	hash ^= hash >> finalizeShifts[0];

	return hash;
}

uint32_t BayesicSpace::murMurHash(const size_t &key, const uint32_t &seed){
	uint32_t hash{murMurHashMixer(key, seed)};
	hash = murMurHashFinalizer(hash);

	return hash;
}

uint32_t BayesicSpace::murMurHash(const std::vector<size_t> &key, const uint32_t &seed) {
	uint32_t hash{seed};
	for (const auto &eachIdx : key){
		hash = murMurHashMixer(eachIdx, hash);
	}
	hash = murMurHashFinalizer(hash);
	return hash;
}

uint32_t BayesicSpace::murMurHash(const size_t &start, const size_t &length, const std::vector<uint16_t> &key, const uint32_t &seed) {
	constexpr size_t keysPerWord{sizeof(size_t) / sizeof(uint16_t)};
	constexpr auto roundMask = static_cast<size_t>( -(keysPerWord) );
	uint32_t hash{seed};
	const size_t end{start + length};
	const size_t wholeEnd{end & roundMask};

	assert( ( end < key.size() ) && "ERROR: length goes past the end of key in murMurHash" );

	size_t keyIdx{start};
	while (keyIdx < wholeEnd){
		size_t keyBlock{0};
		memcpy(&keyBlock, key.data() + keyIdx, keysPerWord);
		hash    = murMurHashMixer(keyBlock, hash);
		keyIdx += keysPerWord;
	}
	if (end > wholeEnd){  // if there is a tail
		size_t keyBlock{0};
		memcpy(&keyBlock, key.data() + keyIdx, keysPerWord);
		hash = murMurHashMixer(keyBlock, hash);
	}
	hash = murMurHashFinalizer(hash);
	return hash;
}

void BayesicSpace::testBedMagicBytes(std::array<char, 3> &bytesToTest) {
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

std::vector< std::pair<size_t, size_t> > BayesicSpace::makeThreadRanges(const size_t &nThreads, const size_t &nElementsPerThread){
	std::vector< std::pair<size_t, size_t> > threadRanges;
	size_t bedInd = 0;
	for (size_t iThread = 0; iThread < nThreads; ++iThread){
		threadRanges.emplace_back(bedInd, bedInd + nElementsPerThread);
		bedInd += nElementsPerThread;
	}
	return threadRanges;
}

void BayesicSpace::saveValues(const std::vector<float> &inVec, std::fstream &outputStream) {
	for (const auto &eachValue : inVec){
		outputStream << eachValue << " ";
	}
}


