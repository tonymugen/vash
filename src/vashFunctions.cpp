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
#include <unordered_map>
#include <utility>  // for std::pair
#include <limits>
#include <algorithm>
#include <immintrin.h>

#include "random.hpp"
#include "gvarHash.hpp"
#include "similarityMatrix.hpp"
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
	while (iByte < nWholeWords) {
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, wordSize);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
		iByte += wordSize;
	}
	if ( nWholeWords < inVec.size() ) {
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, inVec.size() - nWholeWords);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
	}
	return totSet;
}

uint64_t BayesicSpace::countSetBits(const std::vector<uint8_t> &inVec, const LocationWithLength &window) {
	constexpr size_t wordSize{8};
	constexpr uint64_t roundMask{0xfffffffffffffff8};
	uint64_t totSet{0};
	const size_t roundLength{window.length & roundMask};
	const size_t nWholeWords{window.start + roundLength};
	size_t iByte{window.start};
	while (iByte < nWholeWords) {
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, wordSize);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
		iByte += wordSize;
	}
	if (roundLength < window.length) {
		uint64_t chunk{0};
		memcpy(&chunk, inVec.data() + iByte, window.length - roundLength);
		totSet += static_cast<uint64_t>( _mm_popcnt_u64(chunk) );
	}
	return totSet;
}

size_t BayesicSpace::getAvailableRAM() {
	constexpr size_t defaultSize{2147483648};
	if ( std::ifstream("/proc/meminfo").good() ) {
		constexpr size_t bytesInKb{1024};
		constexpr size_t memTokenLength{13};
		std::string memLine;
		std::fstream memInfoStream;
		memInfoStream.open("/proc/meminfo", std::ios::in);
		while ( getline(memInfoStream, memLine) ) {
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

uint32_t BayesicSpace::murMurHashMixer(const std::array<uint32_t, SIZE_OF_SIZET> &key, const uint32_t &seed) {
	constexpr uint32_t const1{0xcc9e2d51};
	constexpr uint32_t const2{0x1b873593};
	constexpr uint32_t hashMultiplier{5};
	constexpr uint32_t hashAdder{0xe6546b64};
	constexpr std::array<uint32_t, 4> blockShifts{15, 17, 13, 19};

	uint32_t hash{seed};

	// body
	for (auto eachBlock : key) {
		eachBlock *= const1;
		eachBlock  = (eachBlock << blockShifts[0]) | (eachBlock >> blockShifts[1]);
		eachBlock *= const2;

		hash ^= eachBlock;
		hash  = (hash << blockShifts[2]) | (hash >> blockShifts[3]);
		hash  = hash * hashMultiplier + hashAdder;
	}
	return hash;
}

uint32_t BayesicSpace::murMurHashFinalizer(const uint32_t &inputHash) {
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

uint32_t BayesicSpace::murMurHash(const std::array<uint32_t, SIZE_OF_SIZET> &key, const uint32_t &seed) {
	uint32_t hash{murMurHashMixer(key, seed)};
	hash = murMurHashFinalizer(hash);

	return hash;
}

uint32_t BayesicSpace::murMurHash(const std::vector<size_t> &key, const uint32_t &seed) {
	uint32_t hash{seed};
	for (const auto eachIdx : key) {
		std::array<uint32_t, SIZE_OF_SIZET> eachKey{};
		memcpy(eachKey.data(), &eachIdx, SIZE_OF_SIZET);
		hash = murMurHashMixer(eachKey, hash);
	}
	hash = murMurHashFinalizer(hash);
	return hash;
}

uint32_t BayesicSpace::murMurHash(const std::vector<uint32_t> &key, const uint32_t &seed) {
	uint32_t hash{seed};
	const size_t lastEven = key.size() & static_cast<uint32_t>(-2); // round down the size to even number
	for (size_t iKey = 0; iKey < lastEven; ++iKey) {
		std::array<uint32_t, SIZE_OF_SIZET> eachKey{};
		memcpy(eachKey.data(), key.data() + iKey, SIZE_OF_SIZET);
		hash = murMurHashMixer(eachKey, hash);
	}
	std::array<uint32_t, SIZE_OF_SIZET> lastKey{};
	memcpy(lastKey.data(), &key.back(), key.size() - lastEven);
	hash = murMurHashMixer(lastKey, hash);

	hash = murMurHashFinalizer(hash);
	return hash;
}

uint32_t BayesicSpace::murMurHash(const std::vector<uint16_t> &key, const LocationWithLength &keyWindow, const uint32_t &seed) {
	constexpr size_t roundMask{-SIZE_OF_SIZET};
	uint32_t hash{seed};
	const size_t end{keyWindow.start + keyWindow.length};
	const size_t wholeEnd{end & roundMask};

	assert( ( end <= key.size() ) && "ERROR: length goes past the end of key in murMurHash" ); // NOLINT

	size_t keyIdx{keyWindow.start};
	while (keyIdx < wholeEnd) {
		std::array<uint32_t, SIZE_OF_SIZET> keyBlock{};
		memcpy(keyBlock.data(), key.data() + keyIdx, SIZE_OF_SIZET);
		hash    = murMurHashMixer(keyBlock, hash);
		keyIdx += SIZE_OF_SIZET;
	}
	if (end > wholeEnd) {  // if there is a tail
		std::array<uint32_t, SIZE_OF_SIZET> keyBlock{};
		memcpy(keyBlock.data(), key.data() + keyIdx, SIZE_OF_SIZET);
		hash = murMurHashMixer(keyBlock, hash);
	}
	hash = murMurHashFinalizer(hash);
	return hash;
}

void BayesicSpace::testBedMagicBytes(const std::array<char, N_BED_TEST_BYTES> &bytesToTest) {
	constexpr std::array<char, N_BED_TEST_BYTES> magicBytes{0x6c, 0x1b, 0x01};   // Leading bytes for .bed files
	if (bytesToTest[0] != magicBytes[0]) {
		throw std::string("ERROR: first magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (bytesToTest[1] != magicBytes[1]) {
		throw std::string("ERROR: second magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (bytesToTest[2] != magicBytes[2]) {
		throw std::string("ERROR: third magic byte in input .bed file is not the expected value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
}

std::vector< std::pair<size_t, size_t> > BayesicSpace::makeThreadRanges(const CountAndSize &threadPoolSizes) {
	std::vector< std::pair<size_t, size_t> > threadRanges;
	size_t bedInd = 0;
	for (size_t iThread = 0; iThread < threadPoolSizes.count; ++iThread) {
		threadRanges.emplace_back(bedInd, bedInd + threadPoolSizes.size);
		bedInd += threadPoolSizes.size;
	}
	return threadRanges;
}

std::vector<size_t> BayesicSpace::makeChunkSizes(const size_t &nElements, const size_t &nChunks) {
	std::vector<size_t> chunkSizes(nChunks, nElements / nChunks);
	// spread the left over elements among chunks
	std::for_each(
		chunkSizes.begin(),
		chunkSizes.begin() + static_cast<std::vector<size_t>::difference_type >(nElements % nChunks),
		[](size_t &currSize) {return ++currSize;}
	);
	return chunkSizes;
}

std::vector< std::pair<RowColIdx, RowColIdx> > BayesicSpace::makeChunkRanges(const LocationWithLength &startAndChunkSize, const size_t nChunks) {
	std::vector<size_t> chunkSizes{makeChunkSizes(startAndChunkSize.length, nChunks)};

	std::vector< std::pair<RowColIdx, RowColIdx> > chunkPairs;
	size_t cumChunkSize{startAndChunkSize.start};
	std::for_each(
		chunkSizes.cbegin(),
		chunkSizes.cend(),
		[&chunkPairs, &cumChunkSize](const size_t &chunkSize) {
			std::pair<RowColIdx, RowColIdx> tmpPair;
			tmpPair.first  = recoverRCindexes(cumChunkSize);
			tmpPair.second = recoverRCindexes(cumChunkSize + chunkSize);
			chunkPairs.emplace_back(tmpPair);
			cumChunkSize += chunkSize;
		}
	);

	return chunkPairs;
}

std::pair<HashGroupItPairCount, HashGroupItPairCount>
	BayesicSpace::makeGroupRanges(const std::vector<HashGroup> &groupVector, const HashGroupItPairCount &startHGPC, const size_t &chunkSize) {
	std::pair<HashGroupItPairCount, HashGroupItPairCount> result;
	result.first = startHGPC;
	const size_t startPairCount{
		startHGPC.hgIterator->cumulativeNpairs - startHGPC.hgIterator->locusIndexes.size() * (startHGPC.hgIterator->locusIndexes.size() - 1) / 2
			+ startHGPC.pairCount
	};
	const size_t chunkCutOff = std::min(startPairCount + chunkSize, groupVector.back().cumulativeNpairs);
	const auto gvIterator = std::lower_bound(
		startHGPC.hgIterator,
		groupVector.cend(),
		chunkCutOff,
		[](const HashGroup &currGrp, const uint64_t &cutoff) {return currGrp.cumulativeNpairs < cutoff;}
	);
	result.second.hgIterator = gvIterator;
	if ( gvIterator == groupVector.end() ) {
		result.second.pairCount = 0;
		return result;
	}
	const size_t lastGroupPairNumber{
		gvIterator->locusIndexes.size() * (gvIterator->locusIndexes.size() - 1) / 2 - (gvIterator->cumulativeNpairs - chunkCutOff)
	};
	result.second.pairCount = lastGroupPairNumber;

	return result;
}

void BayesicSpace::binarizeBedLocus(const LocationWithLength &bedLocusWindow, const std::vector<char> &bedLocus, const size_t &nIndividuals,
														const LocationWithLength &binLocusWindow, std::vector<uint8_t> &binLocus) {
	RanDraw prng;
	constexpr size_t word64size{8};                                                            // size of uint64_t word in bytes
	constexpr size_t word32size{4};                                                            // size of uint32_t word in bytes
	constexpr size_t word32sizeInBits{32};                                                     // size of uint32_t word in bits
	constexpr size_t word64mask{-word64size};                                                  // for rounding down to nearest divisible by 8 number 
	constexpr size_t word32mask{-word32size};                                                  // for rounding down to nearest divisible by 4 number 
	constexpr uint64_t secondBitMask{0xAAAAAAAAAAAAAAAA};                                      // all bed second bits set
	constexpr uint64_t firstBitMask{~secondBitMask};                                           // all bed first bits set
	const size_t nEvenBedBytes{bedLocusWindow.length & word64mask};                            // number of bed bytes that fully fit into 64-bit words
	std::vector<uint32_t> missWords;                                                           // 32-bit words with missing genotype masks
	std::vector<uint32_t> binWords;                                                            // 32-bit words with binarized data
	uint32_t setCount{0};
	uint32_t missingCount{0};
	uint32_t extraHets{0};     // unset heterozygotes; needed to decide whether to bit-flip the result
	size_t iBedByte{bedLocusWindow.start};
	while (iBedByte < bedLocusWindow.start + nEvenBedBytes) {
		uint64_t bedWord{0};
		memcpy(&bedWord, bedLocus.data() + iBedByte, word64size);
		auto binBits{static_cast<uint32_t>( _pext_u64(bedWord, firstBitMask) )};
		const auto secondBedBits{static_cast<uint32_t>( _pext_u64(bedWord, secondBitMask) )};
		const uint32_t setHoms{binBits & secondBedBits};                                       // set bit homozygotes (11 in .bed)
		uint32_t allHets{(~setHoms) & secondBedBits};                                          // set at all het positions
		const uint32_t missing{(~setHoms) & binBits};                                          // set at missing positions
		missingCount += static_cast<uint32_t>( _mm_popcnt_u32(missing) );                      // count missing genotypes
		missWords.push_back(missing);
		const auto ranBits{static_cast<uint32_t>( prng.ranInt() )};
		const auto allHetsCount{static_cast<uint32_t>( _mm_popcnt_u32(allHets) )};
		allHets   &= ranBits;                                                                  // flip some of the het bits randomly
		extraHets += allHetsCount - static_cast<uint32_t>( _mm_popcnt_u32(allHets) );
		binBits    = setHoms | allHets;                                                        // add in the set het bits
		setCount  += static_cast<uint32_t>( _mm_popcnt_u32(binBits) );                         // count the number of set bits
		binWords.push_back(binBits);
		iBedByte += word64size;
	}
	if (bedLocusWindow.length > nEvenBedBytes) {
		const size_t nTrailingBedBytes{bedLocusWindow.length - nEvenBedBytes};
		uint64_t bedWord{0};
		memcpy(&bedWord, bedLocus.data() + iBedByte, nTrailingBedBytes);
		auto binBits{static_cast<uint32_t>( _pext_u64(bedWord, firstBitMask) )};
		const auto secondBedBits{static_cast<uint32_t>( _pext_u64(bedWord, secondBitMask) )};
		const uint32_t setHoms{binBits & secondBedBits};
		uint32_t allHets{(~setHoms) & secondBedBits};
		const uint32_t missing{(~setHoms) & binBits};
		missingCount += static_cast<uint32_t>( _mm_popcnt_u32(missing) );
		missWords.push_back(missing);
		const auto ranBits{static_cast<uint32_t>( prng.ranInt() )};
		const auto allHetsCount{static_cast<uint32_t>( _mm_popcnt_u32(allHets) )};
		allHets   &= ranBits;                                                                  // flip some of the het bits randomly
		extraHets += allHetsCount - static_cast<uint32_t>( _mm_popcnt_u32(allHets) );
		binBits    = setHoms | allHets;                                                        // add in the set het bits
		setCount  += static_cast<uint32_t>( _mm_popcnt_u32(binBits) );
		binWords.push_back(binBits);
	}
	// test if the set bits are minor alleles and flip them if not
	assert( (nIndividuals > missingCount) && "ERROR: fewer individuals than missing values in binarizeBedLocus()" );      // NOLINT
	setCount += extraHets;
	if ( (2UL * setCount) > (nIndividuals - missingCount) ) {
		size_t missWordIdx{0};
		for (auto &eachBinWord : binWords) {
			eachBinWord = (~eachBinWord) & (~missWords[missWordIdx]);                                                     // flip the bits in the binary vector and reset the missing bits to 0
			++missWordIdx;
		}
		const auto padShift{static_cast<uint32_t>( (binWords.size() * word32sizeInBits) % nIndividuals )};
		const uint32_t lastWordMask{std::numeric_limits<uint32_t>::max() >> padShift};                                    // clear the padding bits after the flip
		binWords.back() = binWords.back() & lastWordMask;
	}
	// copy over the binary bits to the locus vector
	const size_t nEvenBinBytes{binLocusWindow.length & word32mask};                                                       // number of bin bytes that fully fit into 32-bit words
	size_t iBinByte{binLocusWindow.start};
	size_t iBinWord{0};
	while (iBinByte < binLocusWindow.start + nEvenBinBytes) {
		memcpy(binLocus.data() + iBinByte, binWords.data() + iBinWord, word32size);
		iBinByte += word32size;
		++iBinWord;
	}
	if (binLocusWindow.length > nEvenBinBytes) {
		const size_t nTrailingBinBytes{binLocusWindow.length - nEvenBinBytes};
		memcpy(binLocus.data() + iBinByte, binWords.data() + iBinWord, nTrailingBinBytes);
	}
}

void BayesicSpace::binarizeMacLocus(const std::vector<int> &macLocus, const LocationWithLength &binLocusWindow, std::vector<uint8_t> &binLocus) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	constexpr uint8_t byteSize{8};
	constexpr uint8_t oneBit{0b00000001};       // One set bit for masking
	constexpr uint8_t middleMask{0b10000011};
	constexpr uint8_t endTwoBitMask{0b00000011};
	RanDraw locPRNG;
	auto remainderInd          = static_cast<uint8_t>( binLocusWindow.length * byteSize - macLocus.size() );
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
	remainderInd               = byteSize - remainderInd;
	// Create a vector to store random bytes for stochastic heterozygote resolution
	const size_t randVecLen{macLocus.size() / sizeof(uint64_t) + static_cast<size_t>( ( macLocus.size() % sizeof(uint64_t) ) > 0 )};
	std::vector<uint64_t> rand(randVecLen);
	auto *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	// Fill the random byte vector
	for (auto &randValue : rand) {
		randValue = locPRNG.ranInt();
	}
	size_t nMissing{0};
	size_t extraHets{0};                                                                                       // unset heterozygotes; needed to decide whether to bit-flip the result
	size_t iIndiv{0};
	size_t iRB{0};                                                                                             // randBytes index
	const size_t begByte{binLocusWindow.start * binLocusWindow.length};
	std::vector<uint8_t> missMasks(binLocusWindow.length, 0);
	size_t i0Byte{0};                                                                                          // to index the missMasks vector
	for (size_t iByte = begByte; iByte < begByte + binLocusWindow.length - 1; ++iByte) {                       // treat the last byte separately
		uint8_t binByte{0};
		for (uint8_t iInByte = 0; iInByte < byteSize; ++iInByte) {
			auto curIndiv{static_cast<uint8_t>(macLocus[iIndiv])};                                             // cramming down to one byte because I do not care what the actual value is
			curIndiv &= middleMask;                                                                            // mask everything in the middle
			const auto missingMask{static_cast<uint8_t>(curIndiv >> 7)};                                       // 0b00000001 iff is missing (negative value)
			nMissing          += missingMask;
			missMasks[i0Byte] |= static_cast<uint8_t>(missingMask << iInByte);
			curIndiv          &= endTwoBitMask;
			const auto randMask{static_cast<uint8_t>( (randBytes[iRB] >> iInByte) & oneBit )};                 // 0b00000000 or 0b00000001 with equal chance
			auto curBitMask{static_cast<uint8_t>( (curIndiv >> 1) ^ (curIndiv & randMask) )};                  // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			const auto allHet{static_cast<uint32_t>( (curIndiv & oneBit) & (~missingMask) )};                  // 1 only if het, regardless of the random bit
			extraHets  += allHet - (randMask & allHet);                                                        // subtract 1 iff the randMask is set and allHet is not 0
			curBitMask &= ~missingMask;                                                                        // zero it out if missing value is set
			binByte    |= static_cast<uint8_t>(curBitMask << iInByte);
			++iIndiv;
		}
		// should be safe: each thread accesses different vector elements
		binLocus[iByte] = binByte;
		++i0Byte;
		++iRB;
	}
	// now deal with the last byte in the individual
	uint8_t lastBinByte{0};
	for (uint8_t iRem = 0; iRem < remainderInd; ++iRem) {
		auto curIndiv{static_cast<uint8_t>(macLocus[iIndiv])};                                                 // cramming down to one byte because I do not care what the actual value is
		curIndiv &= middleMask;                                                                                // mask everything in the middle
		const auto missingMask{static_cast<uint8_t>(curIndiv >> 7)};                                           // 0b00000001 iff is missing (negative value)
		nMissing         += missingMask;
		missMasks.back() |= static_cast<uint8_t>(missingMask << iRem);
		curIndiv         &= endTwoBitMask;                                                
		const auto randMask{static_cast<uint8_t>( (randBytes[binLocusWindow.length - 1] >> iRem) & oneBit )};  // 0b00000000 or 0b00000001 with equal chance
		auto curBitMask{static_cast<uint8_t>( (curIndiv >> 1) ^ (curIndiv & randMask) )};                      // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
		const auto allHet{static_cast<uint32_t>( (curIndiv & oneBit) & (~missingMask) )};                      // 1 only if het, regardless of the random bit
		extraHets   += allHet - (randMask & allHet);                                                           // subtract 1 iff the randMask is set and allHet is not 0
		curBitMask  &= ~missingMask;                                                                           // zero it out if missing value is set
		lastBinByte |= static_cast<uint8_t>(curBitMask << iRem);
		++iIndiv;
	}
	// should be safe: each thread accesses different vector elements
	binLocus[begByte + binLocusWindow.length - 1] = lastBinByte;
	const LocationWithLength genoWindow{begByte, binLocusWindow.length};
	const auto setCount{countSetBits(binLocus, genoWindow) + extraHets};
	if (2UL * setCount > macLocus.size() - nMissing) { // always want the alternative to be the minor allele
		i0Byte = 0;
		for (size_t i = begByte; i < begByte + binLocusWindow.length; ++i) {
			// should be safe: each thread accesses different vector elements
			binLocus[i] = (~binLocus[i]) & (~missMasks[i0Byte]);
			++i0Byte;
		}
		// should be safe: each thread accesses different vector elements
		binLocus[begByte + binLocusWindow.length - 1] &= lastByteMask; // unset the remainder bits
	}
}

std::vector<std::string> BayesicSpace::getLocusNames(const std::string &bimFileName) {
	std::fstream bimStream;
	std::vector<std::string> locusNames;
	bimStream.open(bimFileName, std::ios::in);
	std::string bimLine;
	std::stringstream lineStream;
	while ( std::getline(bimStream, bimLine) ) {
		lineStream.str(bimLine);
		std::string field;
		lineStream >> field >> field;             // locus name is the second column in the .bim file
		locusNames.emplace_back(field);
	}
	bimStream.close();
	return locusNames;
}

void BayesicSpace::parseCL(int &argc, char **argv, std::unordered_map<std::string, std::string> &cli) {
	// set to true after encountering a flag token (the characters after the dash)
	bool val = false;
	// store the token value here
	std::string curFlag;

	for (int iArg = 1; iArg < argc; iArg++) {
		const char *pchar = argv[iArg];
		if ( (pchar[0] == '-') && (pchar[1] == '-') ) { // encountered the double dash, look for the token after it
			if (val) { // A previous flag had no value
				cli[curFlag] = "set";
			}
			// what follows the dash?
			val     = true;
			curFlag = pchar + 2;
		} else {
			if (val) {
				val          = false;
				cli[curFlag] = pchar;
			}
		}
	}
}

void BayesicSpace::extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI, std::unordered_map<std::string, int> &intVariables,
		std::unordered_map<std::string, float> &floatVariables, std::unordered_map<std::string, std::string> &stringVariables) {
	intVariables.clear();
	stringVariables.clear();
	const std::array<std::string, 1> requiredStringVariables{"input-bed"};
	const std::array<std::string, 4> optionalStringVariables{"log-file", "out-file", "only-groups", "add-locus-names"};
	const std::array<std::string, 1> requiredIntVariables{"n-individuals"};
	const std::array<std::string, 3> optionalIntVariables{"hash-size", "threads", "n-rows-per-band"};
	const std::array<std::string, 1> optionalFloatVariables{"min-similarity"};

	const std::unordered_map<std::string, int>         defaultIntValues{ {"hash-size", 0}, {"threads", -1}, {"n-rows-per-band", 0} };
	const std::unordered_map<std::string, float>       defaultFloatValues{ {"min-similarity", 0.0F} };
	const std::unordered_map<std::string, std::string> defaultStringValues{ {"log-file", "ldblocks.log"}, {"out-file", "ldblocksOut.tsv"},
																			{"only-groups", "unset"}, {"add-locus-names", "unset"} };

	if ( parsedCLI.empty() ) {
		throw std::string("No command line flags specified;");
	}
	for (const auto &eachFlag : requiredIntVariables) {
		try {
			intVariables[eachFlag] = stoi( parsedCLI.at(eachFlag));
		} catch(const std::exception &problem) {
			throw std::string("ERROR: ") + eachFlag + std::string(" specification is required and must be an integer");
		}
	}
	for (const auto &eachFlag : optionalIntVariables) {
		try {
			intVariables[eachFlag] = stoi( parsedCLI.at(eachFlag));
		} catch(const std::exception &problem) {
			intVariables[eachFlag] = defaultIntValues.at(eachFlag);
		}
	}
	for (const auto &eachFlag : optionalFloatVariables) {
		try {
			floatVariables[eachFlag] = stof( parsedCLI.at(eachFlag));
		} catch(const std::exception &problem) {
			floatVariables[eachFlag] = defaultFloatValues.at(eachFlag);
		}
	}
	for (const auto &eachFlag : requiredStringVariables) {
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			throw std::string("ERROR: ") + eachFlag + std::string(" specification is required");
		}
	}
	for (const auto &eachFlag : optionalStringVariables) {
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			stringVariables[eachFlag] = defaultStringValues.at(eachFlag);
		}
	}
}
