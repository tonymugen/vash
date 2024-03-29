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
#include <immintrin.h>

#include "gvarHash.hpp"
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

void BayesicSpace::binarizeBedLocus(const LocationWithLength &bedLocusWindow, const std::vector<char> &bedLocus, const size_t &nIndividuals,
														const LocationWithLength &binLocusWindow, std::vector<uint8_t> &binLocus) {
	RanDraw prng;
	constexpr size_t word64size{8};                                                            // size of uint64_t word in bytes
	constexpr size_t word32size{4};                                                            // size of uint32_t word in bytes
	constexpr size_t word32sizeInBits{32};                                                     // size of uint32_t word in bits
	constexpr size_t word64mask{-word64size};                                                  // for rounding down to nearest divisible by 8 number 
	constexpr size_t word32mask{-word32size};                                                  // for rounding down to nearest divisible by 4 number 
	constexpr uint64_t secondBitMask{0xaaaaaaaaaaaaaaaa};                                      // all bed second bits set
	constexpr uint64_t firstBitMask{~secondBitMask};                                           // all bed first bits set
	const size_t nEvenBedBytes{bedLocusWindow.length & word64mask};                            // number of bed bytes that fully fit into 64-bit words
	std::vector<uint32_t> missWords;                                                           // 32-bit words with missing genotype masks
	std::vector<uint32_t> binWords;                                                            // 32-bit words with binarized data
	uint32_t setCount{0};
	uint32_t missingCount{0};
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
		allHets  &= ranBits;                                                                   // flip some of the het bits randomly
		binBits   = setHoms | allHets;                                                         // add in the set het bits
		setCount += static_cast<uint32_t>( _mm_popcnt_u32(binBits) );                          // count the number of set bits
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
		allHets  &= ranBits;
		binBits   = setHoms | allHets;
		setCount += static_cast<uint32_t>( _mm_popcnt_u32(binBits) );
		binWords.push_back(binBits);
	}
	// test if the set bits are minor alleles and flip them if not
	assert( (nIndividuals > missingCount) && "ERROR: fewer individuals than missing values in binarizeBedLocus()" );      // NOLINT
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

std::vector<IndexedPairSimilarity> BayesicSpace::initializeIPSvector(const LocationWithLength &pairSpan, const size_t &totalNelements) {
	std::vector<IndexedPairSimilarity> pairVector;
	const size_t nnElements{totalNelements * (totalNelements - 1) / 2 - 1};
	pairVector.reserve(pairSpan.length);
	size_t iPair{0};
	while (iPair < pairSpan.length) {
		const size_t curPairInd{iPair + pairSpan.start};
		const size_t kpIdx{nnElements - curPairInd};
		const size_t pIdx = (static_cast<size_t>( sqrt( 1.0 + 8.0 * static_cast<double>(kpIdx) ) ) - 1) / 2;
		const auto row    = static_cast<uint32_t>(totalNelements - 2 - pIdx);
		const auto col    = static_cast<uint32_t>(totalNelements - (kpIdx - pIdx * (pIdx + 1) / 2) - 1);

		IndexedPairSimilarity curPair{};
		curPair.groupID         = 0;
		curPair.similarityValue = 0.0F;
		curPair.element1ind     = row;
		curPair.element2ind     = col;
		pairVector.emplace_back(curPair);

		++iPair;
	}
	return pairVector;
}

std::vector<IndexedPairSimilarity> BayesicSpace::vectorizeGroups(const uint32_t &firstGrpIdx,
				const std::vector< std::vector<uint32_t> >::const_iterator grpBlockStart, const std::vector< std::vector<uint32_t> >::const_iterator grpBlockEnd) {
	std::vector<IndexedPairSimilarity> indexedSimilarityVec;
	uint32_t groupID{firstGrpIdx};
	for (auto eachGrpIt = grpBlockStart; eachGrpIt != grpBlockEnd; ++eachGrpIt) {
		assert( (eachGrpIt->size() > 1) && "ERROR: groups must have at least two elements in vectorizeGroups()" ); // NOLINT
		for (size_t iRow = 0; iRow < eachGrpIt->size() - 1; ++iRow) {
			for (size_t jCol = iRow + 1; jCol < eachGrpIt->size(); ++jCol) {
				indexedSimilarityVec.emplace_back(
							IndexedPairSimilarity{0.0, (*eachGrpIt)[iRow], (*eachGrpIt)[jCol], groupID}
						);
			}
		}
		++groupID;
	}
	return indexedSimilarityVec;
}

void BayesicSpace::saveValues(const std::vector<float> &inVec, std::fstream &outputStream) {
	for (const auto &eachValue : inVec) {
		outputStream << eachValue << " ";
	}
}

void BayesicSpace::saveValues(const std::vector<IndexedPairSimilarity> &inVec, std::fstream &outputStream) {
	for (const auto &eachValue : inVec) {
		outputStream << "G" << eachValue.groupID + 1 << "\t" << eachValue.element1ind + 1 << "\t" << eachValue.element2ind + 1 << "\t" << eachValue.similarityValue << "\n";
	}
}

void BayesicSpace::saveValues(const std::vector<IndexedPairSimilarity> &inVec, const std::vector<std::string> &locusNames, std::fstream &outputStream) {
	for (const auto &eachValue : inVec) {
		outputStream << "G" << eachValue.groupID + 1 << "\t" << locusNames[eachValue.element1ind] << "\t" << locusNames[eachValue.element2ind] << "\t" << eachValue.similarityValue << "\n";
	}
}

void BayesicSpace::saveValues(const std::vector<IndexedPairLD> &inVec, std::fstream &outputStream) {
	for (const auto &eachValue : inVec) {
		outputStream << eachValue.element1ind + 1 << "\t" << eachValue.element2ind + 1 << "\t" << eachValue.jaccard << "\t" << eachValue.rSq << "\n";
	}
}

void BayesicSpace::saveValues(const std::vector<IndexedPairLD> &inVec, const std::vector<std::string> &locusNames, std::fstream &outputStream) {
	for (const auto &eachValue : inVec) {
		outputStream << locusNames[eachValue.element1ind] << "\t" << locusNames[eachValue.element2ind] << "\t" << eachValue.jaccard << "\t" << eachValue.rSq << "\n";
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

void BayesicSpace::extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI,
			std::unordered_map<std::string, int> &intVariables, std::unordered_map<std::string, std::string> &stringVariables) {
	intVariables.clear();
	stringVariables.clear();
	const std::array<std::string, 1> requiredStringVariables{"input-bed"};
	const std::array<std::string, 4> optionalStringVariables{"log-file", "out-file", "only-groups", "add-locus-names"};
	const std::array<std::string, 1> requiredIntVariables{"n-individuals"};
	const std::array<std::string, 3> optionalIntVariables{"hash-size", "threads", "n-rows-per-band"};

	const std::unordered_map<std::string, int>         defaultIntValues{ {"hash-size", 0}, {"threads", -1}, {"n-rows-per-band", 0} };
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
