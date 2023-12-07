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

/// Unit tests
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023 Anthony J. Greenberg
 * \version 0.1
 *
 * Unit tests using Catch2.
 *
 */

#include <cstdint>
#include <cstdio>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>

#include <iostream>

#include "random.hpp"
#include "gvarHash.hpp"
#include "vashFunctions.hpp"
#include "similarityMatrix.hpp"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

// Number of times tests of random events will be run
static constexpr uint16_t N_RAN_ITERATIONS{10};
// precision for float comparisons
static constexpr float FPREC{1e-4F};

TEST_CASE("Can count set bits correctly", "[countSetBits]") {
	constexpr uint16_t oneWord{0b11001110'01101001};
	constexpr uint16_t wCorrectCount{9};
	const uint16_t wordBitCount{BayesicSpace::countSetBits(oneWord)};
	REQUIRE(wordBitCount == wCorrectCount);
	constexpr size_t byteArraySize{14};
	constexpr std::array<uint8_t, byteArraySize> byteArray{
		0b11110111, 0b10100111, 0b01011000, 0b11111001, 0b11100110, 0b11111100, 0b00111111,
		0b01101011, 0b01011111, 0b01001001, 0b11001001, 0b11100010, 0b11010101, 0b01110001
	};
	constexpr uint64_t correctFullVectorCount{69};
	std::vector<uint8_t> byteVector( byteArray.begin(), byteArray.end() );
	const uint64_t fullVectorCount{BayesicSpace::countSetBits(byteVector)};
	REQUIRE(fullVectorCount == correctFullVectorCount);
	constexpr BayesicSpace::LocationWithLength byteVectorWindow{2, 4};
	constexpr uint64_t correctVectorWindowCount{20};
	const uint64_t vectorWindowCount{BayesicSpace::countSetBits(byteVector, byteVectorWindow)};
	REQUIRE(vectorWindowCount == correctVectorWindowCount);
}

TEST_CASE("MurMurHash works properly", "[MurMurHash]") {
	constexpr uint32_t mmHashSeed{2153025618};
	constexpr std::array<uint32_t, BayesicSpace::SIZE_OF_SIZET> arrayKey{335636695, 4242517348};
	constexpr uint32_t correctArrayHash{2730141477};
	const uint32_t arrayMMhash{BayesicSpace::murMurHash(arrayKey, mmHashSeed)};
	constexpr std::array<size_t, 11> idxArray{
		2437, 2444, 41116, 42353,
		45949, 58374, 75248, 80113,
		93649, 98640, 99638
	};
	std::vector<size_t> idxVector( idxArray.begin(), idxArray.end() );
	constexpr uint32_t correctIdxVectorHash{2643649892};
	const uint32_t idxVectorHash{BayesicSpace::murMurHash(idxVector, mmHashSeed)};
	constexpr std::array<uint32_t, 11> u32Array{
		2437, 2444, 41116, 42353,
		45949, 58374, 75248, 80113,
		93649, 98640, 99638
	};
	std::vector<size_t> u32Vector( u32Array.begin(), u32Array.end() );
	constexpr uint32_t correctU32VectorHash{2643649892};
	const uint32_t u32VectorHash{BayesicSpace::murMurHash(idxVector, mmHashSeed)};
	constexpr std::array<uint16_t, 13> array16bit{
		1256, 2117, 2866, 7434,
		11737, 16256, 22236, 39883,
		40023, 46299, 58123, 58167, 62187
	};
	std::vector<uint16_t> vector16bit( array16bit.begin(), array16bit.end() );
	constexpr BayesicSpace::LocationWithLength keyWindow{4, 5};
	constexpr uint32_t correct16bitHash{3760365877};
	const uint32_t v16bitHash{BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed)};
	constexpr BayesicSpace::LocationWithLength wholeKeySpan{0, array16bit.size()};
	constexpr uint32_t correctAll16bitHash{2280422248};
	const uint32_t all16bitHash{BayesicSpace::murMurHash(vector16bit, wholeKeySpan, mmHashSeed)};
	SECTION("MurMurHash correctness tests") {
		REQUIRE(arrayMMhash   == correctArrayHash);
		REQUIRE(idxVectorHash == correctIdxVectorHash);
		REQUIRE(u32VectorHash == correctU32VectorHash);
		REQUIRE(u32VectorHash == idxVectorHash);
		REQUIRE(v16bitHash    == correct16bitHash);
		REQUIRE(all16bitHash  == correctAll16bitHash);
	}
	SECTION("MurMurHash sensitivity tests") {
		constexpr uint32_t mmHashSeed2{2153025619};
		constexpr std::array<uint32_t, BayesicSpace::SIZE_OF_SIZET> arrayKey2{335636696, 4242517348};
		REQUIRE(BayesicSpace::murMurHash(arrayKey, mmHashSeed2)  != correctArrayHash);
		REQUIRE(BayesicSpace::murMurHash(arrayKey2, mmHashSeed)  != correctArrayHash);
		REQUIRE(BayesicSpace::murMurHash(idxVector, mmHashSeed2) != correctIdxVectorHash);
		idxVector.at(1)++;
		REQUIRE(BayesicSpace::murMurHash(idxVector, mmHashSeed)  != correctIdxVectorHash);
		REQUIRE(BayesicSpace::murMurHash(u32Vector, mmHashSeed2) != correctU32VectorHash);
		u32Vector.at(1)++;
		REQUIRE(BayesicSpace::murMurHash(u32Vector, mmHashSeed)  != correctU32VectorHash);
		vector16bit.at(2)++;
		REQUIRE(BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed)  == correct16bitHash);
		REQUIRE(BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed2) != correct16bitHash);
		vector16bit.at(keyWindow.start + 2)++;
		REQUIRE(BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed)  != correct16bitHash);
	}
}

TEST_CASE(".bed related file and data parsing works", "[bedData]") {
	SECTION("Thread ranges") {
		constexpr BayesicSpace::CountAndSize threadSizes{4, 13};
		const std::vector< std::pair<size_t, size_t> > correctRanges{ {0, 13}, {13, 26}, {26, 39}, {39, 52} };
		const std::vector< std::pair<size_t, size_t> > threadRanges{BayesicSpace::makeThreadRanges(threadSizes)};
		REQUIRE(threadRanges.size() == threadSizes.count);
		REQUIRE(std::all_of(
				threadRanges.cbegin(),
				threadRanges.cend(),
				[](const std::pair<size_t, size_t> &eachRange){return eachRange.first <= eachRange.second;}
			)
		);
		REQUIRE(std::equal(
				threadRanges.cbegin(),
				threadRanges.cend(),
				correctRanges.cbegin(),
				[](const std::pair<size_t, size_t> &pair1, const std::pair<size_t, size_t>&pair2){
					return (pair1.first == pair2.first) && (pair1.second == pair2.second);
				}
			)
		);
	}

	SECTION("Pair index initialization") {
		constexpr BayesicSpace::LocationWithLength pairSpan{71, 201};
		constexpr size_t allN{397};
		std::vector<BayesicSpace::IndexedPairSimilarity> pairSegment{BayesicSpace::initializeIPSvector(pairSpan, allN)};
		REQUIRE(pairSegment.size() == pairSpan.length);
		REQUIRE(pairSegment[0].element1ind == 0);
		REQUIRE(pairSegment[0].element2ind == pairSpan.start + 1);
		REQUIRE(std::all_of(
				pairSegment.cbegin(),
				pairSegment.cend(),
				[](const BayesicSpace::IndexedPairSimilarity &obj){return obj.element1ind < obj.element2ind;}
			)
		);
		REQUIRE(std::all_of(
				pairSegment.cbegin(),
				pairSegment.cend(),
				[&allN](const BayesicSpace::IndexedPairSimilarity &obj){return (obj.element1ind < allN) || (obj.element2ind < allN);}
			)
		);
	}

	SECTION("Magic byte testing") {
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> correctBedBytes{0x6c, 0x1b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes1{ 0x6d, 0x1b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes2{ 0x6c, 0x0b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes3{ 0x6c, 0x1b, 0x11};
		REQUIRE_NOTHROW( BayesicSpace::testBedMagicBytes(correctBedBytes) );
		REQUIRE_THROWS_WITH(BayesicSpace::testBedMagicBytes(wrongBedBytes1), Catch::Matchers::StartsWith("ERROR: first magic byte in input .bed file") );
		REQUIRE_THROWS_WITH(BayesicSpace::testBedMagicBytes(wrongBedBytes2), Catch::Matchers::StartsWith("ERROR: second magic byte in input .bed file") );
		REQUIRE_THROWS_WITH(BayesicSpace::testBedMagicBytes(wrongBedBytes3), Catch::Matchers::StartsWith("ERROR: third magic byte in input .bed file") );
	}

	SECTION(".bim file reading") {
		const std::string bimFileName("../tests/ind197_397.bim");
		std::vector<std::string> locusNames{BayesicSpace::getLocusNames(bimFileName)};
		REQUIRE(locusNames.at(1)  == "14155618");
		REQUIRE(locusNames.back() == "14168708");
	}

	SECTION("Binarization and similarity groups") {
		constexpr size_t nBitsInByte{8};
		constexpr size_t nIndividuals{17};
		constexpr size_t nIndivPerBedByte{4};
		constexpr size_t nIndivPerBinByte{8};
		constexpr size_t nBedBytes{5};
		constexpr size_t nBinBytes{3};
		constexpr std::array<uint8_t, nBedBytes> bedBytes{0b11001100, 0b00011011, 0b11001100, 0b00111001, 0b00000011};
		for (uint16_t iRanIt = 0; iRanIt < N_RAN_ITERATIONS; ++iRanIt) {
			BayesicSpace::LocationWithLength bedWindow{0, bedBytes.size()};
			std::vector<char> bedByteVec{bedBytes.begin(), bedBytes.end()};
			std::vector<uint8_t> binBytes(nBinBytes, 0);
			BayesicSpace::LocationWithLength binWindow{0, nBinBytes};
			BayesicSpace::binarizeBedLocus(bedWindow, bedByteVec, nIndividuals, binWindow, binBytes);
			REQUIRE(nIndividuals >= BayesicSpace::countSetBits(binBytes) * 2);
			std::vector<BayesicSpace::HashGroup> groups;
			constexpr std::array<size_t, 3> groupSizes{7, 5, 11};
			constexpr size_t correctVGsize{86};
			groups.reserve( groupSizes.size() );
			size_t gStart{0};
			for (const auto &iGrpSize : groupSizes) {
				BayesicSpace::HashGroup currGrp;
				currGrp.locusIndexes.resize(iGrpSize);
				std::iota(currGrp.locusIndexes.begin(), currGrp.locusIndexes.end(), gStart);
				gStart                  += iGrpSize;
				currGrp.cumulativeNpairs = gStart;
				groups.emplace_back(currGrp);
			}
			std::vector<BayesicSpace::IndexedPairSimilarity> vectorizedGroups{BayesicSpace::vectorizeGroups(groups.cbegin(), groups.cend())};
			REQUIRE(vectorizedGroups.size() == correctVGsize);
			REQUIRE(std::all_of(
						vectorizedGroups.cbegin(),
						vectorizedGroups.cend(),
						[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.similarityValue == 0.0F;}
					)
			);
			REQUIRE(std::all_of(
						vectorizedGroups.cbegin(),
						vectorizedGroups.cend(),
						[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.element1ind < eachObj.element2ind;}
					)
			);
		}
	}
}

TEST_CASE("SimilarityMatrix methods work", "[SimilarityMatrix]") {
	constexpr std::array<uint32_t, 7> rowIndexes{4, 5, 5, 6, 6, 8, 8};
	constexpr std::array<uint32_t, 7> colIndexes{3, 1, 2, 2, 4, 1, 7};
	constexpr std::array<uint8_t, 7>  values{54, 81, 84, 141, 124, 199, 146};
	constexpr uint32_t nRows{9};
	constexpr uint64_t previousIdx{7};
	constexpr uint64_t previousIdxTooBig{14};

	std::vector<BayesicSpace::RowColIdx> idxPairs;
	idxPairs.reserve( rowIndexes.size() );
	size_t vecIdx{0};
	while ( vecIdx < rowIndexes.size() ) {
		BayesicSpace::RowColIdx currPair{};
		currPair.iRow = rowIndexes.at(vecIdx);
		currPair.jCol = colIndexes.at(vecIdx);
		idxPairs.emplace_back(currPair);
		++vecIdx;
	}

	// the first chunk in the overall matrix
	BayesicSpace::MatrixInitializer firstInit{};
	firstInit.previousCumulativeIndex = 0;
	firstInit.nRows                   = nRows;
	BayesicSpace::SimilarityMatrix firstMatrix(firstInit);
	constexpr size_t initialSize{28};
	REQUIRE(firstMatrix.size() == initialSize);
	vecIdx = 0;
	while ( vecIdx < values.size() ) {
		firstMatrix.append( idxPairs.at(vecIdx), values.at(vecIdx) );
		++vecIdx;
	}
	REQUIRE( firstMatrix.size() == ( initialSize + sizeof(uint32_t) * idxPairs.size() ) );

	// a later chunk in the overall matrix
	BayesicSpace::MatrixInitializer laterInit{};
	laterInit.previousCumulativeIndex = previousIdx;
	laterInit.nRows                   = nRows;
	BayesicSpace::SimilarityMatrix laterMatrix(laterInit);
	REQUIRE(laterMatrix.size() == initialSize);
	vecIdx = 0;
	while ( vecIdx < values.size() ) {
		laterMatrix.append( idxPairs.at(vecIdx), values.at(vecIdx) );
		++vecIdx;
	}
	REQUIRE( laterMatrix.size() == firstMatrix.size() );
	laterMatrix.save();

	// throwing tests
	BayesicSpace::RowColIdx wrongCombo{};
	wrongCombo.iRow = 1;
	wrongCombo.jCol = 1;
	REQUIRE_THROWS_WITH(laterMatrix.append(wrongCombo, values.at(0)), 
		Catch::Matchers::StartsWith("ERROR: row and column indexes must be different in")
	);

	wrongCombo.iRow = 0;
	REQUIRE_THROWS_WITH(laterMatrix.append(wrongCombo, values.at(0)), 
		Catch::Matchers::StartsWith("ERROR: row index must be > 0 in")
	);

	wrongCombo.iRow = nRows;
	REQUIRE_THROWS_WITH(laterMatrix.append(wrongCombo, values.at(0)), 
		Catch::Matchers::StartsWith("ERROR: row index must be smaller than the number of rows in")
	);

	laterInit.previousCumulativeIndex = previousIdxTooBig;
	BayesicSpace::SimilarityMatrix tooLateMatrix(laterInit);
	REQUIRE_THROWS_WITH(tooLateMatrix.append(idxPairs.at(0), values.at(0)), 
		Catch::Matchers::StartsWith("ERROR: new entry in front of the last element in")
	);
}

TEST_CASE("GenoTableBin methods work", "[gtBin]") {
	const std::string logFileName("../tests/binTest.log");
	const std::string inputBedName("../tests/ind197_397.bed");
	constexpr size_t nIndividuals{197};
	constexpr size_t nThreads{4};
	SECTION("Failed GenoTableBin constructors") {
		constexpr size_t smallNind{1};
		constexpr size_t largeNind{std::numeric_limits<size_t>::max() - 3};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(inputBedName, smallNind, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(inputBedName, largeNind, logFileName, nThreads),
				Catch::Matchers::ContainsSubstring("is too big to make a square relationship matrix") );
		const std::string absentFileName("../tests/noSuchFile.bed");
		const std::string noLociFile("../tests/threeByte.bed");
		const std::string wrongMagicBytes("../tests/wrongMB.bed");
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(absentFileName, nIndividuals, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: failed to open file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(noLociFile, nIndividuals, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: no genotype records in file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(wrongMagicBytes, nIndividuals, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: first magic byte in input .bed file") );
		const std::vector<int> smallMACvec(13, 0);
		const std::vector<int> emptyMACvec{};
		constexpr size_t undivNind{5};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(emptyMACvec, nIndividuals, logFileName),
				Catch::Matchers::StartsWith("ERROR: empty vector of minor allele counts") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(smallMACvec, smallNind, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(smallMACvec, undivNind, logFileName),
				Catch::Matchers::StartsWith("ERROR: length of allele count vector") );
	}
	SECTION("GenoTableBin constructors and methods with correct data") {
		BayesicSpace::GenoTableBin bedGTB(inputBedName, nIndividuals, logFileName, nThreads);
		std::vector<BayesicSpace::IndexedPairLD> bedLD{bedGTB.allJaccardLD()};
		REQUIRE(std::all_of(
					bedLD.cbegin(),
					bedLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.element1ind < eachObj.element2ind;}
				)
		);
		REQUIRE(std::all_of(
					bedLD.cbegin(),
					bedLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard <= 1.0F;}
				)
		);
		REQUIRE(std::all_of(
					bedLD.cbegin(),
					bedLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard >= 0.0F;}
				)
		);
		REQUIRE(std::all_of(
					bedLD.cbegin(),
					bedLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.rSq <= 1.0F;}
				)
		);
		REQUIRE(std::all_of(
					bedLD.cbegin(),
					bedLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.rSq >= 0.0F;}
				)
		);
		const std::string alleleCountsFile("../tests/alleleCounts.txt");
		std::fstream inAlleleCounts;
		std::string eachLine;
		inAlleleCounts.open(alleleCountsFile, std::ios::in);
		std::vector<int> macVector;
		while ( std::getline(inAlleleCounts, eachLine) ) {
			macVector.push_back( std::stoi(eachLine) );
		}
		inAlleleCounts.close();
		BayesicSpace::GenoTableBin macGTB(macVector, nIndividuals, logFileName);
		std::vector<BayesicSpace::IndexedPairLD> macLD{macGTB.allJaccardLD()};
		REQUIRE( bedLD.size() == macLD.size() );
		REQUIRE(std::all_of(
					macLD.cbegin(),
					macLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.element1ind < eachObj.element2ind;}
				)
		);
		REQUIRE(std::all_of(
					macLD.cbegin(),
					macLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard <= 1.0F;}
				)
		);
		REQUIRE(std::all_of(
					macLD.cbegin(),
					macLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard >= 0.0F;}
				)
		);
		REQUIRE(std::all_of(
					macLD.cbegin(),
					macLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.rSq <= 1.0F;}
				)
		);
		REQUIRE(std::all_of(
					macLD.cbegin(),
					macLD.cend(),
					[](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.rSq >= 0.0F;}
				)
		);
		constexpr float upperCutOff{0.9F};
		constexpr uint32_t correctNlargeLD{2400};
		uint32_t nLargeLD = std::count_if(
			bedLD.cbegin(),
			bedLD.cend(), 
			[upperCutOff](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard >= upperCutOff;}
		);
		REQUIRE(nLargeLD >= correctNlargeLD); // cannot test equality b/c of randomness
		constexpr float lowerCutOff{0.1F};
		constexpr uint32_t correctNsmallLD{48600};
		uint32_t nSmallLD = std::count_if(
			bedLD.cbegin(),
			bedLD.cend(), 
			[lowerCutOff](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard <= lowerCutOff;}
		);
		REQUIRE(nSmallLD >= correctNsmallLD); // cannot test equality b/c of randomness
		REQUIRE( nSmallLD + nLargeLD <= bedLD.size() );
		nLargeLD = std::count_if(
			macLD.cbegin(),
			macLD.cend(), 
			[upperCutOff](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard >= upperCutOff;}
		);
		REQUIRE(nLargeLD >= correctNlargeLD); // cannot test equality b/c of randomness
		nSmallLD = std::count_if(
			macLD.cbegin(),
			macLD.cend(), 
			[lowerCutOff](const BayesicSpace::IndexedPairLD &eachObj){return eachObj.jaccard <= lowerCutOff;}
		);
		REQUIRE(nSmallLD >= correctNsmallLD); // cannot test equality b/c of randomness
		REQUIRE( nSmallLD + nLargeLD <= bedLD.size() );
		// test the move constructor
		BayesicSpace::GenoTableBin movedBedGTB = std::move(bedGTB);
		std::vector<BayesicSpace::IndexedPairLD> movedBedLD{movedBedGTB.allJaccardLD()};
		REQUIRE( std::equal(
				bedLD.cbegin(), 
				bedLD.cend(), 
				movedBedLD.cbegin(), 
				[](const BayesicSpace::IndexedPairLD &first, const BayesicSpace::IndexedPairLD &second){return std::fabs(first.jaccard - second.jaccard) <= FPREC; }
			) 
		);
	}
}

TEST_CASE("GenoTableHash methods work", "[gtHash]") {
	const std::string logFileName("../tests/binTest.log");
	const std::string inputBedName("../tests/ind197_397.bed");
	constexpr size_t nIndividuals{197};
	constexpr size_t kSketches{29};
	constexpr size_t nRowsPerBand{5};
	constexpr float invKlowBound{0.034};
	constexpr float invKhighBound{0.966};
	constexpr size_t nThreads{4};
	constexpr size_t nLoci{397};
	constexpr size_t totNpairs{nLoci * (nLoci - 1) / 2};
	const std::string alleleCountsFile("../tests/alleleCounts.txt");
	std::fstream inAlleleCounts;
	std::string eachLine;
	inAlleleCounts.open(alleleCountsFile, std::ios::in);
	std::vector<int> macVector;
	while ( std::getline(inAlleleCounts, eachLine) ) {
		macVector.push_back( std::stoi(eachLine) );
	}
	inAlleleCounts.close();
	constexpr BayesicSpace::IndividualAndSketchCounts sketchParameters{nIndividuals, kSketches};
	SECTION("Failed GenoTableHash constructors") {
		constexpr size_t smallNind{1};
		constexpr size_t smallSketch{1};
		constexpr size_t kGtN{nIndividuals + 1};
		constexpr size_t lgNind{std::numeric_limits<uint16_t>::max() + 3};
		constexpr size_t largeSketch{std::numeric_limits<uint16_t>::max() + 1};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{smallNind, smallSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{nIndividuals, smallSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be at least three") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{nIndividuals, kGtN}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be smaller than the number of individuals") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{lgNind, largeSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: Number of sketches (") );
		const std::string absentFileName("../tests/noSuchFile.bed");
		const std::string noLociFile("../tests/threeByte.bed");
		const std::string wrongMagicBytes("../tests/wrongMB.bed");
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(absentFileName, sketchParameters, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: failed to open file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(noLociFile, sketchParameters, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: no genotype records in file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(wrongMagicBytes, sketchParameters, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: first magic byte in input .bed file") );
		const std::vector<int> smallMACvec(13, 0);
		const std::vector<int> emptyMACvec{};
		constexpr size_t undivNind{5};
		const std::vector<int> tooBig(2UL * std::numeric_limits<uint16_t>::max(), 0);
		constexpr size_t maxNind{2UL * std::numeric_limits<uint16_t>::max()};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(emptyMACvec, sketchParameters, logFileName),
				Catch::Matchers::StartsWith("ERROR: empty vector of minor allele counts") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(smallMACvec, BayesicSpace::IndividualAndSketchCounts{smallNind, kSketches}, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(smallMACvec, BayesicSpace::IndividualAndSketchCounts{undivNind, kSketches}, logFileName),
				Catch::Matchers::StartsWith("ERROR: length of allele count vector") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(macVector, BayesicSpace::IndividualAndSketchCounts{nIndividuals, kGtN}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be smaller than the number of individuals") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(tooBig, BayesicSpace::IndividualAndSketchCounts{maxNind, largeSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: Number of sketches (") );
	}
	SECTION("GenoTableHash .bed file constructor and methods with correct data") {
		BayesicSpace::GenoTableHash bedHSH(inputBedName, sketchParameters, nThreads, logFileName);
		std::vector<BayesicSpace::IndexedPairSimilarity> bedHLD{bedHSH.allHashLD()};
		REQUIRE(bedHLD.size() == totNpairs);
		REQUIRE(std::all_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.element1ind < eachObj.element2ind;}
				)
		);
		REQUIRE(std::all_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.similarityValue >= 0.0F;}
				)
		);
		REQUIRE(std::all_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.similarityValue <= 1.0F;}
				)
		);
		// testing discreteness of the Jaccard values that arises from the chosen sketch size
		REQUIRE(std::none_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return (eachObj.similarityValue > 0.0F) && (eachObj.similarityValue <= invKlowBound);}
				)
		);
		REQUIRE(std::none_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return (eachObj.similarityValue >= invKhighBound) && (eachObj.similarityValue < 1.0F);}
				)
		);
		const std::string tmpJacFile("../tests/tmpJac.tsv");
		BayesicSpace::InOutFileNames tmpFileGrp{};
		tmpFileGrp.outputFileName = tmpJacFile;
		tmpFileGrp.inputFileName  = "";
		constexpr size_t forcedChunks{5};
		bedHSH.allHashLD(tmpFileGrp, forcedChunks);
		std::fstream hashLDfile(tmpJacFile, std::ios::in);
		std::string line;
		std::vector<BayesicSpace::IndexedPairSimilarity> fileLD;
		fileLD.reserve(totNpairs);
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			BayesicSpace::IndexedPairSimilarity curRecord{};
			lineStream >> field;
			curRecord.element1ind = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			curRecord.element2ind = stoi(field) - 1;
			lineStream >> field;
			curRecord.similarityValue = stof(field);
			fileLD.emplace_back(curRecord);
		}
		hashLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE( fileLD.size() == bedHLD.size() );
		REQUIRE(std::equal(
				fileLD.cbegin(),
				fileLD.cend(),
				bedHLD.cbegin(),
				[](const BayesicSpace::IndexedPairSimilarity &obj1, const BayesicSpace::IndexedPairSimilarity &obj2) {
					return (obj1.element1ind == obj2.element1ind)
						&& (obj1.element2ind == obj2.element2ind)
						&& (std::fabs(obj1.similarityValue - obj2.similarityValue) <= invKlowBound);
				}
			)
		);
		std::vector<BayesicSpace::HashGroup> groups{bedHSH.makeLDgroups(nRowsPerBand)};
		auto smallestSizeIt = std::min_element(
			groups.cbegin(),
			groups.cend(),
			[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
				return grp1.locusIndexes.size() < grp2.locusIndexes.size();
			}
		);
		REQUIRE(smallestSizeIt->locusIndexes.size() >= 2);
		auto largestSizeIt = std::max_element(
			groups.cbegin(),
			groups.cend(),
			[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
				return grp1.locusIndexes.size() < grp2.locusIndexes.size();
			}
		);
		REQUIRE(largestSizeIt->locusIndexes.size() <= nIndividuals);
		REQUIRE(std::all_of(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &eachGroup) {
					return std::is_sorted( eachGroup.locusIndexes.cbegin(), eachGroup.locusIndexes.cend() );
				}
			)
		);
		REQUIRE(std::is_sorted(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2){return grp1.locusIndexes.at(0) < grp2.locusIndexes.at(0);}
			)
		);
		std::vector<BayesicSpace::IndexedPairSimilarity> groupLD{bedHSH.ldInGroups(nRowsPerBand)};
		REQUIRE(groupLD.size() <= totNpairs);
		REQUIRE(std::is_sorted(
				groupLD.cbegin(),
				groupLD.cend(),
				[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second){
					return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
				}
			)
		);
		REQUIRE(std::all_of(
				groupLD.cbegin(),
				groupLD.cend(),
				[&bedHLD](const BayesicSpace::IndexedPairSimilarity &eachGrpPair){
					auto findIt = std::find_if(
							bedHLD.cbegin(),
							bedHLD.cend(),
							[&eachGrpPair](const BayesicSpace::IndexedPairSimilarity &allPair){
								return (eachGrpPair.element1ind == allPair.element1ind) && 
										(eachGrpPair.element2ind == allPair.element2ind) &&
										(std::fabs(eachGrpPair.similarityValue - allPair.similarityValue) <= invKlowBound);
							}
						);
					return findIt != bedHLD.end();
				}
			)
		);
		bedHSH.ldInGroups(nRowsPerBand, tmpFileGrp, forcedChunks);
		std::fstream grpLDfile(tmpJacFile, std::ios::in);
		fileLD.clear();
		std::getline(grpLDfile, line);             // get rid of the header
		while ( std::getline(grpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			BayesicSpace::IndexedPairSimilarity curRecord{};
			lineStream >> field;
			curRecord.element1ind = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			curRecord.element2ind = stoi(field) - 1;
			lineStream >> field;
			curRecord.similarityValue = stof(field);
			fileLD.emplace_back(curRecord);
		}
		grpLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE( fileLD.size() >= groupLD.size() );
		std::sort(
			fileLD.begin(),
			fileLD.end(),
			[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second) {
				return (first.element1ind == second.element1ind ? first.element2ind < second.element2ind : first.element1ind < second.element1ind);
			}
		);
		auto lastUniqueIt = std::unique(
			fileLD.begin(),
			fileLD.end(),
			[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second) {
				return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
			}
		);
		fileLD.erase( lastUniqueIt, fileLD.end() );
		fileLD.shrink_to_fit();
		REQUIRE(std::equal(
				fileLD.cbegin(),
				fileLD.cend(),
				groupLD.cbegin(),
				[](const BayesicSpace::IndexedPairSimilarity &obj1, const BayesicSpace::IndexedPairSimilarity &obj2) {
					return (obj1.element1ind == obj2.element1ind)
						&& (obj1.element2ind == obj2.element2ind)
						&& (std::fabs(obj1.similarityValue - obj2.similarityValue) <= invKlowBound);
				}
			)
		);
	}
	SECTION("GenoTableHash mac vector constructor and methods with correct data") {
		BayesicSpace::GenoTableHash vecHSH(macVector, sketchParameters, nThreads, logFileName);
		std::vector<BayesicSpace::IndexedPairSimilarity> bedHLD{vecHSH.allHashLD()};
		REQUIRE(bedHLD.size() == totNpairs);
		REQUIRE(std::all_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.element1ind < eachObj.element2ind;}
				)
		);
		REQUIRE(std::all_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.similarityValue >= 0.0F;}
				)
		);
		REQUIRE(std::all_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.similarityValue <= 1.0F;}
				)
		);
		// testing discreteness of the Jaccard values that arises from the chosen sketch size
		REQUIRE(std::none_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return (eachObj.similarityValue > 0.0F) && (eachObj.similarityValue <= invKlowBound);}
				)
		);
		REQUIRE(std::none_of(
					bedHLD.cbegin(),
					bedHLD.cend(),
					[](const BayesicSpace::IndexedPairSimilarity &eachObj){return (eachObj.similarityValue >= invKhighBound) && (eachObj.similarityValue < 1.0F);}
				)
		);
		const std::string tmpJacFile("../tests/tmpJac.tsv");
		BayesicSpace::InOutFileNames tmpFileGrp{};
		tmpFileGrp.outputFileName = tmpJacFile;
		tmpFileGrp.inputFileName  = "";
		constexpr size_t forcedChunks{5};
		vecHSH.allHashLD(tmpFileGrp, forcedChunks);
		std::fstream hashLDfile(tmpJacFile, std::ios::in);
		std::string line;
		std::vector<BayesicSpace::IndexedPairSimilarity> fileLD;
		fileLD.reserve(totNpairs);
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			BayesicSpace::IndexedPairSimilarity curRecord{};
			lineStream >> field;
			curRecord.element1ind = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			curRecord.element2ind = stoi(field) - 1;
			lineStream >> field;
			curRecord.similarityValue = stof(field);
			fileLD.emplace_back(curRecord);
		}
		hashLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE( fileLD.size() == bedHLD.size() );
		REQUIRE(std::equal(
				fileLD.cbegin(),
				fileLD.cend(),
				bedHLD.cbegin(),
				[](const BayesicSpace::IndexedPairSimilarity &obj1, const BayesicSpace::IndexedPairSimilarity &obj2) {
					return (obj1.element1ind == obj2.element1ind)
						&& (obj1.element2ind == obj2.element2ind)
						&& (std::fabs(obj1.similarityValue - obj2.similarityValue) <= invKlowBound);
				}
			)
		);
		std::vector<BayesicSpace::HashGroup> groups{vecHSH.makeLDgroups(nRowsPerBand)};
		auto smallestSizeIt = std::min_element(
			groups.cbegin(),
			groups.cend(),
			[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
				return grp1.locusIndexes.size() < grp2.locusIndexes.size();
			}
		);
		REQUIRE(smallestSizeIt->locusIndexes.size() >= 2);
		auto largestSizeIt = std::max_element(
			groups.cbegin(),
			groups.cend(),
			[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
				return grp1.locusIndexes.size() < grp2.locusIndexes.size();
			}
		);
		REQUIRE(largestSizeIt->locusIndexes.size() <= nIndividuals);
		REQUIRE(std::all_of(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &eachGroup) {
					return std::is_sorted( eachGroup.locusIndexes.cbegin(), eachGroup.locusIndexes.cend() );
				}
			)
		);
		REQUIRE(std::is_sorted(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2){return grp1.locusIndexes.at(0) < grp2.locusIndexes.at(0);}
			)
		);
		std::vector<BayesicSpace::IndexedPairSimilarity> groupLD{vecHSH.ldInGroups(nRowsPerBand)};
		REQUIRE(groupLD.size() <= totNpairs);
		REQUIRE(std::is_sorted(
				groupLD.cbegin(),
				groupLD.cend(),
				[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second){
					return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
				}
			)
		);
		REQUIRE(std::all_of(
				groupLD.cbegin(),
				groupLD.cend(),
				[&bedHLD](const BayesicSpace::IndexedPairSimilarity &eachGrpPair){
					auto findIt = std::find_if(
							bedHLD.cbegin(),
							bedHLD.cend(),
							[&eachGrpPair](const BayesicSpace::IndexedPairSimilarity &allPair){
								return (eachGrpPair.element1ind == allPair.element1ind) && 
										(eachGrpPair.element2ind == allPair.element2ind) &&
										(std::fabs(eachGrpPair.similarityValue - allPair.similarityValue) <= invKlowBound);
							}
						);
					return findIt != bedHLD.end();
				}
			)
		);
		vecHSH.ldInGroups(nRowsPerBand, tmpFileGrp, forcedChunks);
		std::fstream grpLDfile(tmpJacFile, std::ios::in);
		fileLD.clear();
		std::getline(grpLDfile, line);             // get rid of the header
		while ( std::getline(grpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			BayesicSpace::IndexedPairSimilarity curRecord{};
			lineStream >> field;
			curRecord.element1ind = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			curRecord.element2ind = stoi(field) - 1;
			lineStream >> field;
			curRecord.similarityValue = stof(field);
			fileLD.emplace_back(curRecord);
		}
		grpLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE( fileLD.size() >= groupLD.size() );
		std::sort(fileLD.begin(), fileLD.end(),
					[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second) {
						return (first.element1ind == second.element1ind ? first.element2ind < second.element2ind : first.element1ind < second.element1ind);
					}
				);
		auto lastUniqueIt = std::unique(fileLD.begin(), fileLD.end(),
					[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second) {
						return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
					}
				);
		fileLD.erase( lastUniqueIt, fileLD.end() );
		fileLD.shrink_to_fit();
		REQUIRE(std::equal(
				fileLD.cbegin(),
				fileLD.cend(),
				groupLD.cbegin(),
				[](const BayesicSpace::IndexedPairSimilarity &obj1, const BayesicSpace::IndexedPairSimilarity &obj2) {
					return (obj1.element1ind == obj2.element1ind)
						&& (obj1.element2ind == obj2.element2ind)
						&& (std::fabs(obj1.similarityValue - obj2.similarityValue) <= invKlowBound);
				}
			)
		);
	}
}
