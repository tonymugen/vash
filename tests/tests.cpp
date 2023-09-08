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

#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <fstream>

#include <iostream>
#include <bitset>

#include "random.hpp"
#include "gvarHash.hpp"
#include "vashFunctions.hpp"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

// Number of times tests of random events will be run
static constexpr uint16_t N_RAN_ITERATIONS{10};

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
	constexpr std::array<uint16_t, 13> array16bit{
		1256, 2117, 2866, 7434,
		11737, 16256, 22236, 39883,
		40023, 46299, 58123, 58167, 62187
	};
	std::vector<uint16_t> vector16bit( array16bit.begin(), array16bit.end() );
	constexpr BayesicSpace::LocationWithLength keyWindow{4, 5};
	constexpr uint32_t correct16bitHash{3760365877};
	const uint32_t v16bitHash{BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed)};
	SECTION("MurMurHash correctness tests") {
		REQUIRE(arrayMMhash   == correctArrayHash);
		REQUIRE(idxVectorHash == correctIdxVectorHash);
		REQUIRE(v16bitHash    == correct16bitHash);
	}
	SECTION("MurMurHash sensitivity tests") {
		constexpr uint32_t mmHashSeed2{2153025619};
		constexpr std::array<uint32_t, BayesicSpace::SIZE_OF_SIZET> arrayKey2{335636696, 4242517348};
		REQUIRE(BayesicSpace::murMurHash(arrayKey, mmHashSeed2)  != correctArrayHash);
		REQUIRE(BayesicSpace::murMurHash(arrayKey2, mmHashSeed)  != correctArrayHash);
		REQUIRE(BayesicSpace::murMurHash(idxVector, mmHashSeed2) != correctIdxVectorHash);
		idxVector.at(1)++;
		REQUIRE(BayesicSpace::murMurHash(idxVector, mmHashSeed)  != correctIdxVectorHash);
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

	SECTION("Magic byte testing") {
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> correctBedBytes{0x6c, 0x1b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes1{0x6d, 0x1b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes2{0x6c, 0x0b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes3{0x6c, 0x1b, 0x11};
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
			BayesicSpace::RanDraw prng;
			BayesicSpace::LocationWithLength bedWindow{0, bedBytes.size()};
			std::vector<char> bedByteVec{bedBytes.begin(), bedBytes.end()};
			std::vector<uint8_t> binBytes(nBinBytes, 0);
			BayesicSpace::LocationWithLength binWindow{0, nBinBytes};
			BayesicSpace::binarizeBedLocus(bedWindow, bedByteVec, nIndividuals, prng, binWindow, binBytes);
			REQUIRE(nIndividuals >= BayesicSpace::countSetBits(binBytes) * 2);
			std::vector< std::vector<uint32_t> > groups;
			constexpr std::array<size_t, 3> groupSizes{7, 5, 11};
			constexpr size_t correctVGsize{86};
			groups.reserve( groupSizes.size() );
			for (const auto &iGrpSize : groupSizes) {
				groups.emplace_back(iGrpSize);
			}
			std::vector<BayesicSpace::IndexedPairSimilarity> vectorizedGroups{BayesicSpace::vectorizeGroups(0, groups.begin(), groups.end())};
			REQUIRE(vectorizedGroups.size() == correctVGsize);
			REQUIRE(std::all_of(
						vectorizedGroups.cbegin(),
						vectorizedGroups.cend(),
						[](const BayesicSpace::IndexedPairSimilarity &eachObj){return eachObj.similarityValue == 0.0F;}
					)
			);
			REQUIRE(std::is_sorted(
						vectorizedGroups.cbegin(),
						vectorizedGroups.cend(), 
						[](BayesicSpace::IndexedPairSimilarity obj1, BayesicSpace::IndexedPairSimilarity obj2){return obj1.groupID < obj2.groupID;}
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

TEST_CASE("GenoTableBin methods work", "[gtBin]") {
	const std::string logFileName("../tests/binTest.log");
	const std::string inputBedName("../tests/ind197_397.bed");
	constexpr size_t nIndividuals{197};
	constexpr size_t nThreads{4};
	constexpr size_t oneThread{1};
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
		const std::string alleleCountsFile("../tests/alleleCounts.txt");
		std::fstream inAlleleCounts;
		std::string eachLine;
		inAlleleCounts.open(alleleCountsFile, std::ios::in);
		std::vector<int> macVector;
		while ( std::getline(inAlleleCounts, eachLine) ) {
			macVector.push_back( std::stoi(eachLine) );
		}
		inAlleleCounts.close();
		BayesicSpace::GenoTableBin(macVector, nIndividuals, logFileName);
	}
}
