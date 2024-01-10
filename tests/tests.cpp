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
#include <iterator>
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
#include <chrono>

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
				[](const std::pair<size_t, size_t> &eachRange) {return eachRange.first <= eachRange.second;}
			)
		);
		REQUIRE(std::equal(
				threadRanges.cbegin(),
				threadRanges.cend(),
				correctRanges.cbegin(),
				[](const std::pair<size_t, size_t> &pair1, const std::pair<size_t, size_t>&pair2) {
					return (pair1.first == pair2.first) && (pair1.second == pair2.second);
				}
			)
		);

		constexpr size_t nVecElements{35};
		constexpr size_t nChunks{4};
		constexpr std::array<size_t, nChunks> correctChunkSizes{9, 9, 9, 8};
		const std::vector<size_t> chunkSizes{BayesicSpace::makeChunkSizes(nVecElements, nChunks)};
		REQUIRE(std::equal(
				chunkSizes.cbegin(),
				chunkSizes.cend(),
				correctChunkSizes.cbegin()
			) 
		);
		constexpr size_t smallNelements{3};
		const std::vector<size_t> smallChunkSizes{BayesicSpace::makeChunkSizes(smallNelements, nChunks)};
		constexpr std::array<size_t, nChunks> correctSmallChunkSizes{1, 1, 1, 0};
		REQUIRE(std::equal(
				smallChunkSizes.cbegin(),
				smallChunkSizes.cend(),
				correctSmallChunkSizes.cbegin()
			)
		);

		constexpr std::array<uint32_t, nChunks> correctRowStarts{1, 4, 6, 7};
		constexpr std::array<uint32_t, nChunks> correctRowEnds{4, 6, 7, 8};
		constexpr std::array<uint32_t, nChunks> correctColStarts{0, 3, 3, 6};
		constexpr std::array<uint32_t, nChunks> correctColEnds{3, 3, 6, 7};
		std::vector< std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> > correctRowColPairs;
		size_t iChunk{0};
		while (iChunk < nChunks) {
			std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> tmpPair;
			tmpPair.first.iRow  = correctRowStarts.at(iChunk);
			tmpPair.first.jCol  = correctColStarts.at(iChunk);
			tmpPair.second.iRow = correctRowEnds.at(iChunk);
			tmpPair.second.jCol = correctColEnds.at(iChunk);
			correctRowColPairs.emplace_back(tmpPair);
			++iChunk;
		}
		BayesicSpace::LocationWithLength startAndNelements{};
		startAndNelements.start  = 0;
		startAndNelements.length = nVecElements;
		std::vector< std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> > rowColPairs{BayesicSpace::makeChunkRanges(startAndNelements, nChunks)};
		REQUIRE(std::equal(
				rowColPairs.cbegin(),
				rowColPairs.cend(),
				correctRowColPairs.cbegin(),
				[](const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair1,
							const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair2) {
					return  (pair1.first.iRow  == pair2.first.iRow) &&
							(pair1.first.jCol  == pair2.first.jCol) &&
							(pair1.second.iRow == pair2.second.iRow) &&
							(pair1.second.jCol == pair2.second.jCol);
				}
			)
		);
		BayesicSpace::LocationWithLength smallStartAndNelements{};
		smallStartAndNelements.start  = 0;
		smallStartAndNelements.length = smallNelements;
		std::vector< std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> > smallRowColPairs{BayesicSpace::makeChunkRanges(smallStartAndNelements, nChunks)};
		constexpr std::array<uint32_t, nChunks> smallCorrectRowStarts{1, 2, 2, 3};
		constexpr std::array<uint32_t, nChunks> smallCorrectRowEnds{2, 2, 3, 3};
		constexpr std::array<uint32_t, nChunks> smallCorrectColStarts{0, 0, 1, 0};
		constexpr std::array<uint32_t, nChunks> smallCorrectColEnds{0, 1, 0, 0};
		iChunk = 0;
		while (iChunk < nChunks) {
			std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> tmpPair;
			tmpPair.first.iRow  = smallCorrectRowStarts.at(iChunk);
			tmpPair.first.jCol  = smallCorrectColStarts.at(iChunk);
			tmpPair.second.iRow = smallCorrectRowEnds.at(iChunk);
			tmpPair.second.jCol = smallCorrectColEnds.at(iChunk);
			correctRowColPairs.at(iChunk) = std::move(tmpPair);
			++iChunk;
		}
		REQUIRE(std::equal(
				smallRowColPairs.cbegin(),
				smallRowColPairs.cend(),
				correctRowColPairs.cbegin(),
				[](const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair1,
							const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair2) {
					return  (pair1.first.iRow  == pair2.first.iRow) &&
							(pair1.first.jCol  == pair2.first.jCol) &&
							(pair1.second.iRow == pair2.second.iRow) &&
							(pair1.second.jCol == pair2.second.jCol);
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
				[](const BayesicSpace::IndexedPairSimilarity &obj) {return obj.element1ind < obj.element2ind;}
			)
		);
		REQUIRE(std::all_of(
				pairSegment.cbegin(),
				pairSegment.cend(),
				[&allN](const BayesicSpace::IndexedPairSimilarity &obj) {return (obj.element1ind < allN) || (obj.element2ind < allN);}
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
						[](const BayesicSpace::IndexedPairSimilarity &eachObj) {return eachObj.similarityValue == 0.0F;}
					)
			);
			REQUIRE(std::all_of(
						vectorizedGroups.cbegin(),
						vectorizedGroups.cend(),
						[](const BayesicSpace::IndexedPairSimilarity &eachObj) {return eachObj.element1ind < eachObj.element2ind;}
					)
			);
		}
	}
}

TEST_CASE("SimilarityMatrix methods work", "[SimilarityMatrix]") {
	constexpr std::array<uint32_t, 7> rowIndexes{4, 5, 5, 6, 6, 8, 8};
	constexpr std::array<uint32_t, 7> colIndexes{3, 1, 2, 2, 4, 1, 7};
	// for Jaccard calculation
	constexpr std::array<uint64_t, 7> nIsect{21, 81, 164, 551, 16, 8, 2};
	constexpr std::array<uint64_t, 7> nUnion{101, 255, 502, 1001, 35, 11, 5};
	// to test insertion
	constexpr std::array<uint32_t, 3> addRowIndexes{8, 3, 5};
	constexpr std::array<uint32_t, 3> addColIndexes{6, 0, 1};
	constexpr std::array<uint64_t, 3> addIsect{89, 44, 6};
	constexpr std::array<uint64_t, 3> addUnion{234, 51, 13};

	constexpr std::array<uint32_t, 9> correctRowIndexes{3, 4, 5, 5, 6, 6, 8, 8, 8};
	constexpr std::array<uint32_t, 9> correctColIndexes{0, 3, 1, 2, 2, 4, 1, 6, 7};
	constexpr std::array<float,    9> correctFloatValues{0.8627, 0.2078, 0.3176, 0.3255, 0.5490, 0.4549, 0.7255, 0.3765, 0.4000};
	constexpr uint64_t previousIdx{7};
	constexpr size_t nThreads{2};

	SECTION("Auxiliary functions") {
		constexpr uint32_t byteSize{8};
		constexpr std::array<uint32_t, 4> diffIdxArray{(2 << byteSize), (2 << byteSize), (1 << byteSize), (5 << byteSize)};
		const std::vector<uint32_t> diffIdx( diffIdxArray.cbegin(), diffIdxArray.cend() );
		constexpr uint64_t correctFullIdx{9};
		REQUIRE(BayesicSpace::recoverFullVIdx(diffIdx.cbegin(), previousIdx) == correctFullIdx);
		uint64_t previousIdxMutable{previousIdx};
		BayesicSpace::RowColIdx rowColValues{BayesicSpace::recoverRCindexes(diffIdx.cbegin(), previousIdxMutable)};
		REQUIRE( rowColValues.iRow == rowIndexes.at(0) );
		REQUIRE( rowColValues.jCol == colIndexes.at(0) );
		rowColValues = BayesicSpace::recoverRCindexes(correctFullIdx);
		REQUIRE( rowColValues.iRow == rowIndexes.at(0) );
		REQUIRE( rowColValues.jCol == colIndexes.at(0) );
	}
	SECTION("SimilarityMatrix methods") {
		std::array<BayesicSpace::RowColIdx, rowIndexes.size()> idxPairs{};
		std::array<BayesicSpace::JaccardPair, rowIndexes.size()> jaccPairs{};
		size_t vecIdx{0};
		while ( vecIdx < rowIndexes.size() ) {
			idxPairs.at(vecIdx).iRow = rowIndexes.at(vecIdx);
			idxPairs.at(vecIdx).jCol = colIndexes.at(vecIdx);

			jaccPairs.at(vecIdx).nIntersect = nIsect.at(vecIdx);
			jaccPairs.at(vecIdx).nUnion     = nUnion.at(vecIdx);
			++vecIdx;
		}

		std::array<BayesicSpace::RowColIdx, addRowIndexes.size()> addIdxPairs{};
		std::array<BayesicSpace::JaccardPair, addRowIndexes.size()> addJaccPairs{};
		vecIdx = 0;
		while ( vecIdx < addRowIndexes.size() ) {
			addIdxPairs.at(vecIdx).iRow = addRowIndexes.at(vecIdx);
			addIdxPairs.at(vecIdx).jCol = addColIndexes.at(vecIdx);

			addJaccPairs.at(vecIdx).nIntersect = addIsect.at(vecIdx);
			addJaccPairs.at(vecIdx).nUnion     = addUnion.at(vecIdx);
			++vecIdx;
		}

		BayesicSpace::SimilarityMatrix testMatrix;
		constexpr size_t initialSize{2 * sizeof(uint64_t)};
		REQUIRE(testMatrix.size() == initialSize);
		vecIdx = 0;
		while ( vecIdx < rowIndexes.size() ) {
			testMatrix.insert( idxPairs.at(vecIdx), jaccPairs.at(vecIdx) );
			++vecIdx;
		}
		// test the last value insertion bypass
		testMatrix.insert( idxPairs.back(), jaccPairs.back() );
		REQUIRE( testMatrix.size() == ( initialSize + sizeof(uint32_t) * idxPairs.size() ) );
		vecIdx = 0;
		while ( vecIdx < addRowIndexes.size() ) {
			testMatrix.insert( addIdxPairs.at(vecIdx), addJaccPairs.at(vecIdx) );
			++vecIdx;
		}
		REQUIRE( testMatrix.size() == ( initialSize + sizeof(uint32_t) * correctFloatValues.size() ) );

		// test file save
		const std::string outputFileName("../tests/smallSimilarityMatrix.tsv");
		testMatrix.save(outputFileName, nThreads);
		std::fstream testSMoutfile(outputFileName, std::ios::in);
		std::string line;
		std::array<uint32_t, correctFloatValues.size()> rowsFromFile{};
		std::array<uint32_t, correctFloatValues.size()> colsFromFile{};
		std::array<float,    correctFloatValues.size()> floatsFromFile{};
		size_t arrayIdx{0};
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrayIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrayIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrayIdx) = stof(field);
			++arrayIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correctRowIndexes.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correctColIndexes.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correctFloatValues.cbegin() ) );

		const std::string bimFile("../tests/ind197_397.bim");
		testMatrix.save(outputFileName, nThreads, bimFile);
		arrayIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			floatsFromFile.at(arrayIdx) = stof(field);
			++arrayIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correctFloatValues.cbegin() ) );

		const std::string smallFile("../tests/small.bim");
		REQUIRE_THROWS_WITH(testMatrix.save(outputFileName, nThreads, smallFile),
			Catch::Matchers::StartsWith("ERROR: number of rows exceeds locus name count in") );

		// values that require padding
		constexpr std::array<uint32_t, 3> largeIdxRows{20001, 9999, 40005};
		constexpr std::array<uint32_t, 3> largeIdxCols{10001, 9899, 30005};
		constexpr std::array<uint64_t, 3> largeIsect{3, 246, 39};
		constexpr std::array<uint64_t, 3> largeUnion{255, 255, 255};
		constexpr std::array<float,    3> largeFloats{0.9647F, 0.0118F, 0.1529F};
		constexpr size_t correctLargeSize{46};
		std::array<BayesicSpace::RowColIdx, largeIdxRows.size()> largeIdxPairs{};
		std::array<BayesicSpace::JaccardPair, largeIdxRows.size()> largeJaccPairs{};
		vecIdx = 0;
		while ( vecIdx < largeIdxRows.size() ) {
			largeIdxPairs.at(vecIdx).iRow = largeIdxRows.at(vecIdx);
			largeIdxPairs.at(vecIdx).jCol = largeIdxCols.at(vecIdx);

			largeJaccPairs.at(vecIdx).nIntersect = largeIsect.at(vecIdx);
			largeJaccPairs.at(vecIdx).nUnion     = largeUnion.at(vecIdx);
			++vecIdx;
		}
		vecIdx = 0;
		BayesicSpace::SimilarityMatrix largeMatrix;
		while ( vecIdx < largeIdxRows.size() ) {
			largeMatrix.insert( largeIdxPairs.at(vecIdx), largeJaccPairs.at(vecIdx) );
			++vecIdx;
		}
		REQUIRE(largeMatrix.size() == initialSize + sizeof(uint32_t) * correctLargeSize);

		largeMatrix.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		std::array<uint32_t, largeIdxRows.size()> lgRowsFromFile{};
		std::array<uint32_t, largeIdxRows.size()> lgColsFromFile{};
		std::array<float,    largeIdxRows.size()> lgFloatsFromFile{};
		arrayIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lgRowsFromFile.at(arrayIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			lgColsFromFile.at(arrayIdx) = stoi(field) - 1;
			lineStream >> field;
			lgFloatsFromFile.at(arrayIdx) = stof(field);
			++arrayIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		constexpr std::array<uint32_t, 3> correctLargeIdxRows{9999, 20001, 40005};
		constexpr std::array<uint32_t, 3> correctLargeIdxCols{9899, 10001, 30005};
		REQUIRE(std::equal(lgRowsFromFile.cbegin(),   lgRowsFromFile.cend(),   correctLargeIdxRows.cbegin()));
		REQUIRE(std::equal(lgColsFromFile.cbegin(),   lgColsFromFile.cend(),   correctLargeIdxCols.cbegin()));
		REQUIRE(std::equal(lgFloatsFromFile.cbegin(), lgFloatsFromFile.cend(), largeFloats.cbegin()));

		// throwing tests
		BayesicSpace::RowColIdx wrongCombo{};
		BayesicSpace::SimilarityMatrix wrongMatrix;
		wrongCombo.iRow = 1;
		wrongCombo.jCol = 1;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(wrongCombo, jaccPairs.at(0)), 
			Catch::Matchers::StartsWith("ERROR: row and column indexes must be different in")
		);

		wrongCombo.iRow = 0;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(wrongCombo, jaccPairs.at(0)), 
			Catch::Matchers::StartsWith("ERROR: row index must be non-zero in")
		);
		BayesicSpace::JaccardPair wrongPair{};
		wrongPair.nIntersect = 0;
		wrongPair.nUnion     = 0;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(idxPairs.at(0), wrongPair), 
			Catch::Matchers::StartsWith("ERROR: union count cannot be 0 in")
		);
		wrongPair.nIntersect = 2;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(idxPairs.at(0), wrongPair), 
			Catch::Matchers::StartsWith("ERROR: intersection count cannot be larger than the union count in")
		);
	}
	SECTION("SimilarityMatrix merge") {
		constexpr std::array<uint32_t, 5> rowIndexes1{3, 4, 5, 6, 6};
		constexpr std::array<uint32_t, 5> colIndexes1{0, 3, 2, 2, 4};
		constexpr std::array<uint64_t, 5> nIsect1{87, 8, 77, 68, 96};
		constexpr std::array<uint64_t, 5> nUnion1{255, 255, 255, 255, 255};
		constexpr std::array<uint32_t, 5> rowIndexes2{5, 6, 8, 8, 8};
		constexpr std::array<uint32_t, 5> colIndexes2{1, 5, 1, 6, 7};
		constexpr std::array<uint64_t, 5> nIsect2{217, 160, 228, 176, 167};
		constexpr std::array<uint64_t, 5> nUnion2{255, 255, 255, 255, 255};
		constexpr std::array<uint32_t, 3> rowIndexes3{4, 5, 6};
		constexpr std::array<uint32_t, 3> colIndexes3{1, 0, 1};
		constexpr std::array<uint64_t, 3> nIsect3{140, 138, 136};
		constexpr std::array<uint64_t, 3> nUnion3{255, 255, 255};

		constexpr size_t nThreads{2};
		BayesicSpace::SimilarityMatrix matrix1;
		size_t arrIdx{0};
		while ( arrIdx < rowIndexes1.size() ) {
			BayesicSpace::RowColIdx tmp{};
			tmp.iRow = rowIndexes1.at(arrIdx);
			tmp.jCol = colIndexes1.at(arrIdx);
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion1.at(arrIdx);
			tmpJP.nIntersect = nIsect1.at(arrIdx);
			matrix1.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix matrix2;
		arrIdx = 0;
		while ( arrIdx < rowIndexes2.size() ) {
			BayesicSpace::RowColIdx tmp{};
			tmp.iRow = rowIndexes2.at(arrIdx);
			tmp.jCol = colIndexes2.at(arrIdx);
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion2.at(arrIdx);
			tmpJP.nIntersect = nIsect2.at(arrIdx);
			matrix2.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix matrix3;
		arrIdx = 0;
		while ( arrIdx < rowIndexes3.size() ) {
			BayesicSpace::RowColIdx tmp{};
			tmp.iRow = rowIndexes3.at(arrIdx);
			tmp.jCol = colIndexes3.at(arrIdx);
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion3.at(arrIdx);
			tmpJP.nIntersect = nIsect3.at(arrIdx);
			matrix3.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix tmp1 = matrix1;
		BayesicSpace::SimilarityMatrix tmp2 = matrix2;
		tmp1.merge( std::move(tmp2) );
		const std::string outputFileName("../tests/mergeMatrix.tsv");
		tmp1.save(outputFileName, nThreads);
		constexpr std::array<uint32_t, 10> correct12mergeRow{3, 4, 5, 5, 6, 6, 6, 8, 8, 8};
		constexpr std::array<uint32_t, 10> correct12mergeCol{0, 3, 1, 2, 2, 4, 5, 1, 6, 7};
		constexpr std::array<float,    10> correct12mergeValues{0.3412, 0.0314, 0.8510, 0.3020, 0.2667, 0.3765, 0.6275, 0.8941, 0.6902, 0.6549};

		std::fstream testSMoutfile(outputFileName, std::ios::in);
		std::string line;
		std::array<uint32_t, correct12mergeValues.size()> rowsFromFile{};
		std::array<uint32_t, correct12mergeValues.size()> colsFromFile{};
		std::array<float,    correct12mergeValues.size()> floatsFromFile{};
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// merge of a matrix with identical tail
		tmp2 = matrix2;
		tmp1.merge( std::move(tmp2) );
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// merge of an identical matrix
		tmp2 = tmp1;
		tmp1.merge( std::move(tmp2) );
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// merge of an empty matrix
		tmp1.merge( std::move(tmp2) );
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		tmp2.merge( std::move(tmp1) );
		tmp2.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// matrix with smaller indexes as second
		tmp1 = matrix1;
		matrix2.merge( std::move(tmp1) );
		matrix2.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// a matrix completely within another
		constexpr std::array<uint32_t, 8> correct13mergeRow{3, 4, 4, 5, 5, 6, 6, 6};
		constexpr std::array<uint32_t, 8> correct13mergeCol{0, 1, 3, 0, 2, 1, 2, 4};
		constexpr std::array<float,    8> correct13mergeValues{0.3412, 0.5490, 0.0314, 0.5412, 0.3020, 0.5333, 0.2667, 0.3765};
		tmp1 = matrix1;
		tmp2 = matrix3;
		tmp1.merge( std::move(tmp2) );
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend() - 2,   correct13mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend() - 2,   correct13mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend() - 2, correct13mergeValues.cbegin() ) );

		matrix3.merge( std::move(matrix1) );
		matrix3.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend() - 2,   correct13mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend() - 2,   correct13mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend() - 2, correct13mergeValues.cbegin() ) );

		constexpr std::array<uint32_t, 5> rowIndexes4{3, 4, 5, 5, 6};
		constexpr std::array<uint32_t, 5> colIndexes4{0, 3, 1, 2, 2};
		constexpr std::array<uint64_t, 5> nIsect4{87, 8, 217, 77, 68};
		constexpr std::array<uint64_t, 5> nUnion4{255, 255, 255, 255, 255};
		constexpr std::array<uint32_t, 5> rowIndexes5{6, 6, 8, 8, 8};
		constexpr std::array<uint32_t, 5> colIndexes5{4, 5, 1, 6, 7};
		constexpr std::array<uint64_t, 5> nIsect5{96, 160, 228, 176, 167};
		constexpr std::array<uint64_t, 5> nUnion5{255, 255, 255, 255, 255};

		BayesicSpace::SimilarityMatrix matrix4;
		arrIdx = 0;
		while ( arrIdx < rowIndexes4.size() ) {
			BayesicSpace::RowColIdx tmp{};
			tmp.iRow = rowIndexes4.at(arrIdx);
			tmp.jCol = colIndexes4.at(arrIdx);
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion4.at(arrIdx);
			tmpJP.nIntersect = nIsect4.at(arrIdx);
			matrix4.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix matrix5;
		arrIdx = 0;
		while ( arrIdx < rowIndexes5.size() ) {
			BayesicSpace::RowColIdx tmp{};
			tmp.iRow = rowIndexes5.at(arrIdx);
			tmp.jCol = colIndexes5.at(arrIdx);
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion5.at(arrIdx);
			tmpJP.nIntersect = nIsect5.at(arrIdx);
			matrix5.insert(tmp, tmpJP);
			++arrIdx;
		}
		matrix4.merge( std::move(matrix5) );
		matrix4.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// indexes far from the current (require padding)
		constexpr std::array<uint32_t, 3> largeIdxRows{9999, 20001, 40005};
		constexpr std::array<uint32_t, 3> largeIdxCols{9899, 10001, 30005};
		constexpr std::array<uint64_t, 3> largeNisect{246, 3, 39};
		constexpr std::array<uint64_t, 3> largeNunion{255, 255, 255};
		BayesicSpace::SimilarityMatrix matrixFar;
		arrIdx = 0;
		while ( arrIdx < largeIdxRows.size() ) {
			BayesicSpace::RowColIdx tmp{};
			tmp.iRow = largeIdxRows.at(arrIdx);
			tmp.jCol = largeIdxCols.at(arrIdx);
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = largeNunion.at(arrIdx);
			tmpJP.nIntersect = largeNisect.at(arrIdx);
			matrixFar.insert(tmp, tmpJP);
			++arrIdx;
		}
		matrix4.merge( std::move(matrixFar) );
		matrix4.save(outputFileName, nThreads);
		constexpr std::array<uint32_t, 13> correctFarMergeRow{3, 4, 5, 5, 6, 6, 6, 8, 8, 8, 9999, 20001, 40005};
		constexpr std::array<uint32_t, 13> correctFarMergeCol{0, 3, 1, 2, 2, 4, 5, 1, 6, 7, 9899, 10001, 30005};
		constexpr std::array<float,    13> correctFarMergeValues{0.3412, 0.0314, 0.8510, 0.3020, 0.2667, 0.3765, 0.6275,
																	0.8941, 0.6902, 0.6549, 0.9647, 0.0118, 0.1529};
		testSMoutfile.open(outputFileName, std::ios::in);
		std::array<uint32_t, correctFarMergeValues.size()> rowsFromFileFar{};
		std::array<uint32_t, correctFarMergeValues.size()> colsFromFileFar{};
		std::array<float,    correctFarMergeValues.size()> floatsFromFileFar{};
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFileFar.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFileFar.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFileFar.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFileFar.cbegin(),   rowsFromFileFar.cend(),   correctFarMergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFileFar.cbegin(),   colsFromFileFar.cend(),   correctFarMergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFileFar.cbegin(), floatsFromFileFar.cend(), correctFarMergeValues.cbegin() ) );
	}
}

TEST_CASE("GenoTableBin methods work", "[gtBin]") {
	const std::string logFileName("../tests/binTest.log");
	const std::string inputBedName("../tests/ind197_397.bed");
	constexpr uint32_t nIndividuals{197};
	constexpr size_t nThreads{4};
	SECTION("Failed GenoTableBin constructors") {
		constexpr size_t smallNind{1};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(inputBedName, smallNind, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
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
		constexpr size_t nChunks{3};
		BayesicSpace::GenoTableBin bedGTB(inputBedName, nIndividuals, logFileName, nThreads);
		const std::string ldFileName("../tests/tmpLDfile.tsv");
		std::fstream tmpLDfile;
		std::string line;
		BayesicSpace::InOutFileNames outAndBim{};
		outAndBim.outputFileName = ldFileName;
		bedGTB.allJaccardLD(outAndBim, nChunks);

		tmpLDfile.open(ldFileName, std::ios::in);
		std::vector<float> jaccValues;
		std::getline(tmpLDfile, line); // header
		while ( std::getline(tmpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.push_back( stof(field) );
		}
		tmpLDfile.close();
		std::remove( ldFileName.c_str() ); // NOLINT
		constexpr float upperCutOff{0.9F};
		constexpr uint32_t correctNlargeLD{2400};
		constexpr float lowerCutOff{0.1F};
		constexpr uint32_t correctNsmallLD{48000};
		uint32_t nLargeLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[upperCutOff](float eachValue) {return eachValue >= upperCutOff;}
		);
		REQUIRE(nLargeLD >= correctNlargeLD); // cannot test equality b/c of randomness
		uint32_t nSmallLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[lowerCutOff](float eachValue) {return eachValue <= lowerCutOff;}
		);
		REQUIRE(nSmallLD >= correctNsmallLD); // cannot test equality b/c of randomness
		REQUIRE( nSmallLD + nLargeLD < jaccValues.size() );

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
		outAndBim.outputFileName = ldFileName;
		bedGTB.allJaccardLD(outAndBim, nChunks);

		tmpLDfile.open(ldFileName, std::ios::in);
		jaccValues.clear();
		std::getline(tmpLDfile, line); // header
		while ( std::getline(tmpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.push_back( stof(field) );
		}
		tmpLDfile.close();
		std::remove( ldFileName.c_str() ); // NOLINT
		nLargeLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[upperCutOff](float eachValue) {return eachValue >= upperCutOff;}
		);
		REQUIRE(nLargeLD >= correctNlargeLD); // cannot test equality b/c of randomness
		nSmallLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[lowerCutOff](float eachValue) {return eachValue <= lowerCutOff;}
		);
		REQUIRE(nSmallLD >= correctNsmallLD); // cannot test equality b/c of randomness
		REQUIRE( nSmallLD + nLargeLD < jaccValues.size() );
	}
}

TEST_CASE("GenoTableHash methods work", "[gtHash]") {
	const std::string logFileName("../tests/binTest.log");
	const std::string inputBedName("../tests/ind197_397.bed");
	constexpr uint32_t nIndividuals{197};
	constexpr uint16_t kSketches{29};
	constexpr size_t nRowsPerBand{5};
	constexpr float invKlowBound{0.05};
	constexpr float invKhighBound{0.95};
	constexpr uint32_t lowCountMin{5000};
	constexpr uint32_t highCountMin{1000};
	constexpr size_t nThreads{4};
	constexpr uint32_t nLoci{397};
	constexpr size_t totNpairs{static_cast<size_t>(nLoci) * (static_cast<size_t>(nLoci) - 1) / 2};
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
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{smallNind, smallSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{nIndividuals, smallSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be at least three") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{nIndividuals, kGtN}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be smaller than the number of individuals") );
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
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(emptyMACvec, sketchParameters, logFileName),
				Catch::Matchers::StartsWith("ERROR: empty vector of minor allele counts") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(smallMACvec, BayesicSpace::IndividualAndSketchCounts{smallNind, kSketches}, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(smallMACvec, BayesicSpace::IndividualAndSketchCounts{undivNind, kSketches}, logFileName),
				Catch::Matchers::StartsWith("ERROR: length of allele count vector") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(macVector, BayesicSpace::IndividualAndSketchCounts{nIndividuals, kGtN}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be smaller than the number of individuals") );
	}
	SECTION("GenoTableHash .bed file constructor and methods with correct data") {
		BayesicSpace::GenoTableHash bedHSH(inputBedName, sketchParameters, nThreads, logFileName);
		const std::string tmpJacFile("../tests/tmpJac.tsv");
		BayesicSpace::InOutFileNames tmpFileGrp{};
		tmpFileGrp.outputFileName = tmpJacFile;
		tmpFileGrp.inputFileName  = "";
		constexpr size_t forcedChunks{3};
		bedHSH.allHashLD(tmpFileGrp, forcedChunks);
		std::fstream hashLDfile(tmpJacFile, std::ios::in);
		std::vector<float> jaccValues; 
		std::string line;
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.emplace_back( stof(field) );
		}
		hashLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE(jaccValues.size() == totNpairs);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKlowBound](float value) {return value <= invKlowBound;}
			 ) >= lowCountMin
		);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKhighBound](float value) {return value >= invKhighBound;}
			) >= highCountMin
		);
		std::vector<BayesicSpace::HashGroup> groups{bedHSH.makeLDgroups(nRowsPerBand)};
		/*
		auto chunkSizes{BayesicSpace::makeChunkSizes(groups.back().cumulativeNpairs, forcedChunks)};
		BayesicSpace::HashGroupItPairCount startPair{};
		startPair.hgIterator = groups.cbegin();
		startPair.pairCount  = 0;
		std::vector< std::pair<BayesicSpace::HashGroupItPairCount, BayesicSpace::HashGroupItPairCount> > groupRanges;
		groupRanges.reserve( chunkSizes.size() );
		for (const auto &eachCS : chunkSizes) {
			groupRanges.emplace_back( BayesicSpace::makeGroupRanges(groups, startPair, eachCS) );
			startPair = groupRanges.back().second;
		}
		*/
		//const auto groupRanges{BayesicSpace::makeGroupRanges(groups, chunkSizes)};

		std::fstream stdOut("../tests/groupIdxPairs.tsv", std::ios::out | std::ios::trunc);
		stdOut << "locus1\tlocus2\n";
		for (const auto &eachGroup : groups) {
			for (size_t iRow = 1; iRow < eachGroup.locusIndexes.size(); ++iRow) {
				for (size_t jCol = 0; jCol < iRow; ++jCol) {
					stdOut << eachGroup.locusIndexes.at(iRow) + 1 << "\t" << eachGroup.locusIndexes.at(jCol) + 1 << "\n";
				}
			}
		}
		stdOut.close();

		std::string smFileName("../tests/smTest.tsv");
		tmpFileGrp.outputFileName = smFileName;
		tmpFileGrp.inputFileName  = "";
		bedHSH.ldInGroupsSM(nRowsPerBand, tmpFileGrp, forcedChunks);

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
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {return grp1.locusIndexes.at(0) < grp2.locusIndexes.at(0);}
			)
		);
		/*
		std::vector<BayesicSpace::IndexedPairSimilarity> groupLD{bedHSH.ldInGroups(nRowsPerBand)};
		REQUIRE(groupLD.size() <= totNpairs);
		REQUIRE(std::is_sorted(
				groupLD.cbegin(),
				groupLD.cend(),
				[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second) {
					return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
				}
			)
		);
		REQUIRE(std::all_of(
				groupLD.cbegin(),
				groupLD.cend(),
				[&bedHLD](const BayesicSpace::IndexedPairSimilarity &eachGrpPair) {
					auto findIt = std::find_if(
							bedHLD.cbegin(),
							bedHLD.cend(),
							[&eachGrpPair](const BayesicSpace::IndexedPairSimilarity &allPair) {
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
		*/
	}
	SECTION("GenoTableHash mac vector constructor and methods with correct data") {
		/*
		BayesicSpace::GenoTableHash vecHSH(macVector, sketchParameters, nThreads, logFileName);
		const std::string tmpJacFile("../tests/tmpJac.tsv");
		BayesicSpace::InOutFileNames tmpFileGrp{};
		tmpFileGrp.outputFileName = tmpJacFile;
		tmpFileGrp.inputFileName  = "";
		constexpr size_t forcedChunks{3};
		vecHSH.allHashLD(tmpFileGrp, forcedChunks);
		std::fstream hashLDfile(tmpJacFile, std::ios::in);
		std::vector<float> jaccValues; 
		std::string line;
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.emplace_back( stof(field) );
		}
		hashLDfile.close();
		//std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE(jaccValues.size() == totNpairs);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKlowBound](float value) {return value <= invKlowBound;}
			 ) >= lowCountMin
		);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKhighBound](float value) {return value >= invKhighBound;}
			) >= highCountMin
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
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {return grp1.locusIndexes.at(0) < grp2.locusIndexes.at(0);}
			)
		);
		std::vector<BayesicSpace::IndexedPairSimilarity> groupLD{vecHSH.ldInGroups(nRowsPerBand)};
		REQUIRE(groupLD.size() <= totNpairs);
		REQUIRE(std::is_sorted(
				groupLD.cbegin(),
				groupLD.cend(),
				[](const BayesicSpace::IndexedPairSimilarity &first, const BayesicSpace::IndexedPairSimilarity &second) {
					return (first.element1ind == second.element1ind) && (first.element2ind == second.element2ind);
				}
			)
		);
		REQUIRE(std::all_of(
				groupLD.cbegin(),
				groupLD.cend(),
				[&bedHLD](const BayesicSpace::IndexedPairSimilarity &eachGrpPair) {
					auto findIt = std::find_if(
							bedHLD.cbegin(),
							bedHLD.cend(),
							[&eachGrpPair](const BayesicSpace::IndexedPairSimilarity &allPair) {
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
		*/
	}
}
