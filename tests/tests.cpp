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
#include <vector>
#include <array>

#include <iostream>

#include "gvarHash.hpp"
#include "vashFunctions.hpp"

#include "catch2/catch_test_macros.hpp"

TEST_CASE("Can count set bits correctly", "[countSetBits]") { // NOLINT
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

TEST_CASE("MurMurHash works properly", "[MurMurHash]") { // NOLINT
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
