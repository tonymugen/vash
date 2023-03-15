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

/// Testing bit-parallel bed to binary conversion

#include <cstring>
#include <vector>
#include <array>
#include <string>
#include <bitset>
#include <iostream>
#include <fstream>
#include <x86intrin.h>

#include "random.hpp"

void binarizeBedLocus(const size_t &bedIdx, const size_t &bedLocusSize, const std::vector<char> &bedLocus, const size_t &binIdx, const size_t &binLocusSize, BayesicSpace::RanDraw &prng, std::vector<uint8_t> &binLocus) {
	constexpr size_t word64size{8};                              // size of uint64_t word in bytes
	constexpr size_t word32size{4};                              // size of uint32_t word in bytes
	constexpr auto wordMask{static_cast<uint64_t>(-word64size)}; // for rounding down to nearest divisible by 8 number 
	constexpr uint64_t secondBitMask{0xaaaaaaaaaaaaaaaa};        // all bed second bits set
	constexpr uint64_t firstBitMask{~secondBitMask};             // all bed first bits set
	const size_t nEvenBedBytes{bedLocusSize & wordMask};         // number of bed bytes that fully fit into 64-bit words
	size_t iBedByte{bedIdx};
	size_t iBinByte{binIdx};
	while (iBedByte < nEvenBedBytes){
		uint64_t bedWord{0};
		memcpy(&bedWord, bedLocus.data() + iBedByte, word64size);
		bedWord = ~bedWord;
		auto binBits{static_cast<uint32_t>( _pext_u64(bedWord, firstBitMask) )};
		const auto secondBedBits{static_cast<uint32_t>( _pext_u64(bedWord, secondBitMask) )};
		const uint32_t maHoms{binBits & secondBedBits};
		uint32_t allHets{(~maHoms) & binBits};     // set at all het positions
		const auto ranBits{static_cast<uint32_t>( prng.ranInt() )};
		allHets &= ranBits;
		binBits  = maHoms | allHets;
		memcpy(binLocus.data() + iBinByte, &binBits, word32size);
		iBedByte += word64size;
		iBinByte += word32size;
	}
	// add the trailing word
}

int main() {
	constexpr size_t byteSize{8};                           // size of byte in bits
	constexpr size_t llWordSize{8};                         // size of uint64_t word in bytes
	constexpr size_t lWordInBits{32};                       // uint32_t size in bits
	constexpr uint64_t secondBitMask{0xAAAAAAAAAAAAAAAA};
	constexpr uint64_t firstBitMask{~secondBitMask};
	constexpr size_t seed{1711};
	const std::string bedFile("bedConv.bed");
	BayesicSpace::RanDraw prng(seed);
	std::array<char, 3> magicBytes{0};
	std::vector<char> bedBytes(llWordSize,0);
	std::fstream bedIn;
	bedIn.open(bedFile, std::ios::in | std::ios::binary);
	bedIn.read( magicBytes.data(), magicBytes.size() );
	bedIn.read(bedBytes.data(), llWordSize);
	bedIn.close();
	uint64_t bedWord{0};
	memcpy(&bedWord, bedBytes.data(), llWordSize);
	bedWord = ~bedWord;
	const auto firstBedBits{static_cast<uint32_t>( _pext_u64(bedWord, firstBitMask) )};
	const auto secondBedBits{static_cast<uint32_t>( _pext_u64(bedWord, secondBitMask) )};
	const uint32_t maHoms{firstBedBits & secondBedBits};
	const uint32_t allHets{(~maHoms) & firstBedBits};     // set at all het positions
	const auto ranBits{static_cast<uint32_t>( prng.ranInt() )};
	std::cout << std::bitset<llWordSize * byteSize>(bedWord) << "\n";
	std::cout << std::bitset<llWordSize * byteSize>(firstBitMask) << "\n";
	std::cout << std::bitset<llWordSize * byteSize>( static_cast<uint64_t>(firstBedBits) ) << "\n\n";
	std::cout << std::bitset<llWordSize * byteSize>(bedWord) << "\n";
	std::cout << std::bitset<llWordSize * byteSize>(secondBitMask) << "\n";
	std::cout << std::bitset<llWordSize * byteSize>( static_cast<uint32_t>(secondBedBits) ) << "\n\n";
	std::cout << std::bitset<lWordInBits>(firstBedBits) << "\n";
	std::cout << std::bitset<lWordInBits>(secondBedBits) << "\n";
	std::cout << std::bitset<lWordInBits>(maHoms) << "\n";
	std::cout << std::bitset<lWordInBits>(allHets) << "\n";
	std::cout << std::bitset<lWordInBits>(ranBits) << "\n";
	std::cout << std::bitset<lWordInBits>(allHets & ranBits) << "\n";
	std::cout << std::bitset<lWordInBits>( maHoms | (allHets & ranBits) ) << "\n";
} 

