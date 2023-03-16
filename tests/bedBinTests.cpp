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

void binarizeBedLocus(const size_t &bedIdx, const size_t &bedLocusSize, const std::vector<char> &bedLocus, const size_t &binIdx, BayesicSpace::RanDraw &prng, std::vector<uint8_t> &binLocus) {
	constexpr size_t word64size{8};                                                           // size of uint64_t word in bytes
	constexpr size_t word32size{4};                                                           // size of uint32_t word in bytes
	constexpr auto wordMask{static_cast<uint64_t>(-word64size)};                              // for rounding down to nearest divisible by 8 number 
	constexpr uint64_t secondBitMask{0xaaaaaaaaaaaaaaaa};                                     // all bed second bits set
	constexpr uint64_t firstBitMask{~secondBitMask};                                          // all bed first bits set
	const size_t nEvenBedBytes{bedLocusSize & wordMask};                                      // number of bed bytes that fully fit into 64-bit words
	size_t iBedByte{bedIdx};
	size_t iBinByte{binIdx};
	while (iBedByte < nEvenBedBytes + bedIdx){
		uint64_t bedWord{0};
		memcpy(&bedWord, bedLocus.data() + iBedByte, word64size);
		auto binBits{static_cast<uint32_t>( _pext_u64(bedWord, firstBitMask) )};
		const auto secondBedBits{static_cast<uint32_t>( _pext_u64(bedWord, secondBitMask) )};
		const uint32_t maHoms{binBits & secondBedBits};
		uint32_t allHets{(~maHoms) & secondBedBits};                                          // set at all het positions
		const auto ranBits{static_cast<uint32_t>( prng.ranInt() )};
		allHets &= ranBits;                                                                   // flip some of the het bits randomly
		binBits  = maHoms | allHets;                                                          // add in the set het bits
		memcpy(binLocus.data() + iBinByte, &binBits, word32size);
		iBedByte += word64size;
		iBinByte += word32size;
	}
	if (bedLocusSize > nEvenBedBytes){
		const size_t nTrailingBedBytes{bedLocusSize - nEvenBedBytes};
		const size_t nTrailingBinBytes{nTrailingBedBytes / 2 + (nTrailingBedBytes & 1)};
		uint64_t bedWord{0};
		memcpy(&bedWord, bedLocus.data() + iBedByte, nTrailingBedBytes);
		auto binBits{static_cast<uint32_t>( _pext_u64(bedWord, firstBitMask) )};
		const auto secondBedBits{static_cast<uint32_t>( _pext_u64(bedWord, secondBitMask) )};
		const uint32_t maHoms{binBits & secondBedBits};
		uint32_t allHets{(~maHoms) & secondBedBits};
		const auto ranBits{static_cast<uint32_t>( prng.ranInt() )};
		allHets &= ranBits;
		binBits  = maHoms | allHets;
		memcpy(binLocus.data() + iBinByte, &binBits, nTrailingBinBytes);
	}
}

int main() {
	constexpr size_t byteSize{8};
	constexpr size_t nIndividuals{151};
	constexpr size_t nLoci{6};
	constexpr size_t bedLocusSize{nIndividuals / 4 + static_cast<size_t>( (nIndividuals % 4) > 0 )};
	constexpr size_t binLocusSize{bedLocusSize / 2 + (bedLocusSize & 1)};
	std::cout << bedLocusSize << " " << binLocusSize << "\n";
	constexpr size_t seed{1711};
	const std::string bedFile("bedConv.bed");
	BayesicSpace::RanDraw prng(seed);
	std::array<char, 3> magicBytes{0};
	std::vector<char> bedBytes(bedLocusSize * nLoci, 0);
	std::fstream bedIn;
	bedIn.open(bedFile, std::ios::in | std::ios::binary);
	bedIn.read( magicBytes.data(), magicBytes.size() );
	bedIn.read( bedBytes.data(), static_cast<std::streamsize>( bedBytes.size() ) );
	bedIn.close();
	for (size_t iLocus = 0; iLocus < nLoci; ++iLocus){
		for (size_t iByte = 0; iByte < bedLocusSize; ++iByte){
			std::cout << std::bitset<byteSize>( static_cast<uint8_t>(bedBytes[iLocus * bedLocusSize + iByte]) ) << " ";
		}
		std::cout << "\n";
	}
	for (size_t iLocus = 0; iLocus < nLoci; ++iLocus){
		std::vector<uint8_t> binLocus(binLocusSize, 0);
		binarizeBedLocus(iLocus * bedLocusSize, bedLocusSize, bedBytes, 0, prng, binLocus);
		for (const auto binByte : binLocus){
			std::cout << std::bitset<byteSize>( static_cast<uint8_t>(binByte) ) << " ";
		}
		std::cout << "\n";
	}
} 

