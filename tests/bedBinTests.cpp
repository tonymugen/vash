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
#include <limits>
#include <x86intrin.h>

#include "random.hpp"

void binarizeBedLocus(const size_t &bedIdx, const size_t &bedLocusSize, const std::vector<char> &bedLocus, const size_t &nIndividuals,
							BayesicSpace::RanDraw &prng, const size_t &binIdx, const size_t &binLocusSize, std::vector<uint8_t> &binLocus) {
	constexpr size_t word64size{8};                                                            // size of uint64_t word in bytes
	constexpr size_t word32size{4};                                                            // size of uint32_t word in bytes
	constexpr size_t word32sizeInBits{32};                                                     // size of uint32_t word in bits
	constexpr auto word64mask{static_cast<uint64_t>(-word64size)};                             // for rounding down to nearest divisible by 8 number 
	constexpr auto word32mask{static_cast<uint64_t>(-word32size)};                             // for rounding down to nearest divisible by 4 number 
	constexpr uint64_t secondBitMask{0xaaaaaaaaaaaaaaaa};                                      // all bed second bits set
	constexpr uint64_t firstBitMask{~secondBitMask};                                           // all bed first bits set
	const size_t nEvenBedBytes{bedLocusSize & word64mask};                                     // number of bed bytes that fully fit into 64-bit words
	std::vector<uint32_t> missWords;                                                           // 32-bit words with missing genotype masks
	std::vector<uint32_t> binWords;                                                            // 32-bit words with binarized data
	uint32_t setCount{0};
	uint32_t missingCount{0};
	size_t iBedByte{bedIdx};
	while (iBedByte < nEvenBedBytes + bedIdx){
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
	if (bedLocusSize > nEvenBedBytes){
		const size_t nTrailingBedBytes{bedLocusSize - nEvenBedBytes};
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
	// TODO: add assert() for nIndividuals > missingCount
	if ( (2UL * setCount) > (nIndividuals - missingCount) ){
		size_t missWordIdx{0};
		for (auto &eachBinWord : binWords){
			eachBinWord = (~eachBinWord) & (~missWords[missWordIdx]);                                            // flip the bits in the binary vector and reset the missing bits to 0
			++missWordIdx;
		}
		const auto padShift{static_cast<uint32_t>( (binWords.size() * word32sizeInBits) % nIndividuals )};
		const uint32_t lastWordMask{std::numeric_limits<uint32_t>::max() >> padShift};                           // clear the padding bits after the flip
		binWords.back() = binWords.back() & lastWordMask;
	}
	// copy over the binary bits to the locus vector
	const size_t nEvenBinBytes{binLocusSize & word32mask};                                                       // number of bin bytes that fully fit into 32-bit words
	size_t iBinByte{binIdx};
	size_t iBinWord{0};
	while (iBinByte < nEvenBinBytes){
		memcpy(binLocus.data() + iBinByte, &binWords[iBinWord], word32size);
		iBinByte += word32size;
		++iBinWord;
	}
	if (binLocusSize > nEvenBinBytes){
		const size_t nTrailingBinBytes{binLocusSize - nEvenBinBytes};
		memcpy(binLocus.data() + iBinByte, &binWords[iBinWord], nTrailingBinBytes);
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
		binarizeBedLocus(iLocus * bedLocusSize, bedLocusSize, bedBytes, nIndividuals, prng, 0, binLocusSize, binLocus);
		for (const auto binByte : binLocus){
			std::cout << std::bitset<byteSize>( static_cast<uint8_t>(binByte) ) << " ";
		}
		std::cout << "\n";
	}
} 

