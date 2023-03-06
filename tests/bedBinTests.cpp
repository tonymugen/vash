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


int main() {
	constexpr size_t byteSize{8};             // size of byte in bitse
	constexpr size_t llWordSize{8};           // size of uint64_t word in bytes
	const std::string bedFile("bedConv.bed");
	std::array<char, 3> magicBytes{0};
	std::vector<char> bedBytes(llWordSize,0);
	std::fstream bedIn;
	bedIn.open(bedFile, std::ios::in | std::ios::binary);
	bedIn.read( magicBytes.data(), magicBytes.size() );
	bedIn.read(bedBytes.data(), llWordSize);
	bedIn.close();
	std::cout << std::bitset<byteSize>( static_cast<uint64_t>(bedBytes[0]) ) << "\n";
	uint64_t bedWord{0};
	memcpy(&bedWord, bedBytes.data(), llWordSize);
	std::cout << std::bitset<llWordSize * byteSize>(bedWord) << "\n";
} 

