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

/// Run tests of performance for LD among all loci
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2021 Anthony J. Greenberg
 * \version 0.1
 *
 * Run correctness and timing tests of linkage disequilibrium among all loci using simulated data.
 * Output the results in a space-delimited file, where the first two values are execution times and the rest are the upper triangle of the among-locus similarity matrix.
 *
 */

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <bitset>

#include "../src/gvarHash.hpp"

using std::string;
using std::vector;
using std::array;
using std::cout;
using std::bitset;


using namespace BayesicSpace;

uint16_t simHash(const vector<uint16_t> &in) {
	uint16_t hash = 0;
	array<int16_t, 16> v{};
	for (const auto &i : in){
		for (uint16_t j = 0; j < 16; j++){
			v[j] += -1 + 2 * static_cast<int16_t>( 1 & (i >> j) );
		}
	}
	for (const auto &vv : v){
		cout << vv << " ";
	}
	cout << "\b\n";
	for (uint16_t i = 0; i < 16; i++){
		hash |= (static_cast<uint16_t>(v[i] > 0) << i);
	}
	return hash;
}
int main(){
	const vector<uint16_t> testBin = {0b1100110000110100,
									  0b1100110010110010};
	cout << bitset<16>(testBin[0]) << "\n";
	cout << bitset<16>(testBin[1]) << "\n";
	const uint16_t res = simHash(testBin);
	cout << bitset<16>(res) << "\n";
}
