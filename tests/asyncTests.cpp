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

#include <cstddef>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "../src/gvarHash.hpp"

using std::string;
using std::fstream;
using std::stringstream;
using std::ios;
using std::cerr;


using namespace BayesicSpace;

int main(){
	try {
		const string inFileName("sim1_2000_7500_03.bed");
		const size_t nIndividuals = 2000;
		const string outFileNameJacc("testJaccOld.txt");
		vector<float> result;
		GenoTable asTest(inFileName, nIndividuals);
		asTest.allJaccardLD(result);
		fstream outJacc;
		outJacc.open(outFileNameJacc.c_str(), ios::out | ios::trunc);
		for (const auto &r : result){
			outJacc << " " << r;
		}
		outJacc << "\n";
		outJacc.close();
		const string outFileNamePairs("testJaccPairs.txt");
		result.clear();
		asTest.allJaccardLDasyncPairs(result);
		fstream outPairs;
		outPairs.open(outFileNamePairs.c_str(), ios::out | ios::trunc);
		for (const auto &r : result){
			outPairs << " " << r;
		}
		outPairs << "\n";
		outPairs.close();
		const string outFileNameBlocks("testJaccBlocks.txt");
		result.clear();
		asTest.allJaccardLDasyncPairs(result);
		fstream outBlocks;
		outBlocks.open(outFileNameBlocks.c_str(), ios::out | ios::trunc);
		for (const auto &r : result){
			outBlocks << " " << r;
		}
		outBlocks << "\n";
		outBlocks.close();
	} catch(string problem){
		cerr << problem << "\n";
		exit(1);
	}
}
