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
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include "../src/gvarHash.hpp"

//using std::vector;
using std::string;
using std::to_string;
using std::fstream;
using std::stringstream;
using std::ios;
using std::cerr;

using namespace BayesicSpace;

int main(){
	try {
		const string inFile("testBinaryVar2.tsv");
		string inputLine;
		fstream input;
		vector< vector<int8_t> > genoCodes; // file rows in the inner vector
		input.open(inFile.c_str(), ios::in);
		size_t Ngeno = 0;
		while ( getline(input, inputLine) ){
			stringstream lineSS(inputLine);
			string field;
			vector<int8_t> line;
			while (lineSS >> field){
				if (field == "0"){
					line.push_back(0);
				} else if (field == "1"){
					line.push_back(1);
				} else if (field == "2"){
					line.push_back(2);
				} else if (field == "-9"){
					line.push_back(-9);
				} else {
					throw string("ERROR: invalid genotype value ") + field ;
				}
			}
			if (Ngeno == 0){
				Ngeno = line.size();
				genoCodes.emplace_back(line);
			} else if ( Ngeno != line.size() ){
				throw string("ERROR: new number of genotypes (") + to_string( line.size() ) + string(") not equal to the previous value (") + to_string(Ngeno) + string(")");
			} else {
				genoCodes.emplace_back(line);
			}
		}
		input.close();
		// Convert the vector of vectors to a variant-major genotype vector
		vector<int8_t> genoVec;
		for (size_t jGeno = 0; jGeno < Ngeno; jGeno++){
			for (size_t iInd = 0; iInd < genoCodes.size(); iInd++){
				genoVec.push_back(genoCodes[iInd][jGeno]);
			}
		}
		size_t k = 7;
		GenoTable testTab(genoVec, genoCodes.size(), k);
		vector<double> out;
		testTab.allSimilarity(out);
		fstream outStr;
		outStr.open("allSimOut.txt", ios::out | ios::trunc);
		for (size_t iRow = 0; iRow < genoCodes.size(); iRow++) {
			for (size_t jCol = 0; jCol < genoCodes.size() - 1; jCol++){
				outStr << out[genoCodes.size() * iRow + jCol] << "\t";
			}
			outStr << out[genoCodes.size() * iRow + genoCodes.size() - 1] << "\n";
		}
		outStr.close();
		//const string inFile("testBinary.bed");
		//GenoTable testTab( inFile, 200 );
		//testTab.saveGenoBinary("testBinary.bin");
		/*
		for (size_t i = 0; i < 4; i++){
			std::cout << testTab.getSketchIdx(i) << " ";
		}
		std::cout << "\n";
		*/
	} catch(string problem){
		cerr << problem << "\n";
		exit(1);
	}
}

