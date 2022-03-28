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
		/*
		const string inFile("testBinaryVar.tsv");
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
		*/
		const size_t Ngeno = 2000;
		const size_t Nindv = 2000;
		const size_t k = 13;
		//const string bedFile("sim1test.bed");
		//GenoTable testTab(bedFile, Nindv);
		//testTab.makeIndividualOPH(k);
		//
		// This is a basic and incomplete tped parser just for testing binarization
		const string tpedFileName("sim1.tped");
		string inputLine;
		vector<int> genoCodes;
		fstream input;
		input.open(tpedFileName.c_str(), ios::in);
		while ( getline(input, inputLine) ){
			stringstream lineSS(inputLine);
			string field;
			// metadata
			lineSS >> field;
			lineSS >> field;
			lineSS >> field;
			lineSS >> field;

			string firstGenotype;
			while (lineSS >> field){
				string secondField;
				lineSS >> secondField;
				if ( firstGenotype.empty() ){ // first non-missing genotype has not been seen yet
					if ( (field == "0") || (secondField == "0") ){
						genoCodes.push_back(-9);
					} else {
						firstGenotype = field;
						genoCodes.push_back( static_cast<int>(secondField != field));
					}
				} else {
					if ( (field == "0") || (secondField == "0") ){
						genoCodes.push_back(-9);
					} else {
						genoCodes.push_back( static_cast<int>(field != firstGenotype) + static_cast<int>(secondField != firstGenotype) );
					}
				}
			}
		}
		input.close();
		/*
		for (size_t jG = 0; jG < Ngeno; jG++) {
			uint16_t inByte = 0;
			for (size_t iInd = 0; iInd < Nindv; iInd++){
				std::cout << genoCodes[jG * Nindv + iInd];
				if (inByte == 7){
					inByte = 0;
					std::cout << " | ";
				} else {
					inByte++;
					std::cout << " ";
				}
			}
			std::cout << "\b\n";
		}
		*/
		//GenoTableHash testTab(genoCodes, Nindv, k);
		const string bedFile("sim1.bed");
		GenoTableHash testBed(bedFile, Nindv, k);
		/*
		vector<float> outLD;
		testTab.allJaccardLD(outLD);
		for (const auto &o : outLD){
			std::cout << o << " ";
		}
		std::cout << "\b\n";
		*/
	} catch(string problem){
		cerr << problem << "\n";
		exit(1);
	}
}

