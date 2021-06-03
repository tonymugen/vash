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
 *
 */

#include <bits/types/clock_t.h>
#include <bits/types/time_t.h>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

#include "../src/gvarHash.hpp"

using std::string;
using std::vector;
using std::fstream;
using std::ios;
using std::cerr;

using namespace BayesicSpace;

int main(int argc, char *argv[]){
	const string inFileName(argv[1]);
	const int Nindv = atoi(argv[2]);
	if (Nindv <= 1){
		cerr << "Number of individuals (is " << Nindv << ") must be greater than 1\n";
		exit(1);
	}
	const int kSketches = atoi(argv[3]);
	if (kSketches <= 1){
		cerr << "Number of sketches (is " << kSketches << ") must be greater than 1\n";
		exit(2);
	}
	const string outFileName(argv[4]);
	vector<double> result;
	time_t loadTime = clock();
	GenoTable abaTest(inFileName, Nindv, kSketches);
	loadTime = clock() - loadTime;
	clock_t runTime = clock();
	abaTest.allSimilarity(result);
	runTime = clock() - runTime;
	//cout << 1000.0*static_cast<double>(vTime)/CLOCKS_PER_SEC << " " << 1000.0*static_cast<double>(sTime)/CLOCKS_PER_SEC << endl;
}
