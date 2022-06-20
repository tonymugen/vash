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

#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <cmath>

#include "../src/gvarHash.hpp"

using std::string;
using std::fstream;
using std::ios;
using std::cerr;
using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::milli;

using namespace BayesicSpace;

int main(int argc, char *argv[]){
	const string inFileName(argv[1]);
	const int Nindv = atoi(argv[2]);
	if (Nindv <= 1){
		cerr << "Number of individuals (is " << Nindv << ") must be greater than 1\n";
		exit(1);
	}
	const size_t kSketches = static_cast<size_t>( atoi(argv[3]) );
	const size_t nThreads  = static_cast<size_t>( atoi(argv[4]) );
	const string outFileName(argv[5]);
	const string outTimeFileName(argv[6]);
	try {
		if (kSketches <= 1){
			auto time1 = high_resolution_clock::now();
			GenoTableBin binTest(inFileName, Nindv, nThreads);
			auto time2 = high_resolution_clock::now();
			duration<float, milli> loadTime = time2 - time1;
			duration<float, milli> runTime;
			time1 = high_resolution_clock::now();
			binTest.allJaccardLD(outFileName);
			time2 = high_resolution_clock::now();
			runTime = time2 - time1;

			fstream output;
			output.open(outTimeFileName.c_str(), ios::out | ios::trunc);
			// output times in milliseconds
			output << loadTime.count() << " " << runTime.count() << "\n";
			output.close();
		} else {
			duration<float, milli> hashTime;
			duration<float, milli> runTime;
			auto time1 = high_resolution_clock::now();
			const string logFileName("lgroups.log");
			GenoTableHash hashTest(inFileName, Nindv, kSketches, nThreads, logFileName);
			auto time2 = high_resolution_clock::now();
			hashTime = time2 - time1;
			time1 = high_resolution_clock::now();
			hashTest.allHashLD(outFileName);
			time2 = high_resolution_clock::now();
			runTime = time2 - time1;
			fstream output;
			output.open(outTimeFileName.c_str(), ios::out | ios::trunc);
			// output times in milliseconds
			output << hashTime.count() << " " << runTime.count() << "\n";
			output.close();
		}
	} catch(string problem){
		cerr << problem << "\n";
		exit(3);
	}
}
