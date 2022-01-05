/*
 * Copyright (c) 2022 Anthony J. Greenberg
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

/// Run tests of performance and accuracy for LD among loci within LD groups
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2022 Anthony J. Greenberg
 * \version 0.1
 *
 * Run correctness and timing tests of linkage disequilibrium among loci within groups determined from LD cut-offs using simulated data.
 *
 */

#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "../src/gvarHash.hpp"

int main(int argc, char *argv[]){
	const string inFileName(argv[1]);
	const int Nindv = atoi(argv[2]);
	if (Nindv <= 1){
		std::cerr << "Number of individuals (is " << Nindv << ") must be greater than 1\n";
		exit(1);
	}
	const int kSketches = atoi(argv[3]); // number of OPH sketches to make
	const int nElements = atoi(argv[4]); // number of OPH sketches to consider for similarity
	if (nElements > kSketches){
		std::cerr << "Number of sketches to consider (" << nElements << ") cannot be bigger than number of sketches made (" << kSketches << ")\n";
		exit(2);
	}
	const int hammingCutoff  = atoi(argv[5]);
	const int lookBackNumber = atoi(argv[6]);
	try {
		const std::string outFileName(argv[7]);
		BayesicSpace::GenoTable grpTst(inFileName, Nindv);
		grpTst.makeIndividualOPH(kSketches);
		grpTst.groupByLD(hammingCutoff, kSketches, lookBackNumber, outFileName);
	} catch (std::string problem) {
		std::cerr << problem << "\n";
		exit(3);
	}
}
