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


/// Build local linkage disequilibrium blocks
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2022 Anthony J. Greenberg
 * \version 0.5
 *
 * Uses the variant hashing library to build local LD blocks from _plink_ .bed files.
 *
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <stdexcept>

#include "src/gvarHash.hpp"

void parseCL(int&, char**, std::unordered_map<string, string> &);

/** \brief Command line parser
 *
 * Maps flags to values. Flags assumed to be of the form `--flag-name value`.
 *
 * \param[in] argc size of the `argv` array
 * \param[in] argv command line input array
 * \param[out] cli map of tags to values
 */
void parseCL(int &argc, char **argv, std::unordered_map<string, string> &cli){
	// set to true after encountering a flag token (the characters after the dash)
	bool val = false;
	// store the token value here
	string curFlag;

	for (int iArg = 1; iArg < argc; iArg++) {
		const char *pchar = argv[iArg];

		if ( (pchar[0] == '-') && (pchar[1] == '-') ) { // encountered the double dash, look for the token after it
			if (!pchar[2]) {
				std::cerr << "WARNING: forgot character after dash. Ignoring.\n";
				continue;
			}
			// what follows the dash?
			val     = true;
			curFlag = pchar + 2;

		} else {
			if (val) {
				val = false;
				cli[curFlag] = pchar;
			} else {
				std::cerr << "WARNING: command line value " << pchar << " ignored because it is not preceded by a flag.\n";
			}

		}

	}
}

int main(int argc, char *argv[]){
	string inFileName;
	string logFileName;
	string outFileName;
	float sparsity;
	int inputKsketches{0};
	int inputThreads{-1};
	int Nindv{0};
	int minGroupSize{0};
	int groupingHashSize{0};

	// set usage message
	std::string cliHelp = "Command line flags (in any order):\n  --input-bed file_name (input file name; required)\n  --sparsity sparsity from 0.0 (all by all) to 1.0 (small expected group size, sparse LD matrix)\n  --n-individuals number_of_individuals (must be 3 or more; required)\n  --hash-size hash_size default is 60\n  --min-group-size min_group_size (disregard groups smaller than this; default 5)\n  --grouping-hash-size hash_size_for_grouping (size of hash subset to use for grouping; defaults to full hash size)\n  --threads number_of_threads (maximal number of threads to use; defaults to maximal available)\n  --log-file log_file_name (log file name; default is ldblocks.log; log file not saved if 'none')\n  --out-file output_file_name (output name file; default ldblocksOut.tsv)\n";
	std::unordered_map <string, string> clInfo;
	parseCL(argc, argv, clInfo);

	if ( clInfo.empty() ){
		std::cerr << "Available command line flags:\n";
		std::cerr << cliHelp;
		exit(1);
	}
	auto clIter = clInfo.begin(); // iterator of the command line flags

	clIter = clInfo.find("input-bed");
	if ( clIter == clInfo.end() ){ // if not there, complain
		std::cerr << "ERROR: specification of the input file name is required\n";
		std::cerr << cliHelp;
		exit(2);
	} else {
		inFileName = clIter->second;
	}

	clIter = clInfo.find("sparsity");
	if ( clIter == clInfo.end() ){ // if not there, set to default
		sparsity = 0.0;
	} else {
		sparsity = std::stof(clIter->second);
	}
	if ( (sparsity < 0.0) || (sparsity > 1.0) ){
		std::cerr << "ERROR: sparsity must be between 0.0 and 1.0\n";
		exit(2);
	}

	clIter = clInfo.find("n-individuals");
	if ( clIter == clInfo.end() ){ // if not there, set to default
		std::cerr << "ERROR: specification of the number of individuals is required\n";
		std::cerr << cliHelp;
		exit(2);
	} else {
		try {
			Nindv = stoi(clIter->second);
		} catch(const std::invalid_argument& ia){
			std::cerr << "ERROR: Provided number of individuals is not an integer\n";
			exit(2);
		}
	}
	if (Nindv < 3){
		std::cerr << "ERROR: must have at least 3 individuals\n";
		exit(2);
	}

	clIter = clInfo.find("hash-size");
	if ( clIter == clInfo.end() ){ // if not there, set to default
		inputKsketches = 60;
	} else {
		try {
			inputKsketches = stoi(clIter->second);
		} catch(const std::invalid_argument& ia){
			std::cerr << "ERROR: Provided hash size is not an integer\n";
			exit(2);
		}
	}
	if (inputKsketches <= 3){
		std::cerr << "ERROR: hash length must be 3 or more\n";
		exit(2);
	}

	clIter = clInfo.find("min-group-size");
	if ( clIter == clInfo.end() ){ // if not there, set to default
		minGroupSize = 5;
	} else {
		try {
			minGroupSize = stoi(clIter->second);
		} catch(const std::invalid_argument& ia){
			std::cerr << "ERROR: Provided minimal group size is not an integer\n";
			exit(2);
		}
	}
	if (minGroupSize < 0){
		std::cerr << "ERROR: minimal group size must be non-negative\n";
		exit(2);
	}

	clIter = clInfo.find("grouping-hash-size");
	if ( clIter == clInfo.end() ){ // if not there, set to default
		groupingHashSize = inputKsketches;
	} else {
		try {
			groupingHashSize = stoi(clIter->second);
		} catch(const std::invalid_argument& ia){
			std::cerr << "ERROR: Provided group precision is not an integer\n";
			exit(2);
		}
	}
	if (groupingHashSize < 3){
		std::cerr << "ERROR: must have at least 3 sketches for locus grouping\n";
		exit(2);
	}

	clIter = clInfo.find("threads");
	if ( clIter == clInfo.end() ){ // if not there, set to negative value (will be changed to default later)
		inputThreads = -1;
	} else {
		try {
			inputThreads = stoi(clIter->second);
		} catch(const std::invalid_argument& ia){
			std::cerr << "ERROR: Provided number of threads is not an integer\n";
			exit(2);
		}
	}

	clIter = clInfo.find("log-file");
	if ( clIter == clInfo.end() ){ // if not there, set to default
		logFileName = "ldblocks.log";
	} else {
		logFileName = clIter->second;
	}

	clIter = clInfo.find("out-file");
	if ( clIter == clInfo.end() ){ // if not there, set to default
		outFileName = "ldblocksOut.tsv";
	} else {
		outFileName = clIter->second;
	}
	// Proceed to actual analysis
	try {
		const size_t kSketches = static_cast<size_t>(inputKsketches);
		const size_t nIndiv    = static_cast<size_t>(Nindv);
		BayesicSpace::GenoTableHash groupLD;
		if (inputThreads < 1){
			groupLD = BayesicSpace::GenoTableHash(inFileName, nIndiv, kSketches, logFileName);
		} else {
			const size_t nThreads = static_cast<size_t>(inputThreads);
			groupLD               = BayesicSpace::GenoTableHash(inFileName, nIndiv, kSketches, nThreads, logFileName);
		}
		if (sparsity == 0.0){
			groupLD.allHashLD(outFileName);
			groupLD.saveLogFile();
		} else {
			const size_t kSubs       = static_cast<size_t>(groupingHashSize);
			const size_t smGrpSize   = static_cast<size_t>(minGroupSize);
			groupLD.ldInGroups(kSubs, sparsity, smGrpSize, outFileName);
			if (logFileName != "none"){
				groupLD.saveLogFile();
			}
		}
	} catch(std::string problem){
		std::cerr << problem << "\n";
		exit(4);
	}
}
