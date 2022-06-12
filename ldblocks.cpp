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

#include "src/gvarHash.hpp"

void parseCL(int&, char**, std::unordered_map<char, std::string> &);

/** \brief Command line parser
 *
 * Maps flags to values. Flags assumed to be of the form -x.
 *
 * \param[in] argc size of the `argv` array
 * \param[in] argv command line input array
 * \param[out] cli map of tags to values
 */
void parseCL(int &argc, char **argv, std::unordered_map<char, std::string> &cli){
	// set to true after encountering a flag token (the character after the dash)
	bool val = false;
	// store the token value here
	char curFlag;

	for (int iArg = 1; iArg < argc; iArg++) {
		const char *pchar = argv[iArg];

		if (pchar[0] == '-') { // encountered the dash, look for the token after it
			if (!pchar[1]) {
				std::cerr << "WARNING: forgot character after dash. Ignoring.\n";
				continue;
			}
			// what follows the dash?
			val     = true;
			curFlag = pchar[1];

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
	std::string inFileName;
	std::string sparcityType; // can be hard, soft, or none
	int inputKhashes{0};
	int locality{0};

	// set usage message
	std::string cliHelp = "Command line flags (in any order):\n  -f fileName (input file name, with or without extension; required)\n  -s sparcity_type (can be hard, soft, or none; default is hard)\n  -k number of hashes (default is 60)\n  -l locality degree (how many blocks back to examine; default is 5)\n";
	std::unordered_map <char, std::string> clInfo;
	parseCL(argc, argv, clInfo);
	if ( clInfo.empty() ){
		std::cerr << "Available command line flags:\n";
		std::cerr << cliHelp;
		exit(1);
	}
	auto clIter = clInfo.begin(); // iterator of the command line flags

	clIter = clInfo.find('f');
	if (clIter == clInfo.end()) { // if not there, complain
		std::cerr << "ERROR: specification of the input file name is required\n";
		std::cerr << cliHelp;
		exit(2);
	} else {
		inFileName = clIter->second;
	}

	clIter = clInfo.find('s');
	if (clIter == clInfo.end()) { // if not there, set to default
		sparcityType = "hard";
	} else {
		sparcityType = clIter->second;
	}

	clIter = clInfo.find('k');
	if (clIter == clInfo.end()) { // if not there, set to default
		inputKhashes = 60;
	} else {
		inputKhashes = stoi(clIter->second);
	}

	clIter = clInfo.find('l');
	if (clIter == clInfo.end()) { // if not there, set to default
		locality = 5;
	} else {
		locality = stoi(clIter->second);
	}

}
