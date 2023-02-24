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
#include <array>
#include <unordered_map>
#include <cmath>
#include <stdexcept>

#include "gvarHash.hpp"

/** \brief Command line parser
 *
 * Maps flags to values. Flags assumed to be of the form `--flag-name value`.
 *
 * \param[in] argc size of the `argv` array
 * \param[in] argv command line input array
 * \param[out] cli map of tags to values
 */
void parseCL(int &argc, char **argv, std::unordered_map<std::string, std::string> &cli){
	// set to true after encountering a flag token (the characters after the dash)
	bool val = false;
	// store the token value here
	std::string curFlag;

	for (int iArg = 1; iArg < argc; iArg++) {
		const char *pchar = argv[iArg];
		if ( (pchar[0] == '-') && (pchar[1] == '-') ) { // encountered the double dash, look for the token after it
			if (pchar[2] == '\0') {
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

/** \brief Extract parameters from parsed command line interface flags
 *
 * Extracts needed variable values, indexed by `std::string` encoded variable names.
 *
 * \param[in] parsedCLI flag values parsed from the command line
 * \param[out] intVariables indexed `int` variables for use by `main()`
 * \param[out] stringVariables indexed `std::string` variables for use by `main()`
 */
void extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI, std::unordered_map<std::string, int> &intVariables, std::unordered_map<std::string, std::string> &stringVariables) {
	intVariables.clear();
	stringVariables.clear();
	const std::array<std::string, 1> requiredStringVariables{"input-bed"};
	const std::array<std::string, 2> optionalStringVariables{"log-file", "out-file"};
	const std::array<std::string, 1> requiredIntVariables{"n-individuals"};
	const std::array<std::string, 3> optionalIntVariables{"hash-size", "threads", "n-rows-per-band"};

	const std::unordered_map<std::string, int>         defaultIntValues{ {"hash-size", 0}, {"threads", -1}, {"n-rows-per-band", 0} };
	const std::unordered_map<std::string, std::string> defaultStringValues{ {"log-file", "ldblocks.log"}, {"out-file", "ldblocksOut.tsv"} };

	if ( parsedCLI.empty() ){
		throw std::string("Available command line flags");
	}
	for (const auto &eachFlag : requiredIntVariables){
		try {
			intVariables[eachFlag] = stoi( parsedCLI.at(eachFlag));
		} catch(const std::exception &problem) {
			throw std::string("ERROR: " + eachFlag + " specification is required and must be an integer");
		}
	}
	for (const auto &eachFlag : optionalIntVariables){
		try {
			intVariables[eachFlag] = stoi( parsedCLI.at(eachFlag));
		} catch(const std::exception &problem) {
			intVariables[eachFlag] = defaultIntValues.at(eachFlag);
		}
	}
	for (const auto &eachFlag : requiredStringVariables){
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			throw std::string("ERROR: " + eachFlag + " specification is required");
		}
	}
	for (const auto &eachFlag : optionalStringVariables){
		try {
			stringVariables[eachFlag] = parsedCLI.at(eachFlag);
		} catch(const std::exception &problem) {
			stringVariables[eachFlag] = defaultStringValues.at(eachFlag);
		}
	}
}

int main(int argc, char *argv[]){

	// set usage message
	std::string cliHelp = "Command line flags (in any order):\n" 
		"  --input-bed        file_name (input file name; required).\n"
		"  --n-individuals    number_of_individuals (must be 3 or more; required).\n"
		"  --n-rows-per-band  number_of_rows (number of rows per band in the hashed genotype matrix).\n"
		"                     Set to 0 or omit to get all by all estimates. The value is ignored if hash size is set to 0 (for full Jaccard similarity estimates).\n"
		"                     This parameter must be smaller than hash size and controls the similarity cut-off; larger values lead to sparser similarity matrices.\n"
		"  --hash-size        hash_size; must be smaller than the number of individuals.\n"
		"                     Larger values give better similarity estimates at the expense of speed. Set to 0 or omit to obtain precise Jaccard similarity estimates.\n"
		"  --threads          number_of_threads (maximal number of threads to use; defaults to maximal available).\n"
		"  --log-file         log_file_name (log file name; default is ldblocks.log; log file not saved if 'none').\n"
		"  --out-file         output_file_name (output name file; default ldblocksOut.tsv).\n"
		"  Invalid (e.g., non-integer) flag values are replaced by defaults, if available\n";
	std::unordered_map <std::string, std::string> clInfo;
	std::unordered_map <std::string, std::string> stringVariables;
	std::unordered_map <std::string, int> intVariables;
	parseCL(argc, argv, clInfo);

	try {
		extractCLinfo(clInfo, intVariables, stringVariables);
	} catch(std::string &problem){
		std::cerr << problem << "\n";
		std::cerr << cliHelp << "\n";
		return 1;
	}

	// Proceed to actual analysis
	try {
		const size_t kSketches{static_cast<size_t>(intVariables["hash-size"])};
		const size_t nIndiv{static_cast<size_t>(intVariables["n-individuals"])};
		if (kSketches == 0){
			if (intVariables["threads"] < 1){
				BayesicSpace::GenoTableBin allJaccard(stringVariables["input-bed"], nIndiv, stringVariables["log-file"]);
				allJaccard.allJaccardLD(stringVariables["out-file"]);
				if (stringVariables["log-file"] != "none"){
					allJaccard.saveLogFile();
				}
			} else {
				const auto nThreads = static_cast<size_t>(intVariables["threads"]);
				BayesicSpace::GenoTableBin allJaccard(stringVariables["input-bed"], nIndiv, stringVariables["log-file"], nThreads);
				allJaccard.allJaccardLD(stringVariables["out-file"]);
				if (stringVariables["log-file"] != "none"){
					allJaccard.saveLogFile();
				}
			}
		} else {
			BayesicSpace::GenoTableHash groupLD;
			if (intVariables["threads"] < 1){
				groupLD = BayesicSpace::GenoTableHash(stringVariables["input-bed"], nIndiv, kSketches, stringVariables["log-file"]);
			} else {
				const auto nThreads{static_cast<size_t>(intVariables["threads"])};
				groupLD = BayesicSpace::GenoTableHash(stringVariables["input-bed"], nIndiv, kSketches, nThreads, stringVariables["log-file"]);
			}
			if (intVariables["n-rows-per-band"] == 0){
				groupLD.allHashLD(stringVariables["out-file"]);
				if (stringVariables["log-file"] != "none"){
					groupLD.saveLogFile();
				}
			} else {
				const auto rowsPB{static_cast<size_t>(intVariables["n-rows-per-band"])};
				groupLD.ldInGroups(rowsPB, stringVariables["out-file"]);
				if (stringVariables["log-file"] != "none"){
					groupLD.saveLogFile();
				}
			}
		}
	} catch(std::string &problem) {
		std::cerr << problem << "\n";
		return 4;
	}
	return 0;
}
