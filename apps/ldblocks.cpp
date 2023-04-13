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

#include <chrono>

#include "gvarHash.hpp"
#include "vashFunctions.hpp"

int main(int argc, char *argv[]) {

	// set usage message
	std::string cliHelp = "Available command line flags (in any order):\n" 
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
		"  --only-groups      if set (with no value), only group IDs are saved for each locus pair. Ignored if --n-rows-per-band or --hash-size is 0.\n"
		"  --add-locus-names  if set (with no value) adds locus names from the corresponding .bim file to the output (otherwise base-1 indexes are listed).\n"
		"  Invalid (e.g., non-integer) flag values are replaced by defaults, if available\n";
	try {
		std::unordered_map <std::string, std::string> clInfo;
		std::unordered_map <std::string, std::string> stringVariables;
		std::unordered_map <std::string, int> intVariables;
		BayesicSpace::parseCL(argc, argv, clInfo);
		BayesicSpace::extractCLinfo(clInfo, intVariables, stringVariables);

		const size_t kSketches{static_cast<size_t>(intVariables["hash-size"])};
		const size_t nIndiv{static_cast<size_t>(intVariables["n-individuals"])};
		const size_t dotPos = stringVariables["input-bed"].rfind('.');
		std::string bimFileName(stringVariables["input-bed"], 0, dotPos);
		bimFileName += ".bim";
		if (kSketches == 0) {
			BayesicSpace::GenoTableBin allJaccard;
			if (intVariables["threads"] < 1) {
				allJaccard = BayesicSpace::GenoTableBin(stringVariables["input-bed"], nIndiv, stringVariables["log-file"]);
			} else {
				const auto nThreads = static_cast<size_t>(intVariables["threads"]);
				allJaccard = BayesicSpace::GenoTableBin(stringVariables["input-bed"], nIndiv, stringVariables["log-file"], nThreads);
			}
			if (stringVariables["add-locus-names"] == "set") {
				allJaccard.allJaccardLD(bimFileName, stringVariables["out-file"]);
			} else {
				allJaccard.allJaccardLD(stringVariables["out-file"]);
			}
			if (stringVariables["log-file"] != "none") {
				allJaccard.saveLogFile();
			}
		} else {
			BayesicSpace::GenoTableHash groupLD;
			if (intVariables["threads"] < 1) {
				groupLD = BayesicSpace::GenoTableHash(stringVariables["input-bed"], nIndiv, kSketches, stringVariables["log-file"]);
			} else {
				const auto nThreads{static_cast<size_t>(intVariables["threads"])};
				groupLD = BayesicSpace::GenoTableHash(stringVariables["input-bed"], nIndiv, kSketches, nThreads, stringVariables["log-file"]);
			}
			if (intVariables["n-rows-per-band"] == 0) {
				if (stringVariables["add-locus-names"] == "set") {
					groupLD.allHashLD(bimFileName, stringVariables["out-file"]);
				} else {
					groupLD.allHashLD(stringVariables["out-file"]);
				}
			} else {
				const auto rowsPB{static_cast<size_t>(intVariables["n-rows-per-band"])};
				if (stringVariables["add-locus-names"] == "set") {
					if (stringVariables["only-groups"] == "set") {
						groupLD.makeLDgroups(rowsPB, bimFileName, stringVariables["out-file"]);
					} else {
						groupLD.ldInGroups(rowsPB, bimFileName, stringVariables["out-file"]);
					}
				} else {
					if (stringVariables["only-groups"] == "set") {
						//groupLD.makeLDgroups(rowsPB, stringVariables["out-file"]);
						auto time1 = std::chrono::high_resolution_clock::now();
						std::unordered_map< uint32_t, std::vector<size_t> > res1{groupLD.makeLDgroups(rowsPB)};
						auto time2 = std::chrono::high_resolution_clock::now();
						std::chrono::duration<float, std::milli> mapTime = time2 - time1;
						time1 = std::chrono::high_resolution_clock::now();
						std::vector< std::vector<size_t> > res2{groupLD.makeLDgroupsVec(rowsPB)};
						time2 = std::chrono::high_resolution_clock::now();
						std::chrono::duration<float, std::milli> vecTime = time2 - time1;
						std::cout << mapTime.count() << "\t" << vecTime.count() << "\n";
					} else {
						//groupLD.ldInGroups(rowsPB, stringVariables["out-file"]);
						auto time1 = std::chrono::high_resolution_clock::now();
						groupLD.ldInGroups(rowsPB, stringVariables["out-file"]);
						auto time2 = std::chrono::high_resolution_clock::now();
						std::chrono::duration<float, std::milli> mapTime = time2 - time1;
						time1 = std::chrono::high_resolution_clock::now();
						groupLD.ldInGroupsVec(rowsPB, stringVariables["out-file"]);
						time2 = std::chrono::high_resolution_clock::now();
						std::chrono::duration<float, std::milli> vecTime = time2 - time1;
						std::cout << mapTime.count() << "\t" << vecTime.count() << "\n";
					}
				}
			}
			if (stringVariables["log-file"] != "none") {
				groupLD.saveLogFile();
			}
		}
	} catch(std::string &problem) {
		std::cerr << problem << "\n";
		std::cerr << cliHelp;
		return 4;
	}
	return 0;
}
