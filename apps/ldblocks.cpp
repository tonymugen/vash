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
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <cmath>
#include <stdexcept>

#include "gvarHash.hpp"
#include "vashFunctions.hpp"

namespace BayesicSpace {
	/** \brief Full Jaccard estimates 
	 *
	 * Run the full Jaccard estimates after binary conversion.
	 *
	 * \param[in] stringVariables string-valued input flag variables
	 * \param[in] intVariables integer-valued input flag variables
	 * \param[in] bimFileName .bim file name
	 *
	 */
	void fullJaccard(const std::unordered_map<std::string, std::string> &stringVariables, const std::unordered_map<std::string, int> &intVariables, const std::string &bimFileName) {
		const auto nIndiv{static_cast<uint32_t>( intVariables.at("n-individuals") )};
		BayesicSpace::GenoTableBin allJaccard;
		try {
			if (intVariables.at("threads") < 1) {
				allJaccard = BayesicSpace::GenoTableBin( stringVariables.at("input-bed"), nIndiv, stringVariables.at("log-file") );
			} else {
				const auto nThreads = static_cast<size_t>( intVariables.at("threads") );
				allJaccard = BayesicSpace::GenoTableBin(stringVariables.at("input-bed"), nIndiv, stringVariables.at("log-file"), nThreads);
			}
			if (stringVariables.at("add-locus-names") == "set") {
				InOutFileNames bimAndLD{};
				bimAndLD.inputFileName  = bimFileName;
				bimAndLD.outputFileName = stringVariables.at("out-file");
				allJaccard.allJaccardLD(bimAndLD);
			} else {
				InOutFileNames bimAndLD{};
				bimAndLD.outputFileName = stringVariables.at("out-file");
				allJaccard.allJaccardLD(bimAndLD);
			}
			if (stringVariables.at("log-file") != "none") {
				allJaccard.saveLogFile();
			}
		} catch(std::string &problem) {
			if (stringVariables.at("log-file") != "none") {
				allJaccard.saveLogFile();
			}
			throw std::move(problem);
		}
	}
	/** \brief Hash-based Jaccard estimates 
	 *
	 * Run the hash-based Jaccard estimates after binary conversion, possibly in groups.
	 *
	 * \param[in] stringVariables string-valued input flag variables
	 * \param[in] intVariables integer-valued input flag variables
	 * \param[in] floatVariables float-valued input flag variables
	 * \param[in] bimFileName .bim file name
	 *
	 */
	void hashJaccard(const std::unordered_map<std::string, std::string> &stringVariables, const std::unordered_map<std::string, int> &intVariables,
			const std::unordered_map<std::string, float> &floatVariables, const std::string &bimFileName) {
		IndividualAndSketchCounts indivSketches{};

		indivSketches.nIndividuals   = static_cast<uint32_t>( intVariables.at("n-individuals") );
		indivSketches.kSketches      = static_cast<uint16_t>( intVariables.at("hash-size") );
		const float similarityCutOff = floatVariables.at("min-similarity");
		BayesicSpace::GenoTableHash groupLD;
		try {
			if (intVariables.at("threads") < 1) {
				groupLD = BayesicSpace::GenoTableHash( stringVariables.at("input-bed"), indivSketches, stringVariables.at("log-file") );
			} else {
				const auto nThreads{static_cast<size_t>( intVariables.at("threads") )};
				groupLD = BayesicSpace::GenoTableHash( stringVariables.at("input-bed"), indivSketches, nThreads, stringVariables.at("log-file") );
			}
			if (intVariables.at("n-rows-per-band") == 0) {
				if (stringVariables.at("add-locus-names") == "set") {
					InOutFileNames bimAndLD{};
					bimAndLD.inputFileName  = bimFileName;
					bimAndLD.outputFileName = stringVariables.at("out-file");
					groupLD.allHashLD(similarityCutOff, bimAndLD);
				} else {
					BayesicSpace::InOutFileNames outFile{};
					outFile.outputFileName = stringVariables.at("out-file");
					outFile.inputFileName  = "";
					groupLD.allHashLD(similarityCutOff, outFile);
				}
			} else {
				const auto rowsPB{static_cast<size_t>( intVariables.at("n-rows-per-band") )};
				if (stringVariables.at("add-locus-names") == "set") {
					if (stringVariables.at("only-groups") == "set") {
						InOutFileNames bimAndLD{};
						bimAndLD.inputFileName  = bimFileName;
						bimAndLD.outputFileName = stringVariables.at("out-file");
						groupLD.makeLDgroups(rowsPB, bimAndLD);
					} else {
						InOutFileNames bimAndLD{};
						bimAndLD.inputFileName  = bimFileName;
						bimAndLD.outputFileName = stringVariables.at("out-file");
						groupLD.ldInGroups(rowsPB, similarityCutOff, bimAndLD);
					}
				} else {
					if (stringVariables.at("only-groups") == "set") {
						InOutFileNames bimAndLD{};
						bimAndLD.inputFileName  = "";
						bimAndLD.outputFileName = stringVariables.at("out-file");
						groupLD.makeLDgroups(rowsPB, bimAndLD);
					} else {
						BayesicSpace::InOutFileNames outFile{};
						outFile.outputFileName = stringVariables.at("out-file");
						outFile.inputFileName  = "";
						groupLD.ldInGroups(rowsPB, similarityCutOff, outFile);
					}
				}
			}
			if (stringVariables.at("log-file") != "none") {
				groupLD.saveLogFile();
			}
		} catch(std::string &problem) {
			if (stringVariables.at("log-file") != "none") {
				groupLD.saveLogFile();
			}
			throw std::move(problem);
		}
	}
}


int main(int argc, char *argv[]) {

	// set usage message
	const std::string cliHelp = "Available command line flags (in any order):\n" 
		"  --input-bed        file_name (input file name; required).\n"
		"  --n-individuals    number_of_individuals (must be 3 or more; required).\n"
		"  --n-rows-per-band  number_of_rows (number of rows per band in the hashed genotype matrix).\n"
		"                     Set to 0 or omit to get all by all estimates. The value is ignored if hash size is set to 0 (for full Jaccard similarity estimates).\n"
		"                     This parameter must be smaller than hash size and controls the similarity cut-off; larger values lead to sparser similarity matrices.\n"
		"  --hash-size        hash_size; must be smaller than the number of individuals.\n"
		"                     Larger values give better similarity estimates at the expense of speed. Set to 0 or omit to obtain precise Jaccard similarity estimates.\n"
		"  --threads          number_of_threads (maximal number of threads to use; defaults to maximal available).\n"
		"  --min-similarity   minimal similarity value for pairs to be saved.\n"
		"  --log-file         log_file_name (log file name; default is ldblocks.log; log file not saved if 'none').\n"
		"  --out-file         output_file_name (output name file; default ldblocksOut.tsv).\n"
		"  --only-groups      if set (with no value), only group IDs are saved for each locus pair. Ignored if --n-rows-per-band or --hash-size is 0.\n"
		"  --add-locus-names  if set (with no value) adds locus names from the corresponding .bim file to the output (otherwise base-1 indexes are listed).\n"
		"  Invalid (e.g., non-integer) flag values are replaced by defaults, if available\n";

	try {
		std::unordered_map <std::string, std::string> clInfo;
		std::unordered_map <std::string, std::string> stringVariables;
		std::unordered_map <std::string, int> intVariables;
		std::unordered_map <std::string, float> floatVariables;
		BayesicSpace::parseCL(argc, argv, clInfo);
		BayesicSpace::extractCLinfo(clInfo, intVariables, floatVariables, stringVariables);

		const size_t dotPos = stringVariables.at("input-bed").rfind('.');
		std::string bimFileName(stringVariables.at("input-bed"), 0, dotPos);
		bimFileName += ".bim";
		if (intVariables.at("hash-size") == 0) {
			BayesicSpace::fullJaccard(stringVariables, intVariables, bimFileName);
		} else {
			BayesicSpace::hashJaccard(stringVariables, intVariables, floatVariables, bimFileName);
		}
	} catch(std::string &problem) {
		std::cerr << problem << "\n";
		std::cerr << cliHelp;
		return 1;
	}
	return 0;
}
