/*
 * Copyright (c) 2026 Anthony J. Greenberg
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

/// Benchmark harness for hash-based LD in groups
/** \file
 * \author Anthony J. Greenberg
 * \version 0.1
 *
 * Standalone timing harness for `GenoTableHash::ldInGroups` on real `.bed` data. It exists to
 * locate the reported single-thread stretch: it times construction, the serial `makeLDgroups`
 * phase in isolation, and the full `ldInGroups` run, and prints the LD-group structure (sizes and
 * pair skew). Build with `VASH_BENCHMARK` defined (the `benchmarks` CMake target does this) so the
 * in-library `[vash-bench]` phase markers additionally break `ldInGroups` down into its parallel
 * block-compute, serial consolidation, sink-merge, and finalize phases, and report the per-block
 * wall-time spread that reveals load imbalance.
 *
 * Usage:
 *   ldGroupBench --input-bed FILE --n-individuals N --hash-size K --n-rows-per-band R
 *                [--threads T] [--min-similarity S] [--max-mem GB]
 *                [--out-file OUT] [--log-file LOG] [--group-reps M] [--skip-ld]
 */

#include <cstdlib>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <thread>
#include <iostream>
#include <iomanip>

#include "gvarHash.hpp"

namespace {
	using clockType = std::chrono::steady_clock;

	// Milliseconds between two steady-clock stamps.
	double millisBetween(const clockType::time_point &begin, const clockType::time_point &end) {
		return std::chrono::duration<double, std::milli>(end - begin).count();
	}

	// Fetch a flag value from the parsed map, or a default when absent.
	// NOLINTNEXTLINE(bugprone-easily-swappable-parameters) key and fallback are distinct roles; callers pass literals
	std::string flagOr(const std::unordered_map<std::string, std::string> &flags, const std::string &key, const std::string &fallback) {
		const auto flagIt = flags.find(key);
		return flagIt == flags.cend() ? fallback : flagIt->second;
	}

	// Parse `--flag value` pairs; a flag with no following value (or followed by another flag) is a bare switch set to "set".
	std::unordered_map<std::string, std::string> parseArgs(const std::vector<std::string> &args) {
		std::unordered_map<std::string, std::string> flags;
		for (size_t iArg = 0; iArg < args.size(); ++iArg) {
			if (args[iArg].rfind("--", 0) != 0) {
				continue;
			}
			const std::string key(args[iArg], 2);
			if ( (iArg + 1 < args.size()) && (args[iArg + 1].rfind("--", 0) != 0) ) {
				flags[key] = args[iArg + 1];
				++iArg;
			} else {
				flags[key] = "set";
			}
		}
		return flags;
	}

	// Summarize and print the LD-group structure that feeds the parallel estimation.
	void reportGroups(const std::vector<BayesicSpace::HashGroup> &groups) {
		if ( groups.empty() ) {
			std::cout << "  (no multi-locus LD groups)\n";
			return;
		}
		std::vector<size_t> sizes;
		sizes.reserve( groups.size() );
		std::vector< std::pair<size_t, uint64_t> > sizeAndPairs;   // (loci, pairs) per group
		sizeAndPairs.reserve( groups.size() );
		for (const auto &eachGroup : groups) {
			const size_t nLoci{eachGroup.locusIndexes.size()};
			const uint64_t nPairs{static_cast<uint64_t>(nLoci) * (nLoci - 1) / 2};
			sizes.push_back(nLoci);
			sizeAndPairs.emplace_back(nLoci, nPairs);
		}
		const uint64_t totalPairs{groups.back().cumulativeNpairs};
		std::sort( sizes.begin(), sizes.end() );
		const size_t minLoci{sizes.front()};
		const size_t maxLoci{sizes.back()};
		const size_t medianLoci{sizes[sizes.size() / 2]};
		const double meanLoci{
			static_cast<double>( std::accumulate(sizes.cbegin(), sizes.cend(), static_cast<size_t>(0)) ) / static_cast<double>( sizes.size() )
		};
		std::sort(
			sizeAndPairs.begin(), sizeAndPairs.end(),
			[](const std::pair<size_t, uint64_t> &lhs, const std::pair<size_t, uint64_t> &rhs) { return lhs.second > rhs.second; }
		);
		std::cout << "  groups (>= 2 loci):       " << groups.size() << "\n";
		std::cout << "  total pairs to estimate:  " << totalPairs << "\n";
		std::cout << "  group size (loci):        min " << minLoci << ", median " << medianLoci
			<< ", mean " << std::fixed << std::setprecision(1) << meanLoci << ", max " << maxLoci << "\n";
		const uint64_t largestPairs{sizeAndPairs.front().second};
		const double largestShare{totalPairs == 0 ? 0.0 : (100.0 * static_cast<double>(largestPairs) / static_cast<double>(totalPairs))};
		std::cout << "  largest group:            " << sizeAndPairs.front().first << " loci -> "
			<< largestPairs << " pairs (" << std::setprecision(2) << largestShare << "% of all pairs)\n";
		const size_t topN{std::min<size_t>(5, sizeAndPairs.size())};
		std::cout << "  top " << topN << " groups by pairs:   ";
		for (size_t iTop = 0; iTop < topN; ++iTop) {
			std::cout << sizeAndPairs[iTop].first << "L/" << sizeAndPairs[iTop].second << "P"
				<< (iTop + 1 < topN ? ", " : "\n");
		}
	}
} // anonymous namespace

int main(int argc, char *argv[]) {
	const std::vector<std::string> args(argv + 1, argv + argc);
	const std::unordered_map<std::string, std::string> flags{parseArgs(args)};
	if ( flags.count("input-bed") == 0 || flags.count("n-individuals") == 0
			|| flags.count("hash-size") == 0 || flags.count("n-rows-per-band") == 0 ) {
		std::cerr <<
			"Usage: ldGroupBench --input-bed FILE --n-individuals N --hash-size K --n-rows-per-band R\n"
			"                    [--threads T] [--min-similarity S] [--max-mem GB]\n"
			"                    [--out-file OUT] [--log-file LOG] [--group-reps M] [--skip-ld]\n";
		return 1;
	}
	try {
		const std::string inputBed{flags.at("input-bed")};
		const auto nIndividuals{static_cast<uint32_t>( std::stoul( flags.at("n-individuals") ) )};
		const auto kSketches{static_cast<uint16_t>( std::stoul( flags.at("hash-size") ) )};
		const size_t nRowsPerBand{ std::stoul( flags.at("n-rows-per-band") ) };
		const size_t reqThreads{ std::stoul( flagOr(flags, "threads", "0") ) };
		const size_t nThreads{reqThreads < 1 ? static_cast<size_t>( std::thread::hardware_concurrency() ) : reqThreads};
		const float minSimilarity{std::stof( flagOr(flags, "min-similarity", "0.0") )};
		const double maxMemGB{std::stod( flagOr(flags, "max-mem", "0.0") )};
		const std::string outFile{flagOr(flags, "out-file", "ldGroupBenchOut.tsv")};
		const std::string logFileArg{flagOr(flags, "log-file", "none")};
		const std::string logFile{logFileArg == "none" ? "" : logFileArg};
		const size_t groupReps{ std::stoul( flagOr(flags, "group-reps", "3") ) };
		const bool skipLD{flags.count("skip-ld") != 0};

		BayesicSpace::MemoryParameters memParams{};
		if (maxMemGB > 0.0) {
			constexpr double bytesPerGB{1073741824.0};
			memParams.maxRAMbytes = static_cast<size_t>(maxMemGB * bytesPerGB);
		}
		BayesicSpace::IndividualAndSketchCounts indivSketches{};
		indivSketches.nIndividuals = nIndividuals;
		indivSketches.kSketches    = kSketches;

		std::cout << "=== vash ldInGroups benchmark ===\n";
		std::cout << "input: " << inputBed << " | individuals: " << nIndividuals
			<< " | hash-size: " << kSketches << " | rows/band: " << nRowsPerBand
			<< " | threads: " << nThreads
			<< " | max-mem(GB): " << (maxMemGB > 0.0 ? std::to_string(maxMemGB) : std::string("auto")) << "\n\n";

		const auto tCtor0 = clockType::now();
		BayesicSpace::GenoTableHash hashTable(inputBed, indivSketches, nThreads, logFile, memParams);
		const auto tCtor1 = clockType::now();
		std::cout << "[time] construction (read .bed + OPH sketching): "
			<< std::fixed << std::setprecision(1) << millisBetween(tCtor0, tCtor1) << " ms\n\n";

		// Time makeLDgroups in isolation to confirm the serial grouping cost. The output file is created
		// only inside ldInGroups (after grouping), so a slow single thread seen after the file appears is
		// downstream of this phase; still, measure it to rule it out.
		std::cout << "--- makeLDgroups (serial grouping), " << groupReps << " rep(s) ---\n";
		std::vector<double> groupMillis;
		groupMillis.reserve(groupReps);
		std::vector<BayesicSpace::HashGroup> lastGroups;
		for (size_t iRep = 0; iRep < groupReps; ++iRep) {
			const auto tGrp0 = clockType::now();
			lastGroups = hashTable.makeLDgroups(nRowsPerBand);
			const auto tGrp1 = clockType::now();
			groupMillis.push_back( millisBetween(tGrp0, tGrp1) );
		}
		if ( !groupMillis.empty() ) {
			const double grpMean{std::accumulate(groupMillis.cbegin(), groupMillis.cend(), 0.0) / static_cast<double>(groupMillis.size())};
			std::cout << "[time] makeLDgroups: mean " << grpMean << " ms (min "
				<< *std::min_element(groupMillis.cbegin(), groupMillis.cend()) << ", max "
				<< *std::max_element(groupMillis.cbegin(), groupMillis.cend()) << ")\n";
		}
		reportGroups(lastGroups);
		std::cout << "\n";

		if (skipLD) {
			std::cout << "--skip-ld set: stopping after grouping.\n";
			return 0;
		}

		// Full ldInGroups run. The in-library [vash-bench] markers (enabled when the library is built
		// with VASH_BENCHMARK) print the per-phase and per-block breakdown to stderr during this call.
		std::cout << "--- ldInGroups end-to-end (grouping + parallel estimate + save) ---\n";
		std::cout << "    (per-phase [vash-bench] lines below come from inside the library, on stderr)\n";
		BayesicSpace::SparsityParameters sparsity{};
		sparsity.nRowsPerBand     = nRowsPerBand;
		sparsity.similarityCutOff = minSimilarity;
		BayesicSpace::InOutFileNames ldFiles{};
		ldFiles.outputFileName = outFile;
		ldFiles.inputFileName  = "";

		const auto tLD0 = clockType::now();
		hashTable.ldInGroups(sparsity, ldFiles);
		const auto tLD1 = clockType::now();
		std::cout << "\n[time] ldInGroups total: " << millisBetween(tLD0, tLD1) << " ms\n";
		std::cout << "output written to: " << outFile << "\n";
	} catch (const std::string &problem) {
		std::cerr << "ERROR: " << problem << "\n";
		return 2;
	} catch (const std::exception &problem) {
		std::cerr << "ERROR: " << problem.what() << "\n";
		return 2;
	}
	return 0;
}
