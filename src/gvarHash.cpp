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

/// Summarize variant tables by hashing
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2021 Anthony J. Greenberg
 * \version 0.6
 *
 * Implementation of classes that take binary variant files and generate lossy summaries with hashing.
 *
 */

#include <ctime>
#include <cstring>
#include <cassert>
#include <string>
#include <fstream>
#include <vector>
#include <array>
#include <utility>  // for std::pair
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <limits>
#include <thread>
#include <immintrin.h>

#include "gvarHash.hpp"
#include "vashFunctions.hpp"
#include "vashParallel.hpp"
#include "vashLogging.hpp"
#include "random.hpp"
#include "similarityMatrix.hpp"
#include "vashBenchmark.hpp"    // phase timers; no-ops unless VASH_BENCHMARK is defined

using namespace BayesicSpace;

namespace {
	/** \brief Draw a seed for a run that was not given one
	 *
	 * The drawn value is limited to the range a signed 32-bit integer can hold, so that the seed a run
	 * records in its log can always be handed back through an integer-typed interface (the `--seed`
	 * command line flag) to repeat that run. The reduced width costs nothing here: the seed only has to
	 * be unpredictable across runs, not cryptographically strong.
	 *
	 * \return the drawn seed
	 */
	uint64_t drawReplayableSeed() {
		return RanDraw().ranInt() % static_cast<uint64_t>( std::numeric_limits<int>::max() );
	}

	/** \brief Band key vector
	 *
	 * Storage for the packed band-hash/locus-index words built by `GenoTableHash::makeLDgroups()`.
	 * The buffer is completely overwritten by the fill that follows its sizing, so the elements are
	 * left uninitialized (see `DefaultInitAllocator`).
	 */
	using BandKeyVector = std::vector< uint64_t, BayesicSpace::DefaultInitAllocator<uint64_t> >;

	/** \brief Width of the locus index in a packed band key */
	constexpr uint32_t bandKeyShift{32};

	/** \brief Extract the band hash from a packed band key
	 *
	 * \param[in] packedKey packed band hash and locus index
	 * \return the band hash
	 */
	uint32_t bandHashOf(const uint64_t &packedKey) noexcept {
		return static_cast<uint32_t>(packedKey >> bandKeyShift);
	}

	/** \brief Extract the locus index from a packed band key
	 *
	 * \param[in] packedKey packed band hash and locus index
	 * \return the locus index
	 */
	uint32_t locusIndexOf(const uint64_t &packedKey) noexcept {
		return static_cast<uint32_t>( packedKey & std::numeric_limits<uint32_t>::max() );
	}

	/** \brief Split a span into ranges, one per worker
	 *
	 * Ranges are floor-sized so that no range bound can overrun the span, with the remainder folded
	 * into the last range. A span shorter than the worker count therefore collapses to a single
	 * non-empty range preceded by empty ones, which the parallel loops that consume the ranges skip.
	 *
	 * \param[in] spanLength number of elements to split
	 * \param[in] nRanges number of ranges (must not be 0)
	 * \return half-open ranges covering the span
	 */
	std::vector< std::pair<size_t, size_t> > makeSpanRanges(const size_t &spanLength, const size_t &nRanges) {
		CountAndSize rangeCounts{0, 0};
		rangeCounts.count = nRanges;
		rangeCounts.size  = spanLength / nRanges;
		std::vector< std::pair<size_t, size_t> > ranges{ makeThreadRanges(rangeCounts) };
		ranges.back().second = spanLength;

		return ranges;
	}

	/** \brief Collect locus groups from sorted band keys
	 *
	 * Turns each run of packed keys sharing a band hash into a group of the loci in that run. Runs of
	 * one locus contribute no pairs and are dropped. The keys must be sorted, which also makes every
	 * group ascending in locus index because the locus index occupies the low bits of a key.
	 *
	 * The scan is split across ranges of the key buffer, with each range boundary pushed forward to
	 * the start of the next run so that no run straddles two ranges and the ranges can be scanned
	 * independently. Groups are concatenated in range order, so the result is ordered by band hash and
	 * does not depend on how the scan was split.
	 *
	 * \param[in] bandKeys sorted packed band hash and locus index words
	 * \param[in] nRanges number of ranges to split the scan into (must not be 0)
	 * \return groups of loci sharing a band hash, each with at least two members
	 */
	std::vector< std::vector<uint32_t> > groupsFromBandKeys(const BandKeyVector &bandKeys, const size_t &nRanges) {
		const size_t nBandKeys{ bandKeys.size() };
		std::vector< std::pair<size_t, size_t> > keyRanges{ makeSpanRanges(nBandKeys, nRanges) };
		for (size_t iRange = 1; iRange < keyRanges.size(); ++iRange) {
			const size_t nominalStart{keyRanges[iRange].first};
			size_t snappedStart{nominalStart};
			// A run longer than the nominal range size leaves the ranges it swallows empty.
			if ( (nominalStart > 0) && (nominalStart < nBandKeys) ) {
				const auto runEndIt = std::upper_bound(
					bandKeys.cbegin() + static_cast<BandKeyVector::difference_type>(nominalStart),
					bandKeys.cend(),
					bandKeys[nominalStart - 1UL],
					[](const uint64_t &lhs, const uint64_t &rhs) {
						return bandHashOf(lhs) < bandHashOf(rhs);
					}
				);
				snappedStart = static_cast<size_t>( runEndIt - bandKeys.cbegin() );
			}
			keyRanges[iRange - 1UL].second = snappedStart;
			keyRanges[iRange].first        = snappedStart;
		}
		std::vector< std::vector< std::vector<uint32_t> > > rangeGroups( keyRanges.size() );
		std::vector<size_t> rangeIndexes( keyRanges.size() );
		std::iota( rangeIndexes.begin(), rangeIndexes.end(), static_cast<size_t>(0) );
		std::for_each(
			parallelPolicy,
			rangeIndexes.cbegin(),
			rangeIndexes.cend(),
			[&rangeGroups, &keyRanges, &bandKeys](const size_t iRange) {
				const size_t rangeEnd{keyRanges[iRange].second};
				size_t runStart{keyRanges[iRange].first};
				while (runStart < rangeEnd) {
					const uint32_t runHash{ bandHashOf(bandKeys[runStart]) };
					size_t runEnd{runStart + 1UL};
					while ( (runEnd < rangeEnd) && (bandHashOf(bandKeys[runEnd]) == runHash) ) {
						++runEnd;
					}
					if ( (runEnd - runStart) >= 2UL ) {
						std::vector<uint32_t> members;
						members.reserve(runEnd - runStart);
						std::transform(
							bandKeys.cbegin() + static_cast<BandKeyVector::difference_type>(runStart),
							bandKeys.cbegin() + static_cast<BandKeyVector::difference_type>(runEnd),
							std::back_inserter(members),
							locusIndexOf
						);
						rangeGroups[iRange].emplace_back( std::move(members) );
					}
					runStart = runEnd;
				}
			}
		);
		std::vector< std::vector<uint32_t> > groups;
		groups.reserve(
			std::accumulate(
				rangeGroups.cbegin(),
				rangeGroups.cend(),
				static_cast<size_t>(0),
				[](const size_t &runningTotal, const std::vector< std::vector<uint32_t> > &eachRangeGroups) {
					return runningTotal + eachRangeGroups.size();
				}
			)
		);
		for (auto &eachRangeGroups : rangeGroups) {
			std::move( eachRangeGroups.begin(), eachRangeGroups.end(), std::back_inserter(groups) );
		}

		return groups;
	}

	/** \brief Drop repeated groups from a sorted group vector
	 *
	 * Removes each group whose hash matches that of the last group kept, leaving the survivors in their
	 * original relative order. Equivalent to `std::unique` under the same hash comparison, but the hash
	 * of every group is computed once, in parallel, instead of twice per adjacent comparison in a serial
	 * pass; hashing dominates the cost, so this is where the parallelism has to go.
	 *
	 * Only groups that are adjacent after hash removal collapse, so the vector must already be sorted
	 * such that identical groups are neighbors (the lexicographic sort in `GenoTableHash::makeLDgroups()`
	 * guarantees this).
	 *
	 * \param[in,out] groups sorted groups, de-duplicated in place
	 * \param[in] hashSeed seed for the group hashes
	 */
	void deduplicateGroups(std::vector< std::vector<uint32_t> > &groups, const uint32_t &hashSeed) {
		if ( groups.empty() ) {
			return;
		}
		std::vector<uint32_t> groupHashes( groups.size() );
		std::transform(
			parallelPolicy,
			groups.cbegin(),
			groups.cend(),
			groupHashes.begin(),
			[&hashSeed](const std::vector<uint32_t> &eachGroup) {
				return murMurHash(eachGroup, hashSeed);
			}
		);
		// Comparing against the last group kept rather than the immediate predecessor matches
		// std::unique, which matters when three or more identical groups are adjacent.
		size_t lastKeptIdx{0};
		for (size_t iGroup = 1; iGroup < groups.size(); ++iGroup) {
			if (groupHashes[iGroup] == groupHashes[lastKeptIdx]) {
				continue;
			}
			++lastKeptIdx;
			if (lastKeptIdx != iGroup) {
				groups[lastKeptIdx]      = std::move(groups[iGroup]);
				groupHashes[lastKeptIdx] = groupHashes[iGroup];
			}
		}
		groups.erase(
			groups.begin() + static_cast<std::vector< std::vector<uint32_t> >::difference_type>(lastKeptIdx + 1UL),
			groups.end()
		);
		groups.shrink_to_fit();
	}

	/** \brief Element budget for a streaming similarity-matrix sink
	 *
	 * Three quarters of half the residual RAM budget (the memory left after the resident genotype
	 * table, established at construction) expressed in matrix elements: the sink splits its reserved
	 * memory 3/4 matrix : 1/4 save-string scratch, and the other half of the residual is left for the
	 * parallel block builders and other work. Never coarser than `nPairs / suggestNchunks`, so a
	 * `suggestNchunks` hint still forces at least that many flushes.
	 *
	 * \param[in] nPairs upper bound on the number of pairs to be processed
	 * \param[in] suggestNchunks minimum number of chunks (flushes) to force
	 * \param[in] workingRAMbytes residual RAM budget (bytes) available for the similarity computation
	 * \return element budget (at least one)
	 */
	// NOLINTNEXTLINE(bugprone-easily-swappable-parameters) counts and a byte budget, but callers pass named locals
	size_t sinkElementBudget(const size_t &nPairs, const size_t &suggestNchunks, const size_t &workingRAMbytes) {
		// 3/4 of (workingRAMbytes/2) bytes, in elements: 3 * workingRAMbytes / (8 * elementSize)
		const size_t autoBudget     = (3UL * workingRAMbytes) / ( 8UL * SimilarityMatrix::elementSize() );
		const size_t clampedNchunks = std::max( suggestNchunks, static_cast<size_t>(1) );
		const size_t forcedBudget   = (nPairs + clampedNchunks - 1UL) / clampedNchunks;      // ceil(nPairs / suggestNchunks)
		return std::max( std::min(autoBudget, forcedBudget), static_cast<size_t>(1) );
	}

	/** \brief Resolve the locus name file for saving
	 *
	 * The names themselves are read downstream by `SimilarityMatrix::save()`, which is handed this
	 * file name; all that is decided here is whether to hand it over. A named but absent `.bim` is not
	 * an error, so the name is cleared and the output falls back to base-1 locus indexes rather than
	 * `save()` failing to read the file.
	 *
	 * \param[in] bimAndOutNames `.bim` and output file names as provided by the caller
	 * \param[in] nLoci number of loci in the table, checked against the `.bim` record count
	 * \param[in,out] logMessages log to record the outcome in
	 * \return the file names to pass downstream, with an unusable `.bim` name removed
	 */
	InOutFileNames resolveLocusNameFile(const InOutFileNames &bimAndOutNames, [[maybe_unused]] const size_t &nLoci, VashLog &logMessages) {
		InOutFileNames resolvedNames{bimAndOutNames};
		if ( resolvedNames.inputFileName.empty() ) {
			return resolvedNames;
		}
		std::fstream bimExistenceTest(resolvedNames.inputFileName, std::ios::in);
		const bool bimExists = bimExistenceTest.good();
		bimExistenceTest.close();
		if (bimExists) {
			logMessages.add("Getting locus names from the " + resolvedNames.inputFileName + " .bim file");
			assert( (getLocusNames(resolvedNames.inputFileName).size() == nLoci) // NOLINT
					&& "ERROR: number of loci in the .bim file not the same as nLoci_");
			return resolvedNames;
		}
		logMessages.add("WARNING: no .bim file " + resolvedNames.inputFileName + "; falling back to locus indexes");
		resolvedNames.inputFileName.clear();

		return resolvedNames;
	}
}

// GenoTableBin methods
constexpr size_t   GenoTableBin::nMagicBytes_    = 3;                // number of leading bytes for .bed files
constexpr uint8_t  GenoTableBin::oneBit_         = 0b00000001;       // One set bit for masking
constexpr uint8_t  GenoTableBin::byteSize_       = 8;                // Size of one byte in bits
constexpr uint8_t  GenoTableBin::bedGenoPerByte_ = 4;                // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableBin::llWordSize_     = 8;                // 64 bit word size in bytes

// Constructors
GenoTableBin::GenoTableBin(const std::string &inputFileName, const uint32_t &nIndividuals, const std::string &logFileName, const size_t &nThreads,
					const MemoryParameters &memParams, const std::optional<uint64_t> &ranSeed)
															: nIndividuals_{nIndividuals}, nThreads_{nThreads}, workingRAMbytes_{0}, locusSeed_{0} {
	if ( !logFileName.empty() ) {
		LogFileNameWithMessage lfMessage;
		lfMessage.logFileName    = logFileName;
		lfMessage.initialMessage = "Genotype binarization from the " + inputFileName + " .bed file";
		logMessages_             = VashLog(lfMessage);
	}
	// The seed is resolved before anything else so that it is the first log entry after the header and
	// every stochastic step below can be replayed from it. Drawing sub-seeds from one master stream keeps
	// the purposes independent: reusing the master directly would give, say, the first locus the same
	// stream as the permutation.
	const uint64_t masterSeed{ ranSeed.value_or( drawReplayableSeed() ) };
	logMessages_.add( "Random number generator seed: " + std::to_string(masterSeed) );
	locusSeed_ = RanDraw(masterSeed).ranInt();
	if (nIndividuals <= 1) {
		logMessages_.add("ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting");
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nIndividuals > std::numeric_limits<size_t>::max() / nIndividuals ) { // a square will overflow
		logMessages_.add("ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too big to make a square relationship matrix; aborting");
		throw std::string("ERROR: the number of individuals (") + std::to_string(nIndividuals) + 
			std::string( ") is too big to make a square relationship matrix in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nThreads_ = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_ = std::max(nThreads_, 1UL);
	const uint32_t nBedBytesPerLocus = (nIndividuals_ / bedGenoPerByte_) + static_cast<uint32_t>( (nIndividuals_ % bedGenoPerByte_) > 0);
	std::fstream inStr;
	// Start by measuring file size
	inStr.open(inputFileName, std::ios::in | std::ios::binary | std::ios::ate);
	if ( inStr.fail() ) {
		logMessages_.add("ERROR: failed to open file " + inputFileName);
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const uint64_t endPosition = static_cast<uint64_t>( inStr.tellg() );
	if (endPosition <= nMagicBytes_) {
		logMessages_.add("ERROR: no genotype records in file " + inputFileName);
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t nTotalBedBytes{static_cast<uint64_t>(endPosition) - nMagicBytes_};
	inStr.close();
	if ( nTotalBedBytes > std::numeric_limits<uint32_t>::max() ) {
		logMessages_.add("ERROR: .bed file (" + inputFileName + ") too large");
		throw std::string("ERROR: there must be fewer than 2^32 bytes in the .bed file ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nLoci_ = static_cast<uint32_t>(nTotalBedBytes) / nBedBytesPerLocus;
	logMessages_.add( "Number of individuals: " + std::to_string(nIndividuals_) );
	logMessages_.add( "Number of loci: "        + std::to_string(nLoci_) );
	logMessages_.add( "Number of threads: "     + std::to_string(nThreads_) );

	inStr.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStr.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	binLocusSize_ = (nIndividuals_ / byteSize_) + static_cast<size_t>( (nIndividuals_ % byteSize_) > 0 );
	// Establish the memory budget before allocating the genotype table.
	// The table allocation is then subtracted from the budget; the remainder bounds the .bed read buffer here
	// and the SimilarityMatrix when estimating LD. If the table alone does not fit, all-by-all similarity is
	// impossible, so fail now.
	const size_t tableBytes = static_cast<size_t>(nLoci_) * binLocusSize_;
	const size_t ramBudget  = memParams.maxRAMbytes > 0 ? memParams.maxRAMbytes : ( 3UL * getAvailableRAM() ) / 4UL;
	if (ramBudget <= tableBytes) {
		logMessages_.add("ERROR: the genotype table (" + std::to_string(tableBytes) + " bytes) does not fit the memory budget (" + std::to_string(ramBudget) + " bytes); aborting");
		throw std::string("ERROR: the genotype table does not fit within the memory budget; raise the limit or reduce the data, in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	workingRAMbytes_ = ramBudget - tableBytes;
	logMessages_.add("Memory budget: "              + std::to_string(ramBudget)        + " bytes");
	logMessages_.add("Genotype table: "             + std::to_string(tableBytes)       + " bytes");
	logMessages_.add("RAM for reading/similarity: " + std::to_string(workingRAMbytes_) + " bytes");
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);
	BedDataStats locusGroupAttributes{};
	locusGroupAttributes.nLociToRead    = std::max( std::min( workingRAMbytes_ / nBedBytesPerLocus, static_cast<size_t>(nLoci_) ), 1UL );   // number of .bed loci to read at a time
	if (memParams.maxLociPerChunk > 0) {                                                                          // optional cap (bounds memory; lets tests force multi-chunk reads)
		locusGroupAttributes.nLociToRead = std::min(locusGroupAttributes.nLociToRead, memParams.maxLociPerChunk);
	}
	locusGroupAttributes.nMemChunks     = nLoci_ / locusGroupAttributes.nLociToRead;
	const size_t remainingLoci          = nLoci_ % locusGroupAttributes.nLociToRead;
	const size_t remainingBytes         = remainingLoci * nBedBytesPerLocus;
	locusGroupAttributes.nBytesToRead   = std::min( locusGroupAttributes.nLociToRead * nBedBytesPerLocus,
													static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) );
	locusGroupAttributes.nLociPerThread = std::max(locusGroupAttributes.nLociToRead / nThreads_, 1UL);
	locusGroupAttributes.nBytesPerLocus = (nIndividuals_ / bedGenoPerByte_) + static_cast<size_t>(nIndividuals_ % bedGenoPerByte_ > 0);
	logMessages_.add(".bed file will be read in " + std::to_string(locusGroupAttributes.nMemChunks) + " chunk(s)");
	assert( ( remainingBytes < std::numeric_limits<std::streamsize>::max() ) //NOLINT
			&& "ERROR: remainingBytes larger than maximum streamsize in GenoTableBin constructor");

	locusGroupAttributes.firstLocusIdx = 0;
	locusGroupAttributes.firstLocusIdx = bed2bin_(locusGroupAttributes, inStr);
	if (remainingLoci > 0) {
		assert( ( remainingBytes < std::numeric_limits<std::streamsize>::max() ) // NOLINT
				&& "ERROR: remainingBytes exceeds maximum streamsize in GenoTableBin constructor" );
		locusGroupAttributes.nLociPerThread = std::max(remainingLoci / nThreads_, 1UL);
		locusGroupAttributes.nBytesToRead   = remainingBytes;
		locusGroupAttributes.nLociToRead    = remainingLoci;
		locusGroupAttributes.nMemChunks     = 1;
		bed2bin_(locusGroupAttributes, inStr);
	}
	inStr.close();
	logMessages_.add("Genotype binarization completed");
}

GenoTableBin::GenoTableBin(const std::vector<int> &maCounts, const uint32_t &nIndividuals, const std::string &logFileName, const size_t &nThreads,
					const MemoryParameters &memParams, const std::optional<uint64_t> &ranSeed)
							: nIndividuals_{nIndividuals}, nLoci_{static_cast<uint32_t>( maCounts.size() / static_cast<size_t>(nIndividuals) )}, nThreads_{nThreads}, workingRAMbytes_{0}, locusSeed_{0} {
	if ( !logFileName.empty() ) {
		LogFileNameWithMessage lfMessage;
		lfMessage.logFileName    = logFileName;
		lfMessage.initialMessage = "Genotype binarization from minor allele count vector";
		logMessages_             = VashLog(lfMessage);
	}
	// The seed is resolved before anything else so that it is the first log entry after the header and
	// every stochastic step below can be replayed from it. Drawing sub-seeds from one master stream keeps
	// the purposes independent: reusing the master directly would give, say, the first locus the same
	// stream as the permutation.
	const uint64_t masterSeed{ ranSeed.value_or( drawReplayableSeed() ) };
	logMessages_.add( "Random number generator seed: " + std::to_string(masterSeed) );
	locusSeed_ = RanDraw(masterSeed).ranInt();
	if ( ( maCounts.size() / static_cast<size_t>(nIndividuals) ) > std::numeric_limits<uint32_t>::max() ) {
		logMessages_.add("ERROR: too many loci");
		throw std::string("ERROR: there must be fewer than 2^32 loci in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (nIndividuals <= 1) {
		logMessages_.add("ERROR: the number of individuals (" + std::to_string(nIndividuals) + ") is too small; aborting");
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( (maCounts.size() % nIndividuals) > 0 ) {
		logMessages_.add( "ERROR: length of allele count vector (" + std::to_string( maCounts.size() ) + " is not divisible by the provided number of individuals (" +
			std::to_string(nIndividuals) );
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(nIndividuals) + std::string( ") in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( maCounts.empty() ) {
		logMessages_.add("ERROR: empty vector of minor allele counts");
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nThreads_ = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_ = std::max(nThreads_, 1UL);

	logMessages_.add( "Number of individuals: " + std::to_string(nIndividuals_) );
	logMessages_.add( "Number of loci: "        + std::to_string(nLoci_) );
	logMessages_.add( "Number of threads: "     + std::to_string(nThreads_) );

	binLocusSize_ = (nIndividuals_ / byteSize_) + static_cast<size_t>( (nIndividuals_ % byteSize_) > 0 );
	// Enforce the memory budget: the resident table must fit within it, leaving room for the
	// SimilarityMatrix estimation. The count vector is caller-owned and not counted here.
	const size_t tableBytes = static_cast<size_t>(nLoci_) * binLocusSize_;
	const size_t ramBudget  = memParams.maxRAMbytes > 0 ? memParams.maxRAMbytes : (3UL * getAvailableRAM()) / 4UL;
	if (ramBudget <= tableBytes) {
		logMessages_.add("ERROR: the genotype table (" + std::to_string(tableBytes) + " bytes) does not fit the memory budget (" + std::to_string(ramBudget) + " bytes); aborting");
		throw std::string("ERROR: the genotype table does not fit within the memory budget; raise the limit or reduce the data, in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	workingRAMbytes_ = ramBudget - tableBytes;
	logMessages_.add("Memory budget: " + std::to_string(ramBudget) + " bytes; genotype table: " + std::to_string(tableBytes) + " bytes; RAM for similarity: " + std::to_string(workingRAMbytes_) + " bytes");
	binGenotypes_.resize(nLoci_ * binLocusSize_, 0);

	// Each locus binarizes into a disjoint slice of binGenotypes_, so the loop is
	// data-parallel. The parallel STL chunks the loci across worker threads; the
	// ThreadCeiling caps concurrency to nThreads_ when it is below the hardware
	// maximum, and is a no-op otherwise (and without a TBB backend).
	const ThreadCeiling threadCeiling(nThreads_);
	std::vector<size_t> locusIndices(nLoci_);
	std::iota(locusIndices.begin(), locusIndices.end(), static_cast<size_t>(0));
	std::for_each(
		parallelPolicy,
		locusIndices.cbegin(),
		locusIndices.cend(),
		[this, &maCounts](const size_t iLocus) {
			LocationWithLength binLocusRange{};
			binLocusRange.start  = iLocus;
			binLocusRange.length = binLocusSize_;
			std::vector<int> macLocus(nIndividuals_);
			std::copy_n(
				std::next( maCounts.cbegin(), static_cast<std::vector<int>::difference_type>(iLocus * nIndividuals_) ),
				nIndividuals_,
				macLocus.begin()
			);
			binarizeMacLocus(macLocus, binLocusRange, binGenotypes_, locusSeed_ + iLocus);
		}
	);
	logMessages_.add("Genotype binarization completed");
}

void GenoTableBin::saveGenoBinary(const std::string &outFileName) const {
	std::fstream out;
	assert( ( binGenotypes_.size() < std::numeric_limits<std::streamsize>::max() ) // NOLINT
			&& "ERROR: binGenotypes_ size exceeds maximum streamsize in GenoTableBin.saveGenoBinary()");
	out.open(outFileName, std::ios::out | std::ios::binary | std::ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), static_cast<std::streamsize>( binGenotypes_.size() ) ); // OK because we are casting to const char*
	out.close();
}

void GenoTableBin::allJaccardLD(const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks) const {
	logMessages_.add("Calculating all pairwise LD using full Jaccard similarity estimates");
	const InOutFileNames outputNames{ resolveLocusNameFile(bimAndLDnames, nLoci_, logMessages_) };

	const size_t nPairs      = nLoci_ * ( nLoci_ - static_cast<size_t>(1) ) / static_cast<size_t>(2);
	const size_t maxElements = sinkElementBudget(nPairs, suggestNchunks, workingRAMbytes_);

	logMessages_.add( "Maximum number of locus pairs held in RAM: " + std::to_string(maxElements) );

	// set up the header
	std::fstream output;
	output.open(outputNames.outputFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\n";
	output.close();

	SimilarityMatrixSink sink(outputNames, WorkloadLimits{nThreads_, maxElements});
	size_t cumChunkIdx{0};
	uint32_t base1chunkIdx{1};
	while (cumChunkIdx < nPairs) {
		const size_t batchSize = std::min(maxElements, nPairs - cumChunkIdx);
		LocationWithLength currStartAndSize{};
		currStartAndSize.start  = cumChunkIdx;
		currStartAndSize.length = batchSize;
		// makeChunkRanges limits the block count to batchSize; save() clamps its own thread count,
		// so nThreads_ can be passed straight through (the ThreadCeiling is an upper bound either way).
		std::vector< std::pair<RowColIdx, RowColIdx> > threadRanges{makeChunkRanges(currStartAndSize, nThreads_)};

		SimilarityMatrix block{
			parallelBuild(
				WorkloadLimits{threadRanges.size(), nThreads_},
				[this, &threadRanges](size_t blockIdx) {
					return jaccardBlock_(threadRanges[blockIdx]);
				}
			)
		};
		logMessages_.add( "\testimated similarity matrix for chunk " + std::to_string(base1chunkIdx) );
		sink.add(block);
		cumChunkIdx += batchSize;
		++base1chunkIdx;
	}
	sink.finalize();
	logMessages_.add("Done calculating and saving all-pair LD");
}

void GenoTableBin::bed2binBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const LocationWithLength &locusSpan) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	size_t begByte{locusSpan.start * binLocusSize_};
	// locusSpan.start is this range's first output locus (see bed2binThreaded_), so iLocus tracks the
	// global locus index the seed must key on
	size_t iLocus{locusSpan.start};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus) {
		LocationWithLength bedWindow{0, 0};
		bedWindow.start  = iBedLocus * locusSpan.length;
		bedWindow.length = locusSpan.length;
		LocationWithLength binWindow{0, 0};
		binWindow.start  = begByte;
		binWindow.length = binLocusSize_;
		binarizeBedLocus(bedWindow, bedData, nIndividuals_, binWindow, binGenotypes_, locusSeed_ + iLocus);
		begByte += binLocusSize_;
		++iLocus;
	}
}

size_t GenoTableBin::bed2binThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const LocationWithLength &locusSpan) {
	// The thread ranges are contiguous and start at 0 (see makeThreadRanges), so a
	// range's output locus index is simply locusSpan.start + bedLocusIndRange.first.
	// Each range writes a disjoint slice of binGenotypes_ and reads a disjoint slice
	// of bedData, so the loop is data-parallel; the ThreadCeiling caps concurrency to
	// nThreads_ (a no-op without a TBB backend).
	const ThreadCeiling threadCeiling(nThreads_);
	std::for_each(
		parallelPolicy,
		threadRanges.cbegin(),
		threadRanges.cend(),
		[this, &bedData, &locusSpan](const std::pair<size_t, size_t> &bedLocusIndRange) {
			LocationWithLength currentLocusSpan{0, 0};
			currentLocusSpan.start  = locusSpan.start + bedLocusIndRange.first;
			currentLocusSpan.length = locusSpan.length;
			bed2binBlk_(bedData, bedLocusIndRange, currentLocusSpan);
		}
	);
	// contiguous-from-0 ranges sum to back().second,
	// so the next free locus index is locusSpan.start + that total.
	return locusSpan.start + threadRanges.back().second;
}

size_t GenoTableBin::bed2bin_(const BedDataStats &locusGroupStats, std::fstream &bedStream) {
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = locusGroupStats.nLociPerThread;
	size_t locusInd    = locusGroupStats.firstLocusIdx;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	assert( (locusGroupStats.nLociToRead >= threadRanges.back().second) // NOLINT
								&& "ERROR: nLociToRead smaller than threadRanges.back().second in bed2bin_" );
	// Extend the last thread's range to cover the remainder loci (nLociToRead is not
	// necessarily divisible by nThreads_); bed2binThreaded_ then writes all
	// nLociToRead loci of this chunk contiguously and returns the next free index.
	threadRanges.back().second = locusGroupStats.nLociToRead;
	std::vector<char> bedChunkToRead(locusGroupStats.nBytesToRead, 0);
	for (size_t iChunk = 0; iChunk < locusGroupStats.nMemChunks; ++iChunk) {
		assert( ( locusGroupStats.nBytesToRead < std::numeric_limits<std::streamsize>::max() ) // NOLINT
								&& "ERROR: nBedBytesToRead exceeds maximum streamsize in bed2bin_" );
		bedStream.read( bedChunkToRead.data(), static_cast<std::streamsize>(locusGroupStats.nBytesToRead) );
		LocationWithLength currentLocusSpan{0, 0};
		currentLocusSpan.start  = locusInd;
		currentLocusSpan.length = locusGroupStats.nBytesPerLocus;
		locusInd                = bed2binThreaded_(bedChunkToRead, threadRanges, currentLocusSpan);
	}
	return locusInd;
}

SimilarityMatrix GenoTableBin::jaccardBlock_(const std::pair<RowColIdx, RowColIdx> &blockRange) const {
	SimilarityMatrix result;
	uint32_t iRow{blockRange.first.iRow};
	if (blockRange.first.iRow != blockRange.second.iRow) {
		for (uint32_t jCol = blockRange.first.jCol; jCol < iRow; ++jCol) { // first, possibly incomplete, row
			RowColIdx localRC{};
			localRC.iRow = iRow;
			localRC.jCol = jCol;
			const JaccardPair localJP{makeJaccardPair_(localRC)};
			result.insert(localRC, localJP);
		}
		++iRow;
		while ( iRow < (blockRange.second.iRow) ) { // complete triangle
			for (uint32_t jCol = 0; jCol < iRow; ++jCol) {
				RowColIdx localRC{};
				localRC.iRow = iRow;
				localRC.jCol = jCol;
				const JaccardPair localJP{makeJaccardPair_(localRC)};
				result.insert(localRC, localJP);
			}
			++iRow;
		}
		for (uint32_t jColRem = 0; jColRem < blockRange.second.jCol; ++jColRem) { // last, possibly incomplete, row (starts at column 0)
			RowColIdx localRC{};
			localRC.iRow = iRow;
			localRC.jCol = jColRem;
			const JaccardPair localJP{makeJaccardPair_(localRC)};
			result.insert(localRC, localJP);
		}
		return result;
	}
	// Range confined to a single row: process only [first.jCol, second.jCol); starting at column 0 here
	// would re-emit pairs owned by earlier ranges and duplicate them across sink flushes.
	for (uint32_t jCol = blockRange.first.jCol; jCol < blockRange.second.jCol; ++jCol) {
		RowColIdx localRC{};
		localRC.iRow = iRow;
		localRC.jCol = jCol;
		const JaccardPair localJP{makeJaccardPair_(localRC)};
		result.insert(localRC, localJP);
	}
	return result;
}

JaccardPair GenoTableBin::makeJaccardPair_(const RowColIdx &rowColumn) const {
	std::vector<uint8_t> locus(binLocusSize_);
	JaccardPair localJP{};
	const size_t rowBin = rowColumn.iRow * binLocusSize_;
	const size_t colBin = rowColumn.jCol * binLocusSize_;
	for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc) {
		locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] & binGenotypes_[colBin + iBinLoc];
	}
	localJP.nIntersect = countSetBits(locus);
	for (size_t iBinLoc = 0; iBinLoc < binLocusSize_; ++iBinLoc) {
		locus[iBinLoc] = binGenotypes_[rowBin + iBinLoc] | binGenotypes_[colBin + iBinLoc];
	}
	localJP.nUnion = countSetBits(locus);
	return localJP;
}

// GenoTableHash methods
constexpr size_t   GenoTableHash::nMagicBytes_    = 3;                                    // number of leading bytes for .bed files
constexpr uint8_t  GenoTableHash::oneBit_         = 0b00000001;                           // One set bit for masking 
constexpr uint8_t  GenoTableHash::byteSize_       = 8;                                    // Size of one byte in bits 
constexpr uint8_t  GenoTableHash::bedGenoPerByte_ = 4;                                    // Number of genotypes in a .bed byte
constexpr uint8_t  GenoTableHash::llWordSize_     = 8;                                    // 64 bit word size in bytes 
constexpr uint32_t GenoTableHash::roundMask_      = 0xfffffff8;                           // mask for rounding down to nearest whole-byte value
constexpr uint64_t GenoTableHash::allBitsSet_     = std::numeric_limits<uint64_t>::max(); // 64-bit word with all bits set
constexpr size_t   GenoTableHash::wordSizeInBits_ = 64;                                   // 64-bit word size
constexpr uint16_t GenoTableHash::emptyBinToken_  = std::numeric_limits<uint16_t>::max(); // Value corresponding to an empty token 

// Constructors
GenoTableHash::GenoTableHash(const std::string &inputFileName, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, const std::string &logFileName,
					const MemoryParameters &memParams, const std::optional<uint64_t> &ranSeed) :
					kSketches_{indivSketchCounts.kSketches},
					nLoci_{0},
					nThreads_{nThreads},
					workingRAMbytes_{0},
					emptyBinIdxSeed_{0},
					locusSeed_{0},
					bandHashSeed_{0} {
	if ( !logFileName.empty() ) {
		LogFileNameWithMessage lfMessage;
		lfMessage.logFileName    = logFileName;
		lfMessage.initialMessage = "Genotype hashing from the " + inputFileName + " .bed file";
		logMessages_             = VashLog(lfMessage);
	}
	// The seed is resolved before anything else so that it is the first log entry after the header and
	// every stochastic step below can be replayed from it. Sub-seeds are drawn from one master stream so
	// the purposes stay independent: reusing the master directly would give, say, the first locus the
	// same stream as the permutation.
	const uint64_t masterSeed{ ranSeed.value_or( drawReplayableSeed() ) };
	logMessages_.add( "Random number generator seed: " + std::to_string(masterSeed) );
	RanDraw prng(masterSeed);
	emptyBinIdxSeed_ = prng.ranInt();
	locusSeed_       = prng.ranInt();
	bandHashSeed_    = static_cast<uint32_t>( prng.ranInt() );
	if (indivSketchCounts.nIndividuals <= 1) {
		logMessages_.add("ERROR: the number of individuals (" + std::to_string(indivSketchCounts.nIndividuals) + ") is too small; aborting");
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ < 3) {
		logMessages_.add("ERROR: number of sketches (" + std::to_string(kSketches_) + ") is too small; aborting");
		throw std::string("ERROR: sketch number must be at least three in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ > indivSketchCounts.nIndividuals) {
		logMessages_.add("ERROR: number of sketches (" + std::to_string(kSketches_) + ") is larger than the number of individuals; aborting");
		throw std::string("ERROR: sketch number must be smaller than the number of individuals in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	// Round up the number of individuals to nearest divisible by kSketches_
	sketchSize_   = (indivSketchCounts.nIndividuals / kSketches_) + static_cast<uint32_t>( (indivSketchCounts.nIndividuals % kSketches_) > 0 );
	nIndividuals_ = sketchSize_ * kSketches_;
	if (indivSketchCounts.kSketches >= emptyBinToken_) {
		logMessages_.add("ERROR: sketch size (" + std::to_string(indivSketchCounts.kSketches) + ") is too big; aborting");
		throw std::string("ERROR: Number of sketches (") + std::to_string(indivSketchCounts.kSketches) + std::string(") implies sketch size (") +
			std::to_string(sketchSize_) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const size_t nBedBytes{(indivSketchCounts.nIndividuals / bedGenoPerByte_) + static_cast<size_t>( (indivSketchCounts.nIndividuals % bedGenoPerByte_) > 0 )};
	nThreads_ = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_ = std::max(nThreads_, 1UL);
	logMessages_.add( "Number of threads used: " + std::to_string(nThreads_) );
	std::fstream inStream;
	// Start by measuring file size
	inStream.open(inputFileName, std::ios::in | std::ios::binary | std::ios::ate);
	if ( inStream.fail() ) {
		logMessages_.add("ERROR: failed to open file " + inputFileName + "; aborting");
		throw std::string("ERROR: failed to open file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const auto endPosition{static_cast<size_t>( inStream.tellg() )};
	if (endPosition <= nMagicBytes_) {
		logMessages_.add("ERROR: no loci in the input .bed file " + inputFileName + "; aborting");
		throw std::string("ERROR: no genotype records in file ") + inputFileName + std::string(" in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	inStream.close();
	const size_t fileSize{endPosition - nMagicBytes_};
	const size_t tmpNloci{fileSize / nBedBytes};
	if ( tmpNloci > std::numeric_limits<uint32_t>::max() ) {
		logMessages_.add( "ERROR: too many loci (" + std::to_string(tmpNloci) );
		throw std::string("ERROR: there must be fewer than 2^32 loci in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	nLoci_ = static_cast<uint32_t>(tmpNloci);

	logMessages_.add( "Number of individuals: "         + std::to_string(indivSketchCounts.nIndividuals) );
	logMessages_.add( "Number of individuals to hash: " + std::to_string(nIndividuals_) );
	logMessages_.add( "Number of loci: "                + std::to_string(nLoci_) );
	logMessages_.add( "Hash size: "                     + std::to_string(kSketches_) );

	locusSize_       = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;                    // round up to the nearest multiple of 8
	nFullWordBytes_  = (nIndividuals_ - 1) / byteSize_;
	// Establish the memory budget before allocating the hashed genotype table.
	// The table allocation is then subtracted; the remainder bounds the .bed read buffer here and the
	// SimilarityMatrix size for LD estimation. If the table alone does not fit, fail now, because our
	// implementation does not allow for all potential pairs to be considered is the whole hashed genotype table
	// is not in RAM.
	const size_t tableBytes = static_cast<size_t>(kSketches_) * nLoci_ * sizeof(uint16_t);
	const size_t ramBudget  = memParams.maxRAMbytes > 0 ? memParams.maxRAMbytes : (3UL * getAvailableRAM()) / 4UL;
	if (ramBudget <= tableBytes) {
		logMessages_.add("ERROR: the genotype table (" + std::to_string(tableBytes) + " bytes) does not fit the memory budget (" + std::to_string(ramBudget) + " bytes); aborting");
		throw std::string("ERROR: the genotype table does not fit within the memory budget; raise the limit or reduce the data, in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	workingRAMbytes_ = ramBudget - tableBytes;
	logMessages_.add("Memory budget: "              + std::to_string(ramBudget)        + " bytes");
	logMessages_.add("Genotype table: "             + std::to_string(tableBytes)       + " bytes");
	logMessages_.add("RAM for reading/similarity: " + std::to_string(workingRAMbytes_) + " bytes");
	sketches_.resize(static_cast<size_t>(kSketches_) * nLoci_, emptyBinToken_);
	inStream.open(inputFileName, std::ios::in | std::ios::binary);
	std::array<char, nMagicBytes_> magicBuf{0};
	inStream.read( magicBuf.data(), magicBuf.size() );
	testBedMagicBytes(magicBuf);
	// Generate the binary genotype table while reading the .bed file
	BedDataStats locusGroupAttributes{};
	locusGroupAttributes.nBytesPerLocus = (indivSketchCounts.nIndividuals / bedGenoPerByte_) + static_cast<size_t>(indivSketchCounts.nIndividuals % bedGenoPerByte_ > 0);
	locusGroupAttributes.nLociToRead    = std::max( std::min( workingRAMbytes_ / locusGroupAttributes.nBytesPerLocus, static_cast<size_t>(nLoci_) ), 1UL );   // number of .bed loci to read at a time
	if (memParams.maxLociPerChunk > 0) {                                                                                                                      // optional cap (bounds memory; lets tests force multi-chunk reads)
		locusGroupAttributes.nLociToRead = std::min(locusGroupAttributes.nLociToRead, memParams.maxLociPerChunk);
	}
	const size_t remainingLoci          = nLoci_ % locusGroupAttributes.nLociToRead;
	const size_t remainingBytes         = remainingLoci * locusGroupAttributes.nBytesPerLocus;
	locusGroupAttributes.nMemChunks     = nLoci_ / locusGroupAttributes.nLociToRead;
	locusGroupAttributes.nBytesToRead   = std::min( locusGroupAttributes.nLociToRead * locusGroupAttributes.nBytesPerLocus,
													static_cast<size_t>( std::numeric_limits<std::streamsize>::max() ) );
	locusGroupAttributes.nLociPerThread = locusGroupAttributes.nLociToRead / nThreads_;

	logMessages_.add(".bed file will be read in " + std::to_string(locusGroupAttributes.nMemChunks) + " chunk(s)");

	// Sample with replacement additional individuals to pad out the total
	std::vector< std::pair<size_t, size_t> > addIndv;
	for (size_t iAddIndiv = indivSketchCounts.nIndividuals; iAddIndiv < nIndividuals_; ++iAddIndiv) {
		addIndv.emplace_back( iAddIndiv, prng.sampleInt(indivSketchCounts.nIndividuals) );
	}
	if ( !addIndv.empty() ) {
		std::string addIndexes;
		for (const auto &[originalIdx, sampledIdx] : addIndv) {
			addIndexes += std::to_string(sampledIdx) + " ";
		}
		logMessages_.add("Re-sampled individuals: " + addIndexes);
	}
	// generate the sequence of random integers; each column must be permuted the same
	const std::vector<size_t> ranInts{prng.fyIndexesUp(nIndividuals_)};

	locusGroupAttributes.firstLocusIdx = 0;
	locusGroupAttributes.firstLocusIdx = bed2oph_(locusGroupAttributes, inStream, ranInts, addIndv);
	if (remainingLoci > 0) {
		locusGroupAttributes.nLociPerThread = std::max(remainingLoci / nThreads_, 1UL);
		locusGroupAttributes.nBytesToRead   = remainingBytes;
		locusGroupAttributes.nLociToRead    = remainingLoci;
		locusGroupAttributes.nMemChunks     = 1;
		bed2oph_(locusGroupAttributes, inStream, ranInts, addIndv);
	}
	inStream.close();
	logMessages_.add("Genotype hashing completed");
}

GenoTableHash::GenoTableHash(const std::vector<int> &maCounts, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, const std::string &logFileName,
					const MemoryParameters &memParams, const std::optional<uint64_t> &ranSeed) :
								nIndividuals_{indivSketchCounts.nIndividuals},
								kSketches_{indivSketchCounts.kSketches},
								nLoci_{static_cast<uint32_t>(maCounts.size() / indivSketchCounts.nIndividuals)},
								nThreads_{nThreads},
								workingRAMbytes_{0},
								emptyBinIdxSeed_{0},
								locusSeed_{0},
								bandHashSeed_{0} {
	if ( !logFileName.empty() ) {
		LogFileNameWithMessage lfMessage;
		lfMessage.logFileName    = logFileName;
		lfMessage.initialMessage = "Genotype hashing from a minor allele count vector";
		logMessages_             = VashLog(lfMessage);
	}
	// The seed is resolved before anything else so that it is the first log entry after the header and
	// every stochastic step below can be replayed from it. Sub-seeds are drawn from one master stream so
	// the purposes stay independent: reusing the master directly would give, say, the first locus the
	// same stream as the permutation.
	const uint64_t masterSeed{ ranSeed.value_or( drawReplayableSeed() ) };
	logMessages_.add( "Random number generator seed: " + std::to_string(masterSeed) );
	RanDraw prng(masterSeed);
	emptyBinIdxSeed_ = prng.ranInt();
	locusSeed_       = prng.ranInt();
	bandHashSeed_    = static_cast<uint32_t>( prng.ranInt() );
	if (indivSketchCounts.nIndividuals <= 1) {
		logMessages_.add("ERROR: the number of individuals (" + std::to_string(indivSketchCounts.nIndividuals) + ") is too small; aborting");
		throw std::string("ERROR: number of individuals must be greater than 1 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( (maCounts.size() % indivSketchCounts.nIndividuals) > 0) {
		logMessages_.add("ERROR: minor allele vector size (" + std::to_string( maCounts.size() ) + ") is not evenly divisible by the number of individuals (" +
							std::to_string(nIndividuals_) + "); aborting");
		throw std::string("ERROR: length of allele count vector (") + std::to_string( maCounts.size() ) + std::string(" is not divisible by the provided number of individuals (") +
			std::to_string(indivSketchCounts.nIndividuals) + std::string( ") in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if ( maCounts.empty() ) {
		logMessages_.add("ERROR: minor allele count vector is empty; aborting");
		throw std::string("ERROR: empty vector of minor allele counts in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ < 3) {
		logMessages_.add("ERROR: sketch size (" + std::to_string(kSketches_) + ") is too small; aborting");
		throw std::string("ERROR: sketch size must be at least three in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (kSketches_ > indivSketchCounts.nIndividuals) {
		logMessages_.add("ERROR: number of sketches (" + std::to_string(kSketches_) + ") is larger than the number of individuals; aborting");
		throw std::string("ERROR: sketch number must be smaller than the number of individuals in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}

	nThreads_ = std::min( nThreads_, static_cast<size_t>( std::thread::hardware_concurrency() ) );
	nThreads_ = std::max(nThreads_, 1UL);
	logMessages_.add( "Number of threads used: " + std::to_string(nThreads_) );

	sketchSize_   = (indivSketchCounts.nIndividuals / kSketches_) + static_cast<uint16_t>( (indivSketchCounts.nIndividuals % kSketches_) > 0 );
	nIndividuals_ = sketchSize_ * kSketches_;
	if (indivSketchCounts.kSketches >= emptyBinToken_) {
		logMessages_.add("ERROR: sketch size (" + std::to_string(indivSketchCounts.kSketches) + ") is too small; aborting");
		throw std::string("ERROR: Number of sketches (") + std::to_string(kSketches_) + std::string(") implies sketch size (") +
			std::to_string(indivSketchCounts.kSketches) + std::string(") that is larger than ") + std::to_string(emptyBinToken_) +
			std::string( ", the largest allowed value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	// Sample with replacement additional individuals to pad out the total
	std::vector< std::pair<size_t, size_t> > addIndv;
	for (size_t iAddIndiv = indivSketchCounts.nIndividuals; iAddIndiv < nIndividuals_; ++iAddIndiv) {
		addIndv.emplace_back( iAddIndiv, prng.sampleInt(indivSketchCounts.nIndividuals) );
	}
	if ( !addIndv.empty() ) {
		std::string addIndexes;
		for (const auto &[originalIdx, sampledIdx] : addIndv) {
			addIndexes += std::to_string(sampledIdx) + " ";
		}
		logMessages_.add("Re-sampled individuals: " + addIndexes);
	}
	locusSize_      = ( ( nIndividuals_ + (byteSize_ - 1) ) & roundMask_ ) / byteSize_;   // round up to the nearest multiple of 8
	nFullWordBytes_ = (nIndividuals_ - 1) / byteSize_;
	// Enforce the memory budget: the resident sketch table must fit within it, leaving room for the
	// SimilarityMatrix for LD estimation. The count vector is caller-owned and not counted here.
	const size_t tableBytes = static_cast<size_t>(kSketches_) * nLoci_ * sizeof(uint16_t);
	const size_t ramBudget  = memParams.maxRAMbytes > 0 ? memParams.maxRAMbytes : (3UL * getAvailableRAM()) / 4UL;
	if (ramBudget <= tableBytes) {
		logMessages_.add("ERROR: the genotype table (" + std::to_string(tableBytes) + " bytes) does not fit the memory budget (" + std::to_string(ramBudget) + " bytes); aborting");
		throw std::string("ERROR: the genotype table does not fit within the memory budget; raise the limit or reduce the data, in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	workingRAMbytes_ = ramBudget - tableBytes;
	logMessages_.add("Memory budget: " + std::to_string(ramBudget) + " bytes; genotype table: " + std::to_string(tableBytes) + " bytes; RAM for similarity: " + std::to_string(workingRAMbytes_) + " bytes");
	sketches_.resize(static_cast<size_t>(kSketches_) * nLoci_, emptyBinToken_);
	// generate the sequence of random integers; each column must be permuted the same
	std::vector<size_t> ranInts{prng.fyIndexesUp(nIndividuals_)};

	logMessages_.add( "Number of individuals: "  + std::to_string(nIndividuals_) );
	logMessages_.add( "Number of loci: "         + std::to_string(nLoci_) );
	logMessages_.add( "Hash size: "              + std::to_string(kSketches_) );

	const size_t nLociPerThread = nLoci_ / nThreads_;
	if (nLociPerThread == 0) {
		const std::pair<size_t, size_t> allLoci{0, nLoci_};
		mac2ophBlk_(maCounts, allLoci, ranInts, addIndv);
		return;
	}
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = nLociPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	threadRanges.back().second = nLoci_;
	// Each range processes a disjoint set of loci, writing disjoint sketches_ slices via mac2ophBlk_
	// (local binLocus/macLocus, read-only ranInts/addIndv), so the loop is data-parallel; the
	// ThreadCeiling caps concurrency to nThreads_ (a no-op without a TBB backend).
	const ThreadCeiling threadCeiling(nThreads_);
	std::for_each(
		parallelPolicy,
		threadRanges.cbegin(),
		threadRanges.cend(),
		[this, &maCounts, &ranInts, &addIndv](const std::pair<size_t, size_t> &eachTR) {
			mac2ophBlk_(maCounts, eachTR, ranInts, addIndv);
		}
	);
	logMessages_.add("Genotype hashing completed");
}

void GenoTableHash::allHashLD(const float &similarityCutOff, const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks) const {
	std::vector<uint32_t> allLocusIndexes(nLoci_);
	std::iota(allLocusIndexes.begin(), allLocusIndexes.end(), 0);

	const size_t nPairs      = static_cast<size_t>(nLoci_) * (static_cast<size_t>(nLoci_) - 1UL) / 2UL;
	const size_t maxElements = sinkElementBudget(nPairs, suggestNchunks, workingRAMbytes_);

	const InOutFileNames outputNames{ resolveLocusNameFile(bimAndLDnames, nLoci_, logMessages_) };

	logMessages_.add("Calculating all pairwise LD");
	logMessages_.add( "Maximum number of locus pairs held in RAM: " + std::to_string(maxElements) );

	std::fstream output;
	output.open(outputNames.outputFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\n";
	output.close();

	SimilarityMatrixSink sink(outputNames, WorkloadLimits{nThreads_, maxElements});
	size_t cumChunkIdx{0};
	while (cumChunkIdx < nPairs) {
		const size_t batchSize = std::min(maxElements, nPairs - cumChunkIdx);
		LocationWithLength currStartAndSize{};
		currStartAndSize.start  = cumChunkIdx;
		currStartAndSize.length = batchSize;

		// makeChunkRanges limits the block count to batchSize; save() clamps its own thread count,
		// so nThreads_ can be passed straight through (the ThreadCeiling is an upper bound either way).
		std::vector< std::pair<RowColIdx, RowColIdx> > threadRanges{makeChunkRanges(currStartAndSize, nThreads_)};
		SimilarityMatrix block{
			parallelBuild(
				WorkloadLimits{threadRanges.size(), nThreads_},
				[this, &threadRanges, &allLocusIndexes, &similarityCutOff](size_t blockIdx) {
					return hashJacBlock_(threadRanges[blockIdx], allLocusIndexes, similarityCutOff);
				}
			)
		};
		sink.add(block);

		cumChunkIdx += batchSize;
	}
	sink.finalize();
	logMessages_.add("All pairwise LD calculated and saved");
}

std::vector<HashGroup> GenoTableHash::makeLDgroups(const size_t &nRowsPerBand) const {
	assert( (nRowsPerBand != 0) // NOLINT
			&& "ERROR: nRowsPerBand must not be 0 in makeLDgroups()" );
	assert( (nRowsPerBand < kSketches_) // NOLINT
			&& "ERROR: nRowsPerBand must be less than kSketches_ in makeLDgroups()" );
	const size_t nBands = kSketches_ / nRowsPerBand;                                       // only using full-size bands because smaller ones permit inclusion of low-similarity pairs
	assert( ( nBands <= std::numeric_limits<uint16_t>::max() ) // NOLINT
			&& "ERROR: number of bands cannot exceed uint16_t max in makeLDgroups()" );

	logMessages_.add("Grouping loci");
	logMessages_.add( "Number of rows per band: " + std::to_string(nRowsPerBand) );
	logMessages_.add( "Number of bands: "         + std::to_string(nBands) );

	VASH_BENCH_TP(vashBenchLDgroups);

	// Grouping is nothing more than collecting the loci that share a band hash, so instead of routing
	// the hashes through a shared hash table -- which serializes the pass and pays a node allocation
	// per insertion -- each hash goes into the high half of a 64-bit word with its locus index in the
	// low half, and the words are sorted. Equal hashes then form contiguous runs, and because the
	// locus index occupies the low bits every run comes out ascending in locus order, which is the
	// ordering the group sort, the de-duplication and the std::set_union below all rely on.
	// The band index is part of the hashed key, so one sort separates the bands without the bands ever
	// being handled separately: parallelism comes from the locus count rather than the band count,
	// which matters because there are usually far fewer bands than threads. Eight bytes per
	// (locus, band) pair is also less than the hash-table node it replaces, so no input that could be
	// grouped before can fail to be grouped now.
	const size_t nBandKeys{ static_cast<size_t>(nLoci_) * nBands };
	logMessages_.add( "Band key buffer (bytes): " + std::to_string( nBandKeys * sizeof(uint64_t) ) );
	BandKeyVector bandKeys(nBandKeys);                                                                        // sized but not initialized; the fill below writes every element

	VASH_BENCH_NOTE("makeLDgroups: band keys", nBandKeys);
	VASH_BENCH_LAP("makeLDgroups: band key buffer allocation", vashBenchLDgroups);
	// ThreadCeiling caps concurrency to nThreads_ (a no-op without a TBB backend).
	const ThreadCeiling threadCeiling(nThreads_);
	const std::vector< std::pair<size_t, size_t> > locusRanges{ makeSpanRanges(nLoci_, nThreads_) };
	// Each range writes a disjoint slice of bandKeys and only reads sketches_, so the fill needs no
	// synchronization; the band vector is hoisted out of the loops to keep it allocation-free.
	std::for_each(
		parallelPolicy,
		locusRanges.cbegin(),
		locusRanges.cend(),
		[this, &bandKeys, &nBands, &nRowsPerBand](const std::pair<size_t, size_t> &eachRange) {
			std::vector<uint16_t> bandVec;
			bandVec.reserve(nRowsPerBand + 1UL);
			for (size_t iLocus = eachRange.first; iLocus < eachRange.second; ++iLocus) {
				for (uint16_t iBand = 0; iBand < static_cast<uint16_t>(nBands); ++iBand) {
					bandVec.clear();
					bandVec.push_back(iBand);                                                                 // add the band index to the hash, so that only corresponding bands are compared

					const auto firstSketchIt = sketches_.cbegin()
						+ static_cast<std::vector<uint16_t>::difference_type>( (iLocus * kSketches_) + (iBand * nRowsPerBand) );
					const auto lastSketchIt = firstSketchIt
						+ static_cast<std::vector<uint16_t>::difference_type>(nRowsPerBand);
					std::copy( firstSketchIt, lastSketchIt, std::back_inserter(bandVec) );

					LocationWithLength bandVecWindow{0, 0};
					bandVecWindow.start  = 0;
					bandVecWindow.length = bandVec.size();
					const uint32_t hash  = murMurHash(bandVec, bandVecWindow, bandHashSeed_);
					bandKeys[(iLocus * nBands) + iBand] =
						(static_cast<uint64_t>(hash) << bandKeyShift) | static_cast<uint64_t>(iLocus);
				}
			}
		}
	);
	VASH_BENCH_LAP("makeLDgroups: band key fill (parallel hash)", vashBenchLDgroups);
	parallelSort( bandKeys.begin(), bandKeys.end() );
	VASH_BENCH_LAP("makeLDgroups: band key sort (parallel)", vashBenchLDgroups);

	std::vector< std::vector<uint32_t> > groups{groupsFromBandKeys(bandKeys, nThreads_)};
	bandKeys.clear();
	bandKeys.shrink_to_fit();                                                                                 // the keys are dead from here on and the buffer is the largest thing alive
	VASH_BENCH_NOTE("makeLDgroups: groups before de-duplication", groups.size());
	VASH_BENCH_LAP("makeLDgroups: run scan into groups (parallel)", vashBenchLDgroups);

	// pre-sort the groups by position
	// this carries some overhead, but speeds the downstream pair sorting
	// enough that overall execution timing is comparable.
	// It also enables processing by chunks if the whole sparse table does not fit in RAM
	// The comparison is lexicographic over the whole group.
	// A sort on group subsets (I had first two elements before) would be faster
	// but cannot guarantee that identical groups will end up adjacent.
	parallelSort(
		groups.begin(),
		groups.end(),
		[](const std::vector<uint32_t> &group1, const std::vector<uint32_t> &group2) {
			return std::lexicographical_compare( group1.cbegin(), group1.cend(), group2.cbegin(), group2.cend() );
		}
	);
	VASH_BENCH_LAP("makeLDgroups: sort groups (lexicographic, parallel)", vashBenchLDgroups);
	// de-duplicate the groups
	logMessages_.add( "Number of groups before de-duplication: " + std::to_string( groups.size() ) );
	deduplicateGroups(groups, bandHashSeed_);
	VASH_BENCH_LAP("makeLDgroups: de-duplicate groups", vashBenchLDgroups);
	logMessages_.add( "Number of groups after de-duplication: " + std::to_string( groups.size() ) );

	// Each surviving group is kept as it came out of banding. Consecutive groups sharing a first locus
	// were once merged into their union to de-fragment the grouping, but that emits every cross pair
	// between the two, and those loci share no band: they are outside the candidate set that banding
	// selected, so their similarity was reported on the strength of a shared smallest index alone. The
	// merge was also arbitrary in which groups it reached, since only the first (smallest) index was
	// compared -- two groups sharing any other locus were left apart.
	std::vector<HashGroup> indexedGroups;
	indexedGroups.reserve( groups.size() );
	uint64_t cumulativeNpairs{0};
	for (auto &eachGroup : groups) {
		// running total, so the last element carries the pair count over all groups
		cumulativeNpairs += eachGroup.size() * ( eachGroup.size() - 1 ) / 2;
		indexedGroups.emplace_back( HashGroup{cumulativeNpairs, std::move(eachGroup)} );
	}

	return indexedGroups;
}

void GenoTableHash::makeLDgroups(const size_t &nRowsPerBand, const InOutFileNames &bimAndGroupNames) const {
	const std::vector<HashGroup> ldGroups{this->makeLDgroups(nRowsPerBand)};
	logMessages_.add("Saving group IDs only");
	std::vector<std::string> locusNames{};
	if ( !bimAndGroupNames.inputFileName.empty() ) {
		std::fstream bimExistenceTest(bimAndGroupNames.inputFileName, std::ios::in);
		const bool bimExists = bimExistenceTest.good();
		bimExistenceTest.close();
		// A named but absent .bim is not an error: the output falls back to base-1 indexes. The
		// locus-count check therefore only applies when names were actually read.
		if (bimExists) {
			logMessages_.add("Getting locus names from the " + bimAndGroupNames.inputFileName + " .bim file");
			locusNames = getLocusNames(bimAndGroupNames.inputFileName);
			assert( (locusNames.size() == nLoci_) // NOLINT
					&& "ERROR: number of loci in the .bim file not the same as nLoci_");
		} else {
			logMessages_.add("WARNING: no .bim file " + bimAndGroupNames.inputFileName + "; falling back to locus indexes");
		}
	}

	std::fstream out;
	out.open(bimAndGroupNames.outputFileName, std::ios::out | std::ios::trunc);
	out << "groupID\tlocusIdx\n";
	uint32_t groupID{1};
	if ( locusNames.empty() ) {
		for (const auto &eachGroup : ldGroups) {
			for (const auto &locusIdx : eachGroup.locusIndexes) {
				out << "G" << groupID << "\t" << locusIdx + 1 << "\n";
			}
			++groupID;
		}
		out.close();
		logMessages_.add("Finished saving group IDs");
		return;
	}
	for (const auto &eachGroup : ldGroups) {
		for (const auto &locusIdx : eachGroup.locusIndexes) {
			out << "G" << groupID << "\t" << locusNames[locusIdx] << "\n";
		}
		++groupID;
	}
	out.close();
	logMessages_.add("Finished saving group IDs");
}

void GenoTableHash::ldInGroups(const SparsityParameters &sparsityValues, const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks) const {
	VASH_BENCH_TP(vashLDinGroups);
	std::vector<HashGroup> ldGroups{this->makeLDgroups(sparsityValues.nRowsPerBand)};
	VASH_BENCH_LAP("ldInGroups: makeLDgroups (serial grouping)", vashLDinGroups);

	const InOutFileNames outputNames{ resolveLocusNameFile(bimAndLDnames, nLoci_, logMessages_) };

	// makeLDgroups() returns nothing when no bucket collects two or more loci (e.g. every locus is
	// unique at the chosen band width). There is then no pair to estimate, so emit the header alone;
	// the group traversal below indexes back() and cbegin() unconditionally.
	if ( ldGroups.empty() ) {
		logMessages_.add("No LD groups with more than one locus; saving an empty similarity matrix");
		std::fstream emptyOutput;
		emptyOutput.open(outputNames.outputFileName, std::ios::trunc | std::ios::out);
		emptyOutput << "locus1\tlocus2\tjaccard\n";
		emptyOutput.close();
		return;
	}

	const size_t totalPairNumber{ldGroups.back().cumulativeNpairs};    // total number of pairs
	logMessages_.add("Estimating LD in groups");
	logMessages_.add( "number of pairs in the hash table: " + std::to_string(totalPairNumber) );

	const size_t maxElements = sinkElementBudget(totalPairNumber, suggestNchunks, workingRAMbytes_);
	logMessages_.add( "Maximum number of locus pairs held in RAM: " + std::to_string(maxElements) );

	std::fstream output;
	output.open(outputNames.outputFileName, std::ios::trunc | std::ios::out);
	output << "locus1\tlocus2\tjaccard\n";
	output.close();

	SimilarityMatrixSink sink(outputNames, WorkloadLimits{nThreads_, maxElements});
	HashGroupItPairCount startPair{};
	startPair.hgIterator = ldGroups.cbegin();
	startPair.pairCount  = 0;
	const size_t lastPairNumber{ldGroups.back().locusIndexes.size() * (ldGroups.back().locusIndexes.size() - 1) / 2};
	uint32_t base1chunkIdx{1};
	// Consume the groups in pair-batches capped at the element budget, so each parallelBuild result
	// fits the sink's reserved buffer. The sink accumulates batches and flushes as the budget fills;
	// because saved pairs can no longer be de-duplicated, it flushes as late as possible to limit the
	// cross-flush duplication that overlapping groups can introduce.
	bool done{false};
	while (!done) {
		std::vector< std::pair<HashGroupItPairCount, HashGroupItPairCount> > groupRanges;
		// Over-decompose into more blocks than threads so parallelBuild's work-stealing can
		// balance uneven group sizes: with one block per thread the heaviest block sets the
		// wall time, but finer blocks let idle threads pick up the slack.
		constexpr size_t blockOverDecomposition{4};
		const size_t nBlocks{std::min(blockOverDecomposition * nThreads_, maxElements)};
		const std::vector<size_t> threadSizes{makeChunkSizes( maxElements, nBlocks )};
		groupRanges.reserve( threadSizes.size() );
		for (const auto &eachThrSize : threadSizes) {
			groupRanges.emplace_back( makeGroupRanges(ldGroups, startPair, eachThrSize) );
			startPair = groupRanges.back().second;
		}
#ifdef VASH_BENCHMARK
		std::cerr << "[vash-bench] ldInGroups: chunk " << base1chunkIdx << " over " << groupRanges.size()
			<< " blocks (cap " << nThreads_ << " threads)\n";
#endif
		VASH_BENCH_TP(vashChunk);
		SimilarityMatrix block{
			parallelBuild(
				WorkloadLimits{groupRanges.size(), nThreads_},
				[this, &groupRanges, &sparsityValues](size_t blockIdx) {
					return hashJacBlock_(groupRanges[blockIdx], sparsityValues.similarityCutOff);
				}
			)
		};
		VASH_BENCH_LAP("ldInGroups: chunk parallelBuild (parallel estimate + serial consolidate)", vashChunk);
		logMessages_.add( "\tfinished similarity matrix estimation for chunk " + std::to_string(base1chunkIdx) );
		sink.add(block);
		VASH_BENCH_LAP("ldInGroups: chunk sink.add (serial merge, may flush)", vashChunk);
		++base1chunkIdx;
		done = ( startPair.hgIterator == ldGroups.cend() )
			|| ( ( startPair.hgIterator == std::prev( ldGroups.cend() ) ) && (startPair.pairCount == lastPairNumber) );
	}
	VASH_BENCH_TP(vashFinalize);
	sink.finalize();
	VASH_BENCH_LAP("ldInGroups: sink.finalize (final flush + save)", vashFinalize);
	logMessages_.add("Finished calculating and saving LD in groups");
}

void GenoTableHash::permuteBits_(const std::vector<size_t> &permutationIdx, std::vector<uint8_t> &binLocus) const {
	size_t iIndiv = 0;
	size_t iByte  = 0;
	while(iByte < nFullWordBytes_) {
		for (uint8_t iInLocusByte = 0; iInLocusByte < byteSize_; ++iInLocusByte) {
			auto bytePair            = static_cast<uint16_t>(binLocus[iByte]);
			const size_t perIndiv    = permutationIdx[iIndiv++];                                                           // post-increment to use current value for index first
			const size_t permByteInd = perIndiv / byteSize_;
			const auto permInByteInd = static_cast<uint8_t>( perIndiv - (perIndiv & roundMask_) );
			// Pair the current locus byte with the byte containing the value to be swapped
			// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
			bytePair                |= static_cast<uint16_t>(binLocus[permByteInd] << byteSize_);
			const auto mask          = static_cast<uint16_t>(oneBit_ << iInLocusByte);
			const auto perMask       = static_cast<uint8_t>(oneBit_ << permInByteInd);
			const auto shiftDistance = static_cast<uint16_t>( (byteSize_ - iInLocusByte) + permInByteInd );                // subtraction is safe b/c byteSize is the loop terminator
			const auto temp1         = static_cast<uint16_t>( ( bytePair ^ (bytePair >> shiftDistance) ) & mask );
			const auto temp2         = static_cast<uint16_t>(temp1 << shiftDistance);
			bytePair                ^= temp1 ^ temp2;
			// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus)
			// Must modify the current byte in each loop iteration because permutation indexes may fall into it
			binLocus[iByte]       ^= static_cast<uint8_t>( ( binLocus[iByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask) );
			binLocus[permByteInd] ^= static_cast<uint8_t>( ( binLocus[permByteInd] ^ static_cast<uint8_t>(bytePair >> byteSize_) ) & perMask );
 		}
		++iByte;
	}
	// Finish the individuals in the remaining partial byte, if any
	uint8_t iInLocusByte = 0;
	while (iIndiv < nIndividuals_ - 1) {
		auto bytePair            = static_cast<uint16_t>(binLocus[iByte]);
		const size_t perIndiv    = permutationIdx[iIndiv++];                                   // post-increment to use current value for index first
		const size_t permByteInd = perIndiv / byteSize_;
		const auto permInByteInd = static_cast<uint8_t>( perIndiv - (perIndiv & roundMask_) );
		// Pair the current locus byte with the byte containing the value to be swapped
		// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
		bytePair                |= static_cast<uint16_t>(binLocus[permByteInd] << byteSize_);
		const auto mask          = static_cast<uint16_t>(oneBit_ << iInLocusByte);
		const auto perMask       = static_cast<uint8_t>(oneBit_ << permInByteInd);
		const auto shiftDistance = static_cast<uint16_t>( (byteSize_ - iInLocusByte) + permInByteInd );           // subtraction is safe b/c byteSize is the loop terminator
		const auto temp1         = static_cast<uint16_t>( ( bytePair ^ (bytePair >> shiftDistance) ) & mask );
		const auto temp2         = static_cast<uint16_t>(temp1 << shiftDistance);
		bytePair                ^= temp1 ^ temp2;
		// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus)
		// Must modify the current byte in each loop iteration because permutation indexes may fall into it
		binLocus[iByte]       ^= static_cast<uint8_t>( ( binLocus[iByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask) );
		binLocus[permByteInd] ^= static_cast<uint8_t>( ( binLocus[permByteInd] ^ static_cast<uint8_t>(bytePair >> byteSize_) ) & perMask );
		++iInLocusByte;
	}
}

void GenoTableHash::locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint8_t> &binLocus) {
	// Start with a permutation to make OPH
	permuteBits_(permutation, binLocus);
	// Now make the sketches
	std::vector<size_t> filledIndexes;                                                       // indexes of the non-empty sketches
	size_t iByte{0};
	const size_t sketchBeg{locusInd * kSketches_};
	size_t iSketch{0};
	uint64_t sketchTail{0};                                                                  // left over buts from beyond the last full byte of the previous sketch
	while ( iByte < binLocus.size() ) {
		uint64_t nWordUnsetBits{wordSizeInBits_};
		uint64_t nSumUnsetBits{0};
		while ( (nWordUnsetBits == wordSizeInBits_) && ( iByte < binLocus.size() ) ) {
			uint64_t locusChunk{allBitsSet_};
			assert( ( iByte < binLocus.size() ) // NOLINT
					&& "ERROR: iByte must be less than locus size in bytes in locusOPH_()" );
			const size_t nRemainingBytes{binLocus.size() - iByte};
			const size_t locusChunkSize{ (static_cast<size_t>(nRemainingBytes >= llWordSize_) * llWordSize_) + (static_cast<size_t>(nRemainingBytes < llWordSize_) * nRemainingBytes) };
			memcpy(&locusChunk, binLocus.data() + iByte, locusChunkSize);
			locusChunk    &= allBitsSet_ << sketchTail;
			nWordUnsetBits = _tzcnt_u64(locusChunk);
			nSumUnsetBits += nWordUnsetBits - sketchTail;
			sketchTail     = 0;
			iByte         += locusChunkSize;
		}
		iSketch += nSumUnsetBits / sketchSize_;
		if (iSketch >= kSketches_) {
			break;
		}
		filledIndexes.push_back(iSketch);
		sketches_[sketchBeg + iSketch] = static_cast<uint16_t>(nSumUnsetBits % sketchSize_);
		++iSketch;
		const uint64_t bitsDone{iSketch * sketchSize_};
		iByte      = bitsDone / byteSize_;
		sketchTail = bitsDone % byteSize_;
	}
	assert( (filledIndexes.size() <= kSketches_) // NOLINT
					&& "ERROR: filledIndexes.size() must not be greater than sketch number (kSketches_) in locusOPH_()" );
	densifySketches_(filledIndexes, sketchBeg);
}

void GenoTableHash::densifySketches_(std::vector<size_t> filledIndexes, const size_t &sketchBeg) {
	// The index progression is seeded identically for every locus, so an empty sketch resolves to the
	// same filled one across loci and the sketches stay comparable.
	RanDraw prng(emptyBinIdxSeed_);
	std::vector<uint32_t> seeds{static_cast<uint32_t>( prng.ranInt() )};
	if ( filledIndexes.empty() ) {                              // if the whole locus is monomorphic, pick a random index as filled
		filledIndexes.push_back( prng.sampleInt(kSketches_) );
	}
	size_t iSeed = 0;                                           // index into the seed vector
	size_t emptyCount = kSketches_ - filledIndexes.size();
	while (emptyCount > 0) {
		for (const auto eachFI : filledIndexes) {
			std::array<uint32_t, SIZE_OF_SIZET> key{};
			memcpy(key.data(), &eachFI, SIZE_OF_SIZET);
			auto newIdx = static_cast<uint32_t>( (murMurHash(key, seeds[iSeed]) % kSketches_) + sketchBeg );
			// should be safe: each thread accesses different vector elements
			if (sketches_[newIdx] == emptyBinToken_) {
				sketches_[newIdx] = sketches_[eachFI + sketchBeg];
				--emptyCount;
				break;
			}
		}
		++iSeed;
		if ( iSeed == seeds.size() ) {
			seeds.push_back( static_cast<uint32_t>( prng.ranInt() ) );
		}
	}
}

void GenoTableHash::bed2ophBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t>&bedLocusIndRange,
					const LocationWithLength &bedLocusSpan, const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv) {
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	size_t iLocus{bedLocusSpan.start};
	for (size_t iBedLocus = bedLocusIndRange.first; iBedLocus < bedLocusIndRange.second; ++iBedLocus) {
		std::vector<uint8_t> binLocus(locusSize_, 0);
		LocationWithLength bedWindow{0, 0};
		bedWindow.start  = iBedLocus * bedLocusSpan.length;
		bedWindow.length = bedLocusSpan.length;
		LocationWithLength binWindow{0, 0};
		binWindow.length = locusSize_;
		binarizeBedLocus(bedWindow, bedData, nIndividuals_, binWindow, binLocus, locusSeed_ + iLocus);
		// pad the locus to have a whole number of sketches 
		for (const auto &addI : padIndiv) {
			const size_t iLocByte    = addI.first / byteSize_;
			const auto iInLocByte    = static_cast<uint8_t>(addI.first % byteSize_);
			auto bytePair            = static_cast<uint16_t>(binLocus[iLocByte]);
			const size_t permByteInd = addI.second / byteSize_;
			const auto permInByteInd = static_cast<uint8_t>(addI.second % byteSize_);
			// Pair the current locus byte with the byte containing the value to be swapped
			// Then use the exchanging two fields trick from Hacker's Delight Chapter 2-20
			bytePair                |= static_cast<uint16_t>(binLocus[permByteInd] << byteSize_);
			const auto mask          = static_cast<uint16_t>(1 << iInLocByte);
			const auto shiftDistance = static_cast<uint16_t>( (byteSize_ - iInLocByte) + permInByteInd );                        // subtraction is safe b/c iInLocByte is modulo byteSize
			const auto temp1         = static_cast<uint16_t>( ( bytePair ^ (bytePair >> shiftDistance) ) & mask );
			const auto temp2         = static_cast<uint16_t>(temp1 << shiftDistance);
			bytePair                ^= temp1 ^ temp2;
			// Transfer bits using the trick in Hacker's Delight Chapter 2-20 (do not need the full swap, just transfer from the byte pair to binLocus1)
			// Must modify the current byte in each loop iteration because permutation indexes may fall into it
			binLocus[iLocByte] ^= static_cast<uint8_t>( ( binLocus[iLocByte] ^ static_cast<uint8_t>(bytePair) ) & static_cast<uint8_t>(mask) );
		}
		locusOPH_(iLocus, permutation, binLocus);
		++iLocus;
	}
}

size_t GenoTableHash::bed2ophThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges,
						const LocationWithLength &bedLocusSpan, const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv) {
	// The thread ranges are contiguous and start at 0 (see makeThreadRanges), so a range's
	// output locus index is simply bedLocusSpan.start + bedLocusIndRange.first. Each range reads
	// a disjoint slice of bedData and writes a disjoint set of OPH sketches, so the loop is
	// data-parallel; the ThreadCeiling caps concurrency to nThreads_ (a no-op without a TBB backend).
	const ThreadCeiling threadCeiling(nThreads_);
	std::for_each(
		parallelPolicy,
		threadRanges.cbegin(),
		threadRanges.cend(),
		[this, &bedData, &bedLocusSpan, &permutation, &padIndiv](const std::pair<size_t, size_t> &bedLocusIndRange) {
			LocationWithLength currentBedLocusSpan{0, 0};
			currentBedLocusSpan.start  = bedLocusSpan.start + bedLocusIndRange.first;
			currentBedLocusSpan.length = bedLocusSpan.length;
			bed2ophBlk_(bedData, bedLocusIndRange, currentBedLocusSpan, permutation, padIndiv);
		}
	);
	// contiguous-from-0 ranges sum to back().second,
	// so the next free locus index is bedLocusSpan.start + that total.
	return bedLocusSpan.start + threadRanges.back().second;
}

size_t GenoTableHash::bed2oph_(const BedDataStats &locusGroupStats, std::fstream &bedStream, const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv) {
	CountAndSize threadCounts{0, 0};
	threadCounts.count = nThreads_;
	threadCounts.size  = locusGroupStats.nLociPerThread;
	std::vector< std::pair<size_t, size_t> > threadRanges{makeThreadRanges(threadCounts)};
	assert( (locusGroupStats.nLociToRead >= threadRanges.back().second) // NOLINT
				&& "ERROR: nLociToRead smaller than threadRanges.back().second in bed2oph_()");
	assert( ( locusGroupStats.nBytesToRead < std::numeric_limits<std::streamsize>::max() ) // NOLINT
				&& "ERROR: amount to read larger than maximum streamsize in bed2oph_()");
	// Extend the last thread's range to cover the remainder loci; bed2ophThreaded_
	// then writes all nLociToRead loci of this chunk contiguously and returns the
	// next free index (see the matching fix in bed2bin_).
	threadRanges.back().second = locusGroupStats.nLociToRead;
	std::vector<char> bedChunkToRead(locusGroupStats.nBytesToRead, 0);
	size_t locusInd{locusGroupStats.firstLocusIdx};
	for (size_t iChunk = 0; iChunk < locusGroupStats.nMemChunks; ++iChunk) {
		bedStream.read( bedChunkToRead.data(), static_cast<std::streamsize>(locusGroupStats.nBytesToRead) );
		LocationWithLength bedLocusSpan{0, 0};
		bedLocusSpan.start  = locusInd;
		bedLocusSpan.length = locusGroupStats.nBytesPerLocus;
		locusInd            = bed2ophThreaded_(bedChunkToRead, threadRanges, bedLocusSpan, permutation, padIndiv);
	}
	return locusInd;
}

void GenoTableHash::mac2ophBlk_(const std::vector<int> &macData, const std::pair<size_t, size_t> &blockRange,
				const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv) {
	for (size_t iLocus = blockRange.first; iLocus < blockRange.second; ++iLocus) {
		std::vector<uint8_t> binLocus(locusSize_, 0);
		LocationWithLength binLocusRange{};
		binLocusRange.start  = 0;
		binLocusRange.length = locusSize_;
		std::vector<int> macLocus(nIndividuals_);
		const size_t nIndivUnpadded{nIndividuals_ - padIndiv.size()};
		std::copy_n(
			std::next( macData.cbegin(), static_cast<std::vector<int>::difference_type>(iLocus * nIndivUnpadded) ),
			nIndivUnpadded,
			macLocus.begin()
		);
		size_t iPad = nIndivUnpadded;
		std::for_each(
			padIndiv.cbegin(),
			padIndiv.cend(),
			[&iPad, &macLocus](const std::pair<size_t, size_t> &padIdxPair) {
				macLocus[iPad] = macLocus[padIdxPair.second];
				++iPad;
			}
		);
		binarizeMacLocus(macLocus, binLocusRange, binLocus, locusSeed_ + iLocus);
		locusOPH_(iLocus, permutation, binLocus);
	}
}

SimilarityMatrix GenoTableHash::hashJacBlock_(const std::pair<RowColIdx, RowColIdx> &blockRange, const std::vector<uint32_t> &locusIndexes, const float &similarityCutOff) const {
	SimilarityMatrix result;
	uint32_t iRow{blockRange.first.iRow};
	if (blockRange.first.iRow != blockRange.second.iRow) {
		for (uint32_t jCol = blockRange.first.jCol; jCol < iRow; ++jCol) { // first, possibly incomplete, row
			RowColIdx localRC{};
			localRC.iRow = locusIndexes[iRow];
			localRC.jCol = locusIndexes[jCol];
			const JaccardPair localJP{makeJaccardPair_(localRC)};
			if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
				result.insert(localRC, localJP);
			}
		}
		++iRow;
		while ( iRow < (blockRange.second.iRow) ) { // complete triangle
			for (uint32_t jCol = 0; jCol < iRow; ++jCol) {
				RowColIdx localRC{};
				localRC.iRow = locusIndexes[iRow];
				localRC.jCol = locusIndexes[jCol];
				const JaccardPair localJP{makeJaccardPair_(localRC)};
				if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
					result.insert(localRC, localJP);
				}
			}
			++iRow;
		}
		for (uint32_t jColRem = 0; jColRem < blockRange.second.jCol; ++jColRem) {  // last, possibly incomplete, row
			RowColIdx localRC{};
			localRC.iRow = locusIndexes[iRow];
			localRC.jCol = locusIndexes[jColRem];
			const JaccardPair localJP{makeJaccardPair_(localRC)};
			if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
				result.insert(localRC, localJP);
			}
		}
		return result;
	}
	for (uint32_t jCol = blockRange.first.jCol; jCol < blockRange.second.jCol; ++jCol) {  // last, possibly incomplete, row
		RowColIdx localRC{};
		localRC.iRow = locusIndexes[iRow];
		localRC.jCol = locusIndexes[jCol];
		const JaccardPair localJP{makeJaccardPair_(localRC)};
		if (static_cast<float>(localJP.nIntersect) / static_cast<float>(localJP.nUnion) >= similarityCutOff) {
			result.insert(localRC, localJP);
		}
	}

	return result;
}

SimilarityMatrix GenoTableHash::hashJacBlock_(const std::pair<HashGroupItPairCount, HashGroupItPairCount> &blockRange, const float &similarityCutOff) const {
	if (blockRange.first.hgIterator == blockRange.second.hgIterator) { // block falls entirely within a group
		std::pair<RowColIdx, RowColIdx> rowColumnPair{};
		rowColumnPair.first  = recoverRCindexes(blockRange.first.pairCount);
		rowColumnPair.second = recoverRCindexes(blockRange.second.pairCount);

		return hashJacBlock_(rowColumnPair, blockRange.first.hgIterator->locusIndexes, similarityCutOff);
	}
	std::pair<RowColIdx, RowColIdx> rowColumnPair{};
	rowColumnPair.first = recoverRCindexes(blockRange.first.pairCount);
	const size_t nPairs{blockRange.first.hgIterator->locusIndexes.size() * (blockRange.first.hgIterator->locusIndexes.size() - 1) / 2};
	rowColumnPair.second = recoverRCindexes(nPairs);
	SimilarityMatrix result{hashJacBlock_(rowColumnPair, blockRange.first.hgIterator->locusIndexes, similarityCutOff)};
	// Process complete groups, accumulating each group's (individually sorted) block by unordered
	// append rather than merging it in place. A per-group merge is a full-vector set_union rebuild,
	// so merging G groups one at a time is O(nElements * G); appending is O(1) amortized per group,
	// with a single sortAndDeduplicate() below restoring the sorted, de-duplicated invariant.
	std::for_each(
		blockRange.first.hgIterator + 1,
		blockRange.second.hgIterator,
		[this, &result, &similarityCutOff](const HashGroup &eachGroup) {
			std::pair<RowColIdx, RowColIdx> localRCPair{};
			localRCPair.first.iRow = 1;
			localRCPair.first.jCol = 0;
			const size_t locNpairs{eachGroup.locusIndexes.size() * (eachGroup.locusIndexes.size() - 1) / 2};
			localRCPair.second = recoverRCindexes(locNpairs);
			SimilarityMatrix tmp{hashJacBlock_(localRCPair, eachGroup.locusIndexes, similarityCutOff)};
			result.append(tmp);
		}
	);
	// last, possibly incomplete, group
	rowColumnPair.first.iRow = 1;
	rowColumnPair.first.jCol = 0;
	rowColumnPair.second     = recoverRCindexes(blockRange.second.pairCount);
	SimilarityMatrix tmp     = hashJacBlock_(rowColumnPair, blockRange.second.hgIterator->locusIndexes, similarityCutOff);
	result.append(tmp);
	result.sortAndDeduplicate();

	return result;
}

JaccardPair GenoTableHash::makeJaccardPair_(const RowColIdx &rowColumn) const noexcept {
	const auto kSkDst{static_cast<std::vector<uint16_t>::difference_type>(kSketches_)};
	const auto start1 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(rowColumn.iRow) * kSketches_;
	const auto start2 = sketches_.begin() + static_cast<std::vector<uint16_t>::difference_type>(rowColumn.jCol) * kSketches_;
	// count equal elements using the inner_product idiom
	const int32_t simVal = std::inner_product( start1, start1 + kSkDst, start2, 0, std::plus<>(), std::equal_to<>() );
	JaccardPair localJP{};
	// these values calculated from the hash comparison correspond to
	// intersection and union values that would have come from a
	// full binary vector comparison
	localJP.nIntersect = static_cast<uint32_t>(simVal);
	localJP.nUnion     = static_cast<uint32_t>(kSketches_);
	return localJP;
}
