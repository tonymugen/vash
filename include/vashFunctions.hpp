/*
 * Copyright (c) 2023 Anthony J. Greenberg
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

/// Auxiliary functions for variant hashing
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023 Anthony J. Greenberg
 * \version 0.5
 *
 * Definitions of class-external functions needed by hashing classes.
 *
 */

#pragma once

#include <cstdint>
#include <cstddef>
#include <vector>
#include <array>
#include <unordered_map>

#include "gvarHash.hpp"
#include "similarityMatrix.hpp"

namespace BayesicSpace {
	/** \brief Number of 32-bit values in `size_t` */
	constexpr size_t SIZE_OF_SIZET{sizeof(size_t) / sizeof(uint32_t)};
	/** \brief Number of test bytes in a _.bed_ file */
	constexpr size_t N_BED_TEST_BYTES{3};
	/** \brief Count set bits in a 16-bit word
	 *
	 * Counting the set bits using Karnigan's method. Passing by value to modify the copy and also because the address is much bigger than 16 bits.
	 *
	 * \param[in] inVal input value
	 * \return number of bits set
	 */
	[[gnu::warn_unused_result]] uint16_t countSetBits(uint16_t inVal);
	/** \brief Count set bits in a vector
	 *
	 * Counting the set bits in a vector of bytes using Karnigan's method.
	 *
	 * \param[in] inVec input vector
	 * \return number of bits set
	 */
	[[gnu::warn_unused_result]] uint64_t countSetBits(const std::vector<uint8_t> &inVec);
	/** \brief Count set bits in a range within a vector
	 *
	 * Counting the set bits in a range within a vector of bytes using Karnigan's method.
	 *
	 * \param[in] inVec input vector
	 * \param[in] window vector window in bytes
	 * \return number of bits set
	 */
	[[gnu::warn_unused_result]] uint64_t countSetBits(const std::vector<uint8_t> &inVec, const LocationWithLength &window);
	/** \brief Get available RAM
	 *
	 * Estimates available RAM. If `procfs` is mounted, uses information from there. Otherwise, sets available RAM to 2 GiB.
	 *
	 * \return estimated available RAM in bytes
	 */
	[[gnu::warn_unused_result]] size_t getAvailableRAM();
	/** \brief MurMurHash mixer module of an index value
	 *
	 * Generates a 32-bit an unfinalized hash of an index value using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	[[gnu::warn_unused_result]] uint32_t murMurHashMixer(const std::array<uint32_t, SIZE_OF_SIZET> &key, const uint32_t &seed);
	/** \brief MurMurHash finalizer
	 *
	 * MurMurHash3 finalizer for a hash value.
	 *
	 * \param[in] inputHash input unfinlized hash value
	 * 
	 * \return finalized hash value
	 */
	[[gnu::warn_unused_result]] uint32_t murMurHashFinalizer(const uint32_t &inputHash);
	/** \brief MurMurHash of an index value
	 *
	 * Generates a 32-bit hash of an index value using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	[[gnu::warn_unused_result]] uint32_t murMurHash(const std::array<uint32_t, SIZE_OF_SIZET> &key, const uint32_t &seed);
	/** \brief MurMurHash of a vector of indexes
	 *
	 * Generates a 32-bit hash of an index value vector using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key vector to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	[[gnu::warn_unused_result]] uint32_t murMurHash(const std::vector<size_t> &key, const uint32_t &seed);
	/** \brief MurMurHash of a vector of 32-bit unsigned integers
	 *
	 * Generates a 32-bit hash of a vector of unsigned 32-bit integers using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key vector to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	[[gnu::warn_unused_result]] uint32_t murMurHash(const std::vector<uint32_t> &key, const uint32_t &seed);
	/** \brief MurMurHash of a vector of indexes
	 *
	 * Generates a 32-bit hash of a vector of `uint16_t` values using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the vector to be hashed
	 * \param[in] keyWindow the range of elements in the key to hash
	 * \param[in] seed the hash seed
	 *
	 * \return hash value
	 */
	[[gnu::warn_unused_result]] uint32_t murMurHash(const std::vector<uint16_t> &key, const LocationWithLength &keyWindow, const uint32_t &seed);
	/** \brief Test .bed magic bytes
	 *
	 * Throws if one of the input bytes does not match the three magic values in `plink` .bed files.
	 *
	 * \param[in] bytesToTest the byte set to test
	 */
	void testBedMagicBytes(const std::array<char, N_BED_TEST_BYTES> &bytesToTest);
	/** \brief Build thread ranges
	 *
	 * Build index ranges to use within each thread.
	 *
	 * \param[in] threadPoolSizes number of threads and number of loci per thread
	 * \return vector of index ranges
	 */
	[[gnu::warn_unused_result]] std::vector< std::pair<size_t, size_t> > makeThreadRanges(const CountAndSize &threadPoolSizes);
	/** \brief Build chunk sizes 
	 *
	 * Build a vector of chunk sizes. If the number of elements is not evenly divisible by the
	 * number of chunks, the remainder is spread among the first number of elements modulo number of chunks.
	 *
	 * \param[in] nElements number of elements
	 * \param[in] nChunks number of chunks
	 * \return vector of chunk sizes
	 */
	[[gnu::warn_unused_result]] std::vector<size_t> makeChunkSizes(const size_t &nElements, const size_t &nChunks);
	/** \brief Build chunk ranges 
	 *
	 * Build ranges of row/column index pairs for each chunk for a given length of a vectorized similarity matrix.
	 *
	 * \param[in] startAndChunkSize full matrix index and chunk size
	 * \param[in] nChunks number of chunks
	 * \return vector of row/column pair ranges, an element per chunk
	 */
	[[gnu::warn_unused_result]] std::vector< std::pair<RowColIdx, RowColIdx> > makeChunkRanges(const LocationWithLength &startAndChunkSize, const size_t nChunks);
	/** \brief Delimit a chunk of indexes
	 *
	 * Identifies the start and end hash table buckets (groups of loci) and indexes within them
	 * for chunked pair-wise LD estimation.
	 *
	 * \param[in] groupVector vector of hash table buckets
	 * \param[in] startHGPC start group iterator with index into the group
	 * \param[in] chunkSize size of the chunk to be processed
	 * \return pair of hash group (bucket) iterators with indexes into the first and last group
	 */
	[[gnu::warn_unused_result]] std::pair<HashGroupItPairCount, HashGroupItPairCount>
		makeGroupRanges(const std::vector<HashGroup> &groupVector, const HashGroupItPairCount &startHGPC, const size_t &chunkSize);
	/** \brief Convert a locus from _.bed_ to binary format
	 *
	 * Convert the _.bed_ two-bit format to one-bit binary.
	 *
	 *  - 00 (homozygous first _.bim_ allele) becomes 0 if it is the major allele, 1 otherwise
	 *  - 11 (homozygous second _.bim_ allele) becomes 1 if it is the minor allele, 0 otherwise
	 *  - 01 (missing genotype) becomes 0
	 *  - 10 (heterozygous) becomes 1 with probability 0.5
	 *
	 * If the number of individuals is not divisible by eight, the last binary byte is padded with 0s.
	 *
	 * \param[in] bedLocusWindow _.bed_ locus window
	 * \param[in] bedLocus vector of _.bed_ format bytes
	 * \param[in] nIndividuals number of individuals
	 * \param[in] binLocusWindow binary locus window
	 * \param[out] binLocus vector of binary format bytes
	 */
	void binarizeBedLocus(const LocationWithLength &bedLocusWindow, const std::vector<char> &bedLocus, const size_t &nIndividuals,
													const LocationWithLength &binLocusWindow, std::vector<uint8_t> &binLocus);
	/** \brief Convert a locus from a vector of minor allele counts
	 *
	 * Convert minor allele counts to one-bit binary.
	 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
	 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
	 * The counts are checked and re-coded if necessary so that set bits represent the minor allele. This function should run faster if the 0 is the major allele homozygote.
	 * While the above values are the norm, any negative number will be interpreted as missing, any odd number as 1, and any (non-0) even number as 2.
	 *
	 * If the number of individuals is not divisible by eight, the last binary byte is padded with 0s.
	 *
	 * \param[in] macLocus vector of minor allele counts at a locus
	 * \param[in] binLocusWindow window covering the binary locus
	 * \param[out] binLocus vector of binary format bytes
	 */
	void binarizeMacLocus(const std::vector<int> &macLocus, const LocationWithLength &binLocusWindow, std::vector<uint8_t> &binLocus);
	/** \brief Extract locus names 
	 *
	 * Extract locus names from a _.bim_ file.
	 *
	 * \param[in] bimFileName _.bim_ file name
	 * \return vector of locus names
	 */
	[[gnu::warn_unused_result]] std::vector<std::string> getLocusNames(const std::string &bimFileName);
	/** \brief Command line parser
	 *
	 * Maps flags to values. Flags assumed to be of the form `--flag-name value`.
	 *
	 * \param[in] argc size of the `argv` array
	 * \param[in] argv command line input array
	 * \param[out] cli map of tags to values
	 */
	void parseCL(int &argc, char **argv, std::unordered_map<std::string, std::string> &cli);
	/** \brief Extract parameters from parsed command line interface flags
	 *
	 * Extracts needed variable values, indexed by `std::string` encoded variable names.
	 *
	 * \param[in] parsedCLI flag values parsed from the command line
	 * \param[out] intVariables indexed `int` variables for use by `main()`
	 * \param[out] floatVariables indexed `float` variables for use by `main()`
	 * \param[out] stringVariables indexed `std::string` variables for use by `main()`
	 */
	void extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI, std::unordered_map<std::string, int> &intVariables,
			std::unordered_map<std::string, float> &floatVariables, std::unordered_map<std::string, std::string> &stringVariables);
}
