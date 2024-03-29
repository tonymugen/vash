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
#include <fstream>

#include "random.hpp"
#include "gvarHash.hpp"

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
	uint16_t countSetBits(uint16_t inVal);
	/** \brief Count set bits in a vector
	 *
	 * Counting the set bits in a vector of bytes using Karnigan's method.
	 *
	 * \param[in] inVec input vector
	 * \return number of bits set
	 */
	uint64_t countSetBits(const std::vector<uint8_t> &inVec);
	/** \brief Count set bits in a range within a vector
	 *
	 * Counting the set bits in a range within a vector of bytes using Karnigan's method.
	 *
	 * \param[in] inVec input vector
	 * \param[in] window vector window in bytes
	 * \return number of bits set
	 */
	uint64_t countSetBits(const std::vector<uint8_t> &inVec, const LocationWithLength &window);
	/** \brief Get available RAM
	 *
	 * Estimates available RAM. If `procfs` is mounted, uses information from there. Otherwise, sets available RAM to 2 GiB.
	 *
	 * \return estimated available RAM in bytes
	 */
	size_t getAvailableRAM();
	/** \brief MurMurHash mixer module of an index value
	 *
	 * Generates a 32-bit an unfinalized hash of an index value using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	uint32_t murMurHashMixer(const std::array<uint32_t, SIZE_OF_SIZET> &key, const uint32_t &seed);
	/** \brief MurMurHash finalizer
	 *
	 * MurMurHash3 finalizer for a hash value.
	 *
	 * \param[in] inputHash input unfinlized hash value
	 * 
	 * \return finalized hash value
	 */
	uint32_t murMurHashFinalizer(const uint32_t &inputHash);
	/** \brief MurMurHash of an index value
	 *
	 * Generates a 32-bit hash of an index value using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	uint32_t murMurHash(const std::array<uint32_t, SIZE_OF_SIZET> &key, const uint32_t &seed);
	/** \brief MurMurHash of a vector of indexes
	 *
	 * Generates a 32-bit hash of an index value vector using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key vector to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	uint32_t murMurHash(const std::vector<size_t> &key, const uint32_t &seed);
	/** \brief MurMurHash of a vector of 32-bit unsigned integers
	 *
	 * Generates a 32-bit hash of a vector of unsigned 32-bit integers using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key vector to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	uint32_t murMurHash(const std::vector<uint32_t> &key, const uint32_t &seed);
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
	uint32_t murMurHash(const std::vector<uint16_t> &key, const LocationWithLength &keyWindow, const uint32_t &seed);
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
	std::vector< std::pair<size_t, size_t> > makeThreadRanges(const CountAndSize &threadPoolSizes);
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
	/** \brief Initialize an `IndexedPairSimilarity` vector 
	 *
	 * Creates a vector of `IndexedPairSimilarity` objects.
	 * Index pairs are set to reflect the possibility that this vector is a part of a larger vectorized lower triangle of a similarity matrix.
	 * All group IDs and similarity values are set to 0.
	 *
	 * \param[in] pairSpan the overall start and length (in elements) of the segment
	 * \param[in] totalNelements total number of elements
	 * \return vector of `IndexedPairSimilarity` objects
	 */
	std::vector<IndexedPairSimilarity> initializeIPSvector(const LocationWithLength &pairSpan, const size_t &totalNelements);
	/** \brief Create pair vector from groups 
	 *
	 * Create a vector of paired indexes within provided groups, in preparation for estimating Jaccard similarities.
	 *
	 * \param[in] firstGrpIdx index of the first group in the range
	 * \param[in] grpBlockStart iterator of the first group in the block
	 * \param[in] grpBlockEnd iterator of (one past the) end group in the block
	 * \return vector of pair indexes with group IDs
	 */
	std::vector<IndexedPairSimilarity> vectorizeGroups(const uint32_t &firstGrpIdx,
				const std::vector< std::vector<uint32_t> >::const_iterator grpBlockStart, const std::vector< std::vector<uint32_t> >::const_iterator grpBlockEnd);
	/** \brief Save values 
	 *
	 * Saves each value from the vector to the provided `fstream` with space as the delimiter.
	 *
	 * \param[in] inVec vector of floats to save
	 * \param[in, out] outputStream `fstream` to save to
	 */
	void saveValues(const std::vector<float> &inVec, std::fstream &outputStream);
	/** \brief Save indexed values 
	 *
	 * Saves each value with group IDs and the pair of indexes, tab delimited.
	 *
	 * \param[in] inVec vector of indexed floats to save
	 * \param[in, out] outputStream `fstream` to save to
	 */
	void saveValues(const std::vector<IndexedPairSimilarity> &inVec, std::fstream &outputStream);
	/** \brief Save named values 
	 *
	 * Saves each value with group IDs and names of locus pairs, tab delimited.
	 *
	 * \param[in] inVec vector of indexed floats to save
	 * \param[in] locusNames vector of locus names
	 * \param[in, out] outputStream `fstream` to save to
	 */
	void saveValues(const std::vector<IndexedPairSimilarity> &inVec, const std::vector<std::string> &locusNames, std::fstream &outputStream);
	/** \brief Save indexed LD values 
	 *
	 * Saves LD estimates (Jaccard and \f$r^2\f$) with locus pair indexes, tab delimited.
	 *
	 * \param[in] inVec vector of indexed LD values to save
	 * \param[in, out] outputStream `fstream` to save to
	 */
	void saveValues(const std::vector<IndexedPairLD> &inVec, std::fstream &outputStream);
	/** \brief Save named LD values 
	 *
	 * Saves LD estimates (Jaccard and \f$r^2\f$) with locus pair names, tab delimited.
	 *
	 * \param[in] inVec vector of indexed LD values to save
	 * \param[in] locusNames vector of locus names
	 * \param[in, out] outputStream `fstream` to save to
	 */
	void saveValues(const std::vector<IndexedPairLD> &inVec, const std::vector<std::string> &locusNames, std::fstream &outputStream);
	/** \brief Extract locus names 
	 *
	 * Extract locus names from a _.bim_ file.
	 *
	 * \param[in] bimFileName _.bim_ file name
	 * \return vector of locus names
	 */
	std::vector<std::string> getLocusNames(const std::string &bimFileName);
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
	 * \param[out] stringVariables indexed `std::string` variables for use by `main()`
	 */
	void extractCLinfo(const std::unordered_map<std::string, std::string> &parsedCLI, std::unordered_map<std::string, int> &intVariables, std::unordered_map<std::string, std::string> &stringVariables);
}
