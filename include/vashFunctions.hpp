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
#include <fstream>

namespace BayesicSpace {
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
	 * \param[in] start staring index
	 * \param[in] length number of bytes to process
	 * \return number of bits set
	 */
	uint64_t countSetBits(const std::vector<uint8_t> &inVec, const size_t &start, const size_t &length);
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
	uint32_t murMurHashMixer(const size_t &key, const uint32_t &seed);
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
	uint32_t murMurHash(const size_t &key, const uint32_t &seed);
	/** \brief MurMurHash of a vector of indexes
	 *
	 * Generates a 32-bit hash of an index value using the MurMurHash3 algorithm.
	 *
	 * \param[in] key the key vector to be hashed
	 * \param[in] seed the seed
	 *
	 * \return the hash value
	 */
	uint32_t murMurHash(const std::vector<size_t> &key, const uint32_t &seed);
	/** \brief MurMurHash of a vector of indexes
	 *
	 * Generates a 32-bit hash of a vector of `uint16_t` values using the MurMurHash3 algorithm.
	 *
	 * \param[in] start start index
	 * \param[in] length number of elements to hash
	 * \param[in] key the vector to be hashed
	 * \param[in] seed the hash seed
	 *
	 * \return hash value
	 */
	uint32_t murMurHash(const size_t &start, const size_t &length, const std::vector<uint16_t> &key, const uint32_t &seed);
	/** \brief Test .bed magic bytes
	 *
	 * Throws if one of the input bytes does not match the three magic values in `plink` .bed files.
	 *
	 * \param[in] bytesToTest the byte set to test
	 */
	void testBedMagicBytes(std::array<char, 3> &bytesToTest);
	/** \brief Build thread ranges
	 *
	 * Build index ranges to use within each thread.
	 *
	 * \param[in] nThreads number of threads
	 * \param[in] nElementsPerThread number of elements per thread
	 * \return vector of index ranges
	 */
	std::vector< std::pair<size_t, size_t> > makeThreadRanges(const size_t &nThreads, const size_t &nElementsPerThread);
	/** \brief Save values 
	 *
	 * Saves each value from the vector to the provided `fstream` with space as the delimiter.
	 *
	 * \param[in] inVec vector of floats to save
	 * \param[in, out] outputStream `fstream` to save to
	 */
	void saveValues(const std::vector<float> &inVec, std::fstream &outputStream);
}
