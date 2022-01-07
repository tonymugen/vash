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
 * \version 0.1
 *
 * Definitions and interface documentation for classes that take binary variant files and generate lossy summaries with hashing.
 *
 */

#ifndef gvhash_hpp
#define gvhash_hpp

#include <algorithm>
#include <cstddef>
#include <vector>
#include <array>
#include <string>

#include "bayesicUtilities/random.hpp"

using std::vector;
using std::array;
using std::string;

namespace BayesicSpace {
	class GenoTable;

	/** \brief Class to store compressed genotype tables
	 *
	 * Provides facilities to store and manipulate compressed genotype tables.
	 * Genotypes are stored in a one-bit format: bit set for the minor allele, unset for the major.
	 * Bits corresponding to missing data are unset (this is the same as mean imputation), heterozygotes are set with a 50% probability.
	 */
	class GenoTable {
	public:
		/** \brief Default constructor */
		GenoTable(){};
		/** \brief Constructor with input file name
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed).
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * If necessary, alleles are re-coded so that the set bit is always the minor allele.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] nIndividuals number of genotyped individuals
		 */
		GenoTable(const string &inputFileName, const size_t &nIndividuals);
		/** \brief Constructor with count vector
		 *
		 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * The counts are checked and re-coded if necessary so that set bits represent the minor allele. This function should run faster if the 0 is the major allele homozygote.
		 * While the above values are the norm, any negative number will be interpreted as missing, any odd number as 1, and any (non-0) even number as 2.
		 *
		 * \param[in] maCounts vector of minor allele numbers
		 * \param[in] nIndividuals number of genotyped individuals
		 */
		GenoTable(const vector<int> &maCounts, const size_t &nIndividuals);

		/** \brief Copy constructor (deleted) */
		GenoTable(const GenoTable &in) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTable operator=(const GenoTable &in) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] in object to move
		 */
		GenoTable(GenoTable &&in) noexcept;
		/** \brief Move assignment operator
		 *
		 * \param[in] in object to be moved
		 * \return `GenoTable` object
		 */
		GenoTable& operator=(GenoTable &&in) noexcept;

		/** \brief Save the binary genotype file
		 *
		 * Saves the binary approximate genotype data to a binary file.
		 *
		 * \param[in] outFileName output file name
		 */
		void saveGenoBinary(const string &outFileName) const;
		/** \brief Make OPH sketches of individuals
		 *
		 * Make one-permutation sketches of individuals (one sketch per locus).
		 * The permutation of bits is using the Fisher-Yates-Durstenfeld algorithm.
		 * Filling in empty bins using the Mai _et al._ (2020) algorithm.
		 * Running this function modifies the internally stored genotype table in place.
		 *
		 * \param[in] kSketches the number of sketches per locus
		 */
		void makeIndividualOPH(const size_t &kSketches);
		/** \brief All by all LD from hashes
		 *
		 * Calculates linkage disequilibrium among all loci using a modified OPH.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * Expected similarities (\f$p_i \times p_j\f$) are subtracted from OPH similarities.
		 * This function must be run after a sketch-generating function (e.g., `makeIndividualOPH()`). This is checked and an exception thrown if the sketch vector is empty.
		 *
		 * \return lower triangle of the LD matrix
		 */
		vector<float> allHashLD() const;
		/** \brief All by all Jaccad similarity LD
		 *
		 * Calculates linkage disequilibrium among all loci using a corrected Jaccard similarity as the statistic.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * Expected similarities (\f$p_i \times p_j\f$) are subtracted from Jaccard similarities.
		 *
		 * \return lower triangle of the LD matrix
		 */
		vector<float> allJaccardLD() const;
		/** \brief Assign groups from OPH portions
		 *
		 * Use sketch portions to assign loci to groups.
		 * The number of elements to use must be smaller than the number of sketches (or an exception is thrown).
		 * If element number is odd, the next smallest even number is used.
		 *
		 * \param[in] nElements number of OPH elements to use
		 * \return vector of group IDs for each locus
		 */
		vector<uint16_t> assignGroups(const size_t &nElements) const;
		/** \brief Assign groups from simHashing OPH
		 *
		 * Use a 16 bit simHash of the whole OPH to assign loci to groups.
		 *
		 * \return group IDs for each locus
		 */
		vector<uint16_t> assignGroups() const;
		/** \brief Group loci by linkage disequilibrium (LD)
		 *
		 * Group loci by LD along the genome. The algorithm is
		 * Start by using simHash on the first `kSketches` of the first locus OPH. Proceed along the genome, for each locus
		 *  - simHash `kSketches` of the OPH
		 *  - compare to the latest group simHash
		 *  - if the Hamming distance from the latest group is less than `hammingCutoff`, add the locus index to the group
		 *  - if not, compare to up to `lookBackNumber` of groups back along the genome, adding to the first group that meets the cut-off
		 *  - if none of the previous groups are close enough, start a new group, labeling it with the current simHash.
		 *
		 *  \param[in] hammingCutoff the maximum Hamming distance for group inclusion
		 *  \param[in] kSketches number of OPH sketches to use for simHash
		 *  \param[in] lookBackNumber number of previous groups to consider
		 *  \param[in] outFileName name of the output file
		 */
		void groupByLD(const uint16_t &hammingCutoff, const size_t &kSketches, const size_t &lookBackNumber, const string &outFileName) const;
	protected:
		/** \brief Binarized genotype table
		 *
		 * Stores one bit per genotype. Heterozygotes are randomly assigned, missing data are assigned 0.
		 */
		vector<uint8_t> binGenotypes_;
		/** \brief Vector of sketches
		 *
		 * A sketch is the position of the first set bit in a bin of permuted bits.
		 */
		vector<uint16_t> sketches_;
		/** \brief Alternative allele frequencies
		 *
		 * One value per locus. This is typically the minor allele frequency.
		 */
		vector<float> aaf_;
		/** \brief Number of individuals */
		size_t nIndividuals_;
		/** \brief Number of loci */
		size_t nLoci_;
		/** \brief Locus size in bytes */
		size_t locusSize_;
		/** \brief Random number generator */
		RanDraw rng_;
		/** \brief Leading bytes for .bed files */
		static const array<char, 3> magicBytes_;
		/** \brief One set bit for masking */
		static const uint8_t oneBit_;
		/** \brief Size of one byte in bits */
		static const uint8_t byteSize_;
		/** \brief 64 bit word size in bytes */
		static const uint8_t llWordSize_;
		/** \brief MurMurHash number of blocks */
		static const size_t nblocks_;
		/** \brief MurMurHash key length */
		static const uint32_t mmhKeyLen_;
		/** \brief Value corresponding to an empty token */
		static const uint16_t emptyBinToken_;
		/** \brief MurMurHash c1 constant */
		static const uint32_t c1_;
		/** \brief MurMurHash c2 constant */
		static const uint32_t c2_;
		/** \brief Generate binary genotypes
		 *
		 * Generate binary genotypes from the genotype table.
		 */
		void generateBinGeno_();
		/** \brief MurMurHash to fill in empty bins
		 *
		 * Generates a 32-bit hash of an index value using the MurMurHash3 algorithm.
		 *
		 * \param[in] key the key to be hashed
		 * \param[in] seed the seed
		 *
		 * \return the hash value
		 */
		uint32_t murMurHash_(const size_t &key, const uint32_t &seed) const;
		/** \brief 16 bit MurMurHash
		 *
		 * Generates a 16-bit hash of a 16 bit value using the MurMurHash3 algorithm.
		 *
		 * \param[in] key the key to be hashed
		 * \param[in] seed the seed
		 *
		 * \return the hash value
		 */
		uint16_t murMurHash_(const uint16_t &key, const uint32_t &seed) const;
		/** \brief MurMurHash of a sketch portion
		 *
		 * Generates a 32-bit hash of a portion of an OPH sketch using the MurMurHash3 algorithm.
		 * If the number of elements provided is odd, it is rounded down to the next even number.
		 *
		 * \param[in] startInd index of the first sketch element (to the `sketches_` vector)
		 * \param[in] nElements number of elements to hash
		 * \param[in] seed the seed
		 *
		 * \return the hash value
		 */
		uint32_t murMurHash_(const size_t &startInd, const size_t &nElements, const uint32_t &seed) const;
		/** \brief SimHash of OPH sketches
		 *
		 * Takes the collection of OPH sketches for each locus and returns a simHash.
		 *
		 * \param[in] startInd index of the first sketch in the `sketches_` vector
		 * \param[in] kSketches number of sketches to hash
		 * \param[in] seed seed value to use with murMurHash on each sketch element
		 * \return hash value
		 */
		uint16_t simHash_(const size_t &startInd, const size_t &kSketches, const uint32_t &seed) const;
		/** \brief Count set bits in a 16-bit word
		 *
		 * Counting the set bits using Karnigan's method. Passing by value to modify the copy and also because the address is much bigger than 16 bits.
		 *
		 * \param[in] inVal input value
		 * \return number of bits set
		 */
		uint16_t countSetBits_(uint16_t inVal) const;
		/** \brief Count set bits in a vector
		 *
		 * Counting the set bits in a vector of bytes using Karnigan's method.
		 *
		 * \param[in] inVec input vector
		 * \return number of bits set
		 */
		uint32_t countSetBits_(const vector<uint8_t> &inVec) const;
		/** \brief Count set bits in a range within a vector
		 *
		 * Counting the set bits in a range within a vector of bytes using Karnigan's method.
		 *
		 * \param[in] inVec input vector
		 * \param[in] start staring index
		 * \param[in] length number of bytes to process
		 * \return number of bits set
		 */
		uint32_t countSetBits_(const vector<uint8_t> &inVec, const size_t &start, const size_t &length) const;
		/** \brief Jaccard similarity in a block of loci
		 *
		 * \param[in] iLocus first locus index
		 * \param[in] blockInd index (in `jaccardVec`) of the first element in the block
		 * \param[out] jaccardVec vector of Jaccard similarities
		 */
		void jaccardBlock_(const size_t &iLocus, const size_t &blockInd, vector<float> &jaccardVec) const;
		/** \brief Hash-based similarity in a continuous block of loci
		 *
		 * \param[in] iLocus first locus index
		 * \param[in] blockInd index (in `jaccardVec`) of the first element in the block
		 * \param[in] kSketches number of sketches
		 * \param[in] invK 1/K (where K is the number of sketches)
		 * \param[out] hashJacVec vector of hash-estimated Jaccard similarities
		 */
		void hashJacBlock_(const size_t &iLocus, const size_t &blockInd, const size_t &kSketches, const float &invK, vector<float> &hashJacVec) const;
		/** \brief Hash-based similarity in a (possibly discontinuous) block of loci
		 *
		 * \param[in] iLocus first locus index
		 * \param[in] jLocus vector of second locus indexes
		 * \param[in] kSketches number of sketches
		 * \param[in] invK 1/K (where K is the number of sketches)
		 * \param[out] hashJacVec vector of hash-estimated Jaccard similarities
		 */
		void hashJacBlock_(const size_t &iLocus, const vector<size_t> &jLocus, const size_t &kSketches, const float &invK, vector<float> &hashJacVec) const;
		/** \brief Hamming distance
		 *
		 * Calculates the bit-wise Hamming distance between two 16-bit variables. Passing the variables by value since they are much smaller than addresses.
		 *
		 * \param[in] first first value
		 * \param[in] second second value
		 * \return Hamming distance
		 *
		 */
		uint16_t inline hammingDistance_(uint16_t first, uint16_t second) const {return countSetBits_(first ^ second); }
	};
}

#endif // gvhash_hpp
