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
	 * Provides facilities to store and manipulate compressed genotype tables. Genotypes are stored in a two-bit format as in plink .bed files.
	 */
	class GenoTable {
	public:
		/** \brief Default constructor */
		GenoTable(){};
		/** \brief Constructor with input file name
		 *
		 * The suggested number of sketches is modified so that the number of individuals per sketch is divisible by 8.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] kSketches suggested number of sketches per locus
		 */
		GenoTable(const string &inputFileName, const size_t &nIndividuals, const size_t &kSketches);
		/** \brief Constructor with count vector
		 *
		 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
		 * The suggested number of sketches is modified so that the number of individuals per sketch is divisible by 8.
		 *
		 * \param[in] maCounts vector of minor allele numbers
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] kSketches suggested number of sketches per locus
		 */
		GenoTable(const vector<int8_t> &maCounts, const size_t &nIndividuals, const size_t &kSketches);

		/** \brief Copy constructor (deleted) */
		GenoTable(const GenoTable &in) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTable operator=(const GenoTable &in) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] in object to move
		 */
		GenoTable(GenoTable &&in);
		/** \brief Move assignment operator
		 *
		 * \param[in] in object to be moved
		 * \return `GenoTable` object
		 */
		GenoTable& operator=(GenoTable &&in);

		/** \brief Save .bed genotype file
		 *
		 * Saves the raw genotype data in the _plink_ .bed format.
		 *
		 * \param[in] outFileName output file name
		 */
		void saveGenoBed(const string &outFileName) const;
		/** \brief Save the binary genotype file
		 *
		 * Saves the binary approximate genotype data to a binary file.
		 *
		 * \param[in] outFileName output file name
		 */
		void saveGenoBinary(const string &outFileName) const;
		/** \brief Output bits
		 *
		 * Output bits as groups of eight 1s and 0s for debugging.
		 *
		 * \param[in] binVec binary vector
		 * \param[out] bitString string with bits
		 *
		 */
		void outputBits(const vector<uint8_t> &binVec, string &bitString) const;
		/** \brief Output _i_-th non-zero element ID
		 *
		 * \param[in] i index
		 *
		 * \return index of the first non-zero bit in the bin
		 */
		size_t getSketchIdx(const size_t &i) const { return static_cast<size_t>(sketches_[i]);};
		/** \brief All by all linkage disequilibrium
		 *
		 * Calculate similarities among all loci using a modified OPH.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 *
		 * \param[out] LDmat lower triangle of the LD matrix
		 */
		void allSimilarity(vector<float> &LDmat) const;
	protected:
		/** \brief Genotype table
		 *
		 * May store all or part of the genotype file, depending on memory availability.
		 */
		vector<uint8_t> genotypes_;
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
		 * One value per locus.
		 */
		vector<float> aaf_;
		/** \brief Number of individuals */
		size_t nIndividuals_;
		/** \brief Number of loci */
		size_t nLoci_;
		/** \brief Locus size in bytes */
		size_t locusSize_;
		/** \brief Number of sketches per locus */
		uint16_t kSketches_;
		/** \brief Sketch size in bits */
		uint16_t sketchSize_;
		/** \brief Radom number generator */
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
		/** \brief Generate sketches
		 *
		 * Generate sketches from binary genotypes using modified one-permutation hash.
		 * The permutation of bits is using the Fisher-Yates-Durstenfeld algorithm.
		 * The vector of seeds is updated as new seeds are required and re-used for subsequent loci.
		 *
		 * \param[in] locusIdx locus index
		 * \param[in] ranInts permutation integer vector (must be common among loci)
		 * \param[in,out] seeds updateable vector of seeds
		 */
		void makeSketches_(const size_t &locusIdx, const vector<size_t> &ranInts, vector<uint32_t> &seeds);
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
	};
}

#endif // gvhash_hpp
