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
		GenoTable() {};
		/** \brief Constructor with input file name
		 *
		 * \param[in] inputFileName input file name
		 */
		GenoTable(const string &inputFileName);
		/** \brief Constructor with count vector
		 *
		 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
		 *
		 * \param[in] maCounts vector of minor allele numbers
		 * \param[in] nIndividuals number of genotyped individuals
		 */
		GenoTable(const vector<int8_t> &maCounts, const size_t &nIndividuals);

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
	protected:
		/**\brief Genotype table
		 *
		 * May store all or part of the genotype file, depending on memory availability.
		 */
		vector<uint8_t> genotypes_;
		/** \brief Number of individuals */
		size_t nIndividuals_;
		/** \brief Number of loci */
		size_t nLoci_;
		/** \brief Radom number generator */
		RanDraw rng_;
		/** \brief Seeded PRNG for hashes */
		RanDraw pRNG_;
		/** \brief Leading bytes for .bed files */
		static const array<uint8_t, 3> magicBytes_;
	};
}

#endif // gvhash_hpp
