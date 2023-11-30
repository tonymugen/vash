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

/// Binned OPH class
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023 Anthony J. Greenberg
 * \version 0.5
 *
 * Class and friend function definitions for the binned one-permutation hash.
 *
 */

#pragma once

#include <cstddef>
#include <vector>
#include <cstdint>

#include "gvarHash.hpp"

namespace BayesicSpace {
	class BinnedOPH;

	/** \brief Binned OPH class
	 *
	 * Encapsulates a binned one-permutation hash.
	 * Provides facilities to used the class as a key for `unordered_map`.
	 */
	class BinnedOPH {
	public:
		/** \brief Default constructor */
		BinnedOPH() = default;
		/** \brief Constructor form OPH
		 *
		 * Uses a one-permutation hash from a locus and generates bins of provided size.
		 * Each bin is hashed to `uint32_t` using `murMurHash`.
		 * The `murMurHash` seed is provided because the same bands must be hashed to the same values across loci.
		 *
		 * \param[in] ophVector the vector of OPH sketches
		 * \param[in] mmHashSeed `murMurHash` seed to hash bins
		 * \param[in] locusCoordinates index and length (in sketch number) of the locus to hash
		 * \param[in] nRowsPerBand number of rows per bin
		 */
		BinnedOPH(const std::vector<uint16_t> &ophVector, const uint32_t &mmHashSeed, const LocationWithLength &locusCoordinates, const size_t &nRowsPerBin);
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		BinnedOPH(const BinnedOPH &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `BinnedOPH` object
		 */
		BinnedOPH& operator=(const BinnedOPH &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		BinnedOPH( BinnedOPH &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `BinnedOPH` object
		 */
		BinnedOPH& operator=( BinnedOPH &&toMove) noexcept = default;
		/** \brief Destructor */
		~BinnedOPH() = default;
	private:
		std::vector<uint32_t> binnedHash_;
	};
}
