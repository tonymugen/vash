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
 * Class and friend function implementation for the binned one-permutation hash.
 *
 */

#include <vector>
#include <cassert>

#include "binnedOPH.hpp"
#include "vashFunctions.hpp"

using namespace BayesicSpace;

BinnedOPH::BinnedOPH(const std::vector<uint16_t> &ophVector, const uint32_t &mmHashSeed, const LocationWithLength &locusCoordinates, const size_t &nRowsPerBin) {
	binnedHash_.reserve(locusCoordinates.length / nRowsPerBin);
	size_t iSketch = 0;
	while (iSketch < locusCoordinates.length) {
		std::vector<uint16_t> bandVec;

		auto firstSketchIt = ophVector.cbegin()
			+ static_cast<std::vector<uint16_t>::difference_type>(iSketch + locusCoordinates.start * locusCoordinates.length);         // iSketch tracks band IDs
		auto lastSketchIt = firstSketchIt
			+ static_cast<std::vector<uint16_t>::difference_type>(nRowsPerBin);
		std::copy( firstSketchIt, lastSketchIt, std::back_inserter(bandVec) );

		LocationWithLength bandVecWindow{0, 0};
		bandVecWindow.start  = 0;
		bandVecWindow.length = bandVec.size();
		const uint32_t hash  = murMurHash(bandVec, bandVecWindow, mmHashSeed);
		binnedHash_.push_back(hash);

		iSketch += nRowsPerBin;
	}
}

