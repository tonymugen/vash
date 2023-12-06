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

/// Similarity matrix
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2023 Anthony J. Greenberg
 * \version 0.1
 *
 * Method implementation for a compact representation of a (possibly sparse) similarity matrix.
 *
 */

#include <vector>
#include <array>
#include <string>
#include <iterator>
#include <cstdint>
#include <algorithm>
#include <limits>

#include "similarityMatrix.hpp"

using namespace BayesicSpace;

constexpr std::array<float, 256> SimilarityMatrix::floatLookUp_{
	0.0000F, 0.0039F, 0.0078F, 0.0118F, 0.0157F, 0.0196F, 0.0235F, 0.0275F, 0.0314F, 0.0353F, 0.0392F, 
	0.0431F, 0.0471F, 0.0510F, 0.0549F, 0.0588F, 0.0627F, 0.0667F, 0.0706F, 0.0745F, 0.0784F, 0.0824F, 
	0.0863F, 0.0902F, 0.0941F, 0.0980F, 0.1020F, 0.1059F, 0.1098F, 0.1137F, 0.1176F, 0.1216F, 0.1255F,
	0.1294F, 0.1333F, 0.1373F, 0.1412F, 0.1451F, 0.1490F, 0.1529F, 0.1569F, 0.1608F, 0.1647F, 0.1686F,
	0.1725F, 0.1765F, 0.1804F, 0.1843F, 0.1882F, 0.1922F, 0.1961F, 0.2000F, 0.2039F, 0.2078F, 0.2118F,
	0.2157F, 0.2196F, 0.2235F, 0.2275F, 0.2314F, 0.2353F, 0.2392F, 0.2431F, 0.2471F, 0.2510F, 0.2549F,
	0.2588F, 0.2627F, 0.2667F, 0.2706F, 0.2745F, 0.2784F, 0.2824F, 0.2863F, 0.2902F, 0.2941F, 0.2980F,
	0.3020F, 0.3059F, 0.3098F, 0.3137F, 0.3176F, 0.3216F, 0.3255F, 0.3294F, 0.3333F, 0.3373F, 0.3412F,
	0.3451F, 0.3490F, 0.3529F, 0.3569F, 0.3608F, 0.3647F, 0.3686F, 0.3725F, 0.3765F, 0.3804F, 0.3843F,
	0.3882F, 0.3922F, 0.3961F, 0.4000F, 0.4039F, 0.4078F, 0.4118F, 0.4157F, 0.4196F, 0.4235F, 0.4275F,
	0.4314F, 0.4353F, 0.4392F, 0.4431F, 0.4471F, 0.4510F, 0.4549F, 0.4588F, 0.4627F, 0.4667F, 0.4706F,
	0.4745F, 0.4784F, 0.4824F, 0.4863F, 0.4902F, 0.4941F, 0.4980F, 0.5020F, 0.5059F, 0.5098F, 0.5137F,
	0.5176F, 0.5216F, 0.5255F, 0.5294F, 0.5333F, 0.5373F, 0.5412F, 0.5451F, 0.5490F, 0.5529F, 0.5569F,
	0.5608F, 0.5647F, 0.5686F, 0.5725F, 0.5765F, 0.5804F, 0.5843F, 0.5882F, 0.5922F, 0.5961F, 0.6000F,
	0.6039F, 0.6078F, 0.6118F, 0.6157F, 0.6196F, 0.6235F, 0.6275F, 0.6314F, 0.6353F, 0.6392F, 0.6431F,
	0.6471F, 0.6510F, 0.6549F, 0.6588F, 0.6627F, 0.6667F, 0.6706F, 0.6745F, 0.6784F, 0.6824F, 0.6863F,
	0.6902F, 0.6941F, 0.6980F, 0.7020F, 0.7059F, 0.7098F, 0.7137F, 0.7176F, 0.7216F, 0.7255F, 0.7294F,
	0.7333F, 0.7373F, 0.7412F, 0.7451F, 0.7490F, 0.7529F, 0.7569F, 0.7608F, 0.7647F, 0.7686F, 0.7725F,
	0.7765F, 0.7804F, 0.7843F, 0.7882F, 0.7922F, 0.7961F, 0.8000F, 0.8039F, 0.8078F, 0.8118F, 0.8157F,
	0.8196F, 0.8235F, 0.8275F, 0.8314F, 0.8353F, 0.8392F, 0.8431F, 0.8471F, 0.8510F, 0.8549F, 0.8588F,
	0.8627F, 0.8667F, 0.8706F, 0.8745F, 0.8784F, 0.8824F, 0.8863F, 0.8902F, 0.8941F, 0.8980F, 0.9020F,
	0.9059F, 0.9098F, 0.9137F, 0.9176F, 0.9216F, 0.9255F, 0.9294F, 0.9333F, 0.9373F, 0.9412F, 0.9451F,
	0.9490F, 0.9529F, 0.9569F, 0.9608F, 0.9647F, 0.9686F, 0.9725F, 0.9765F, 0.9804F, 0.9843F, 0.9882F,
	0.9922F, 0.9961F, 1.0000F
};

void SimilarityMatrix::append(const RowColIdx &rowColPair, uint8_t quantSimilarity) {
	if (rowColPair.iRow == rowColPair.jCol) {
		throw std::string("ERROR: row and column indexes must be different in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (rowColPair.iRow == 0) {
		throw std::string("ERROR: row index must be > 0 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	// make sure we are addressing the lower triangle
	const auto iRow{static_cast<uint64_t>( std::max(rowColPair.iRow, rowColPair.jCol) )};
	const auto jCol{static_cast<uint64_t>( std::min(rowColPair.iRow, rowColPair.jCol) )};
	const uint64_t newVecIndex = (iRow - 1UL) * iRow / 2 + jCol;

	if (newVecIndex <= lastCumulativeIndex_) {
		throw std::string("ERROR: new entry in front of the last element in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	uint64_t newDiff = newVecIndex - lastCumulativeIndex_;

	// if the difference overflows uint32_t, pad with 0-valued entries
	const uint64_t n32bits = newDiff / std::numeric_limits<uint32_t>::max();
	IndexedSimilarity padding{};
	padding.diffIndex  = std::numeric_limits<uint32_t>::max();
	padding.similarity = 0;
	std::vector<IndexedSimilarity> whole32bits(n32bits, padding);
	std::move( whole32bits.begin(), whole32bits.end(), std::back_inserter(matrix_) );
	newDiff = newDiff % std::numeric_limits<uint32_t>::max();

	IndexedSimilarity newEntry{};
	newEntry.similarity = quantSimilarity;
	newEntry.diffIndex  = static_cast<uint32_t>(newDiff);
	matrix_.emplace_back(newEntry);
	lastCumulativeIndex_ = newVecIndex;
}
