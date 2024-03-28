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
#include <utility>  // for std::pair
#include <iterator>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <fstream>
#include <future>

#include "similarityMatrix.hpp"
#include "vashFunctions.hpp"

using namespace BayesicSpace;

void BayesicSpace::chunkedAppend(std::vector<uint64_t> &source, std::vector<uint64_t> &target) {
	const auto sqrtNelementsToMove = static_cast<std::vector<uint64_t>::difference_type>( std::sqrt( static_cast<float>( source.size() ) ) );
	std::vector<uint64_t>::difference_type iBlock{0};
	while (iBlock < sqrtNelementsToMove) {
		const auto endChunkIt = std::next(source.begin(), sqrtNelementsToMove);
		std::copy( source.begin(), endChunkIt, std::back_inserter(target) );
		// erasing from the front seems to be reasonably fast
		// reversing the vector and erasing from the end is not worth it: ~52X slower
		source.erase(source.begin(), endChunkIt);
		++iBlock;
	}
	std::copy( source.cbegin(), source.cend(), std::back_inserter(target) );
	source.clear();
}

uint64_t BayesicSpace::recoverFullVIdx(std::vector<uint64_t>::const_iterator packedElementIt) noexcept {
	return static_cast<uint64_t>( (*packedElementIt) >> SimilarityMatrix::valueSize_ );
};

RowColIdx BayesicSpace::recoverRCindexes(const uint64_t &packedElement) noexcept {
	constexpr double tfiCoeff{8.0};
	RowColIdx result{};
	const uint64_t vectorizedIdx = packedElement >> SimilarityMatrix::valueSize_;

	const auto row = static_cast<uint64_t>((1.0 + sqrt(1.0 + tfiCoeff * static_cast<double>(vectorizedIdx))) / 2.0);
	result.jCol    = static_cast<uint32_t>(vectorizedIdx - row * (row - 1) / 2);
	result.iRow    = static_cast<uint32_t>(row);

	return result;
}

RowColIdx BayesicSpace::recoverRCindexes(std::vector<uint64_t>::const_iterator packedElementIt) noexcept {
	const uint64_t thisFullIdx{recoverFullVIdx(packedElementIt)};
	RowColIdx result{recoverRCindexes(thisFullIdx)};
	return result;
}

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

constexpr std::array<const char*, 256> SimilarityMatrix::stringLookUp_{
	"0.0000", "0.0039", "0.0078", "0.0118", "0.0157", "0.0196", "0.0235", "0.0275", "0.0314", "0.0353", "0.0392",
	"0.0431", "0.0471", "0.0510", "0.0549", "0.0588", "0.0627", "0.0667", "0.0706", "0.0745", "0.0784", "0.0824",
	"0.0863", "0.0902", "0.0941", "0.0980", "0.1020", "0.1059", "0.1098", "0.1137", "0.1176", "0.1216", "0.1255",
	"0.1294", "0.1333", "0.1373", "0.1412", "0.1451", "0.1490", "0.1529", "0.1569", "0.1608", "0.1647", "0.1686",
	"0.1725", "0.1765", "0.1804", "0.1843", "0.1882", "0.1922", "0.1961", "0.2000", "0.2039", "0.2078", "0.2118",
	"0.2157", "0.2196", "0.2235", "0.2275", "0.2314", "0.2353", "0.2392", "0.2431", "0.2471", "0.2510", "0.2549",
	"0.2588", "0.2627", "0.2667", "0.2706", "0.2745", "0.2784", "0.2824", "0.2863", "0.2902", "0.2941", "0.2980",
	"0.3020", "0.3059", "0.3098", "0.3137", "0.3176", "0.3216", "0.3255", "0.3294", "0.3333", "0.3373", "0.3412",
	"0.3451", "0.3490", "0.3529", "0.3569", "0.3608", "0.3647", "0.3686", "0.3725", "0.3765", "0.3804", "0.3843",
	"0.3882", "0.3922", "0.3961", "0.4000", "0.4039", "0.4078", "0.4118", "0.4157", "0.4196", "0.4235", "0.4275",
	"0.4314", "0.4353", "0.4392", "0.4431", "0.4471", "0.4510", "0.4549", "0.4588", "0.4627", "0.4667", "0.4706",
	"0.4745", "0.4784", "0.4824", "0.4863", "0.4902", "0.4941", "0.4980", "0.5020", "0.5059", "0.5098", "0.5137",
	"0.5176", "0.5216", "0.5255", "0.5294", "0.5333", "0.5373", "0.5412", "0.5451", "0.5490", "0.5529", "0.5569",
	"0.5608", "0.5647", "0.5686", "0.5725", "0.5765", "0.5804", "0.5843", "0.5882", "0.5922", "0.5961", "0.6000",
	"0.6039", "0.6078", "0.6118", "0.6157", "0.6196", "0.6235", "0.6275", "0.6314", "0.6353", "0.6392", "0.6431",
	"0.6471", "0.6510", "0.6549", "0.6588", "0.6627", "0.6667", "0.6706", "0.6745", "0.6784", "0.6824", "0.6863",
	"0.6902", "0.6941", "0.6980", "0.7020", "0.7059", "0.7098", "0.7137", "0.7176", "0.7216", "0.7255", "0.7294",
	"0.7333", "0.7373", "0.7412", "0.7451", "0.7490", "0.7529", "0.7569", "0.7608", "0.7647", "0.7686", "0.7725",
	"0.7765", "0.7804", "0.7843", "0.7882", "0.7922", "0.7961", "0.8000", "0.8039", "0.8078", "0.8118", "0.8157",
	"0.8196", "0.8235", "0.8275", "0.8314", "0.8353", "0.8392", "0.8431", "0.8471", "0.8510", "0.8549", "0.8588",
	"0.8627", "0.8667", "0.8706", "0.8745", "0.8784", "0.8824", "0.8863", "0.8902", "0.8941", "0.8980", "0.9020",
	"0.9059", "0.9098", "0.9137", "0.9176", "0.9216", "0.9255", "0.9294", "0.9333", "0.9373", "0.9412", "0.9451",
	"0.9490", "0.9529", "0.9569", "0.9608", "0.9647", "0.9686", "0.9725", "0.9765", "0.9804", "0.9843", "0.9882",
	"0.9922", "0.9961", "1.0000"
};

constexpr uint64_t SimilarityMatrix::maxIdxBitfield_{0x00FFFFFFFFFFFFFF};
constexpr uint32_t SimilarityMatrix::maxRowColValue_{759250125};
constexpr uint64_t SimilarityMatrix::valueMask_{0x00000000000000FF};
constexpr uint64_t SimilarityMatrix::valueSize_{8};
constexpr uint64_t SimilarityMatrix::maxValueIdx_{0x00000000000000FF};

void SimilarityMatrix::insert(const RowColIdx &rowColPair, const JaccardPair &jaccardCounts) {
	if (rowColPair.iRow == rowColPair.jCol) {
		throw std::string("ERROR: row and column indexes must be different in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (rowColPair.iRow == 0) {
		throw std::string("ERROR: row index must be non-zero in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (jaccardCounts.nIntersect > jaccardCounts.nUnion) {
		throw std::string("ERROR: intersection count cannot be larger than the union count in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	if (jaccardCounts.nUnion == 0) {
		throw std::string("ERROR: union count cannot be 0 in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	// make sure we are addressing the lower triangle
	const auto rowColOrdered = std::minmax(rowColPair.iRow, rowColPair.jCol);
	if (rowColOrdered.first > maxRowColValue_) {
		throw std::string("ERROR: row or column index exceeds maximal allowable value in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
	}
	const auto jCol{static_cast<uint64_t>(rowColOrdered.first)};
	const auto iRow{static_cast<uint64_t>(rowColOrdered.second)};
	const uint64_t newVecIndex = (iRow - 1UL) * iRow / 2UL + jCol;

	const auto quantSimilarity = static_cast<uint8_t>( (jaccardCounts.nIntersect * maxValueIdx_) / jaccardCounts.nUnion );
	FullIdxValue tmp{};
	tmp.fullIdx         = newVecIndex;
	tmp.quantSimilarity = quantSimilarity;
	this->insert_(tmp);
}

void SimilarityMatrix::merge(SimilarityMatrix &&toMerge) {
	// trivial cases
	if ( matrix_.empty() ) {
		*this = std::move(toMerge);
		return;
	}
	if ( toMerge.matrix_.empty() ) {
		return;
	}

	// We will copy things in sqrt(size of the two matrices) chunks to conserve memory
	const float chunkSizeFloat = std::sqrt( static_cast<float>( matrix_.size() + toMerge.matrix_.size() ) );
	const auto  chunkSizeDiff  = static_cast<std::vector<uint64_t>::difference_type>(chunkSizeFloat);

	std::vector<uint64_t> mergedMatrix;
	std::vector<uint64_t> unionBuffer;

	auto idxComparison = [](const uint64_t &packedIdx1, const uint64_t &packedIdx2) {
		return (packedIdx1 >> valueSize_) < (packedIdx2 >> valueSize_);
	};

	// Initialize the merge matrix with the first pair of chunks
	const auto thisChunk1size{std::min( chunkSizeDiff, std::distance( matrix_.cbegin(), matrix_.cend() ) )};
	const auto thisChunk1endIt{matrix_.cbegin() + thisChunk1size};
	const auto toMergeChunk1size{std::min( chunkSizeDiff, std::distance( toMerge.matrix_.cbegin(), toMerge.matrix_.cend() ) )};
	const auto toMergeChunk1endIt{toMerge.matrix_.cbegin() + toMergeChunk1size};
	std::set_union(
		matrix_.cbegin(), thisChunk1endIt,
		toMerge.matrix_.cbegin(), toMergeChunk1endIt,
		std::back_inserter(mergedMatrix),
		idxComparison
	);
	matrix_.erase(matrix_.begin(), thisChunk1endIt);
	toMerge.matrix_.erase(toMerge.matrix_.begin(), toMergeChunk1endIt);
	while (std::max( matrix_.size(), toMerge.matrix_.size() ) > 0) {
		const auto thisChunkSize{
			std::min( chunkSizeDiff, std::distance( matrix_.cbegin(), matrix_.cend() ) )
		};
		const auto thisChunkEndIt{matrix_.cbegin() + thisChunkSize};
		std::set_union(
			mergedMatrix.cbegin(), mergedMatrix.cend(),
			matrix_.cbegin(), thisChunkEndIt,
			std::back_inserter(unionBuffer),
			idxComparison
		);
		mergedMatrix.clear();
		matrix_.erase(matrix_.begin(), thisChunkEndIt);
		const auto toMergeChunkSize{
			std::min( chunkSizeDiff, std::distance( toMerge.matrix_.cbegin(), toMerge.matrix_.cend() ) )
		};
		const auto toMergeChunkEndIt{toMerge.matrix_.cbegin() + toMergeChunkSize};
		std::set_union(
			unionBuffer.cbegin(), unionBuffer.cend(),
			toMerge.matrix_.cbegin(), toMergeChunkEndIt,
			std::back_inserter(mergedMatrix),
			idxComparison
		);
		unionBuffer.clear();
		toMerge.matrix_.erase(toMerge.matrix_.begin(), thisChunkEndIt);
	}

	assert( ( !matrix_.empty() && !toMerge.matrix_.empty() )  //NOLINT
			&& "ERROR: at least one matrix must be empty in merge" );
	// move over any straggling elements
	if ( !matrix_.empty() ) {
		chunkedAppend(matrix_, mergedMatrix);
	}
	if ( !toMerge.matrix_.empty() ) {
		chunkedAppend(toMerge.matrix_, mergedMatrix);
	}
	std::swap(mergedMatrix, matrix_);
}

void SimilarityMatrix::save(const std::string &outFileName, const size_t &nThreads, const std::string &locusNameFile) const {
	if ( matrix_.empty() ) {
		return;
	}
	// determine how many line strings can fit in RAM
	// using the last entry because it will lave longest indexes
	const RowColIdx lastIdxPair{recoverRCindexes( matrix_.back() )};

	std::vector<std::string> locusNames;
	if ( !locusNameFile.empty() ) {
		locusNames = getLocusNames(locusNameFile);
		// testing only the row index, since the column index must be smaller
		if ( lastIdxPair.iRow >= locusNames.size() ) {
			throw std::string("ERROR: number of rows exceeds locus name count in the ") + locusNameFile + std::string(" file in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
		}
	}

	const std::string lastEntry{std::to_string(lastIdxPair.iRow + 1) + "\t" + std::to_string(lastIdxPair.jCol + 1) + "\t" + stringLookUp_[0] + "\n"};
	const size_t maxInRAM = getAvailableRAM() / ( 2UL * lastEntry.size() );      // use half to leave resources for other operations
	const size_t nChunks  = std::max(matrix_.size() / maxInRAM, 1UL);
	std::vector<size_t> chunkSizes{makeChunkSizes(matrix_.size(), nChunks)};

	std::fstream outStream;
	outStream.open(outFileName, std::ios::out | std::ios::binary | std::ios::app);
	std::vector<uint32_t>::difference_type cumChunkSize{0};
	for (const auto &eachChunkSize : chunkSizes) {
		const size_t actualNthreads = std::min(eachChunkSize, nThreads);
		std::vector<size_t> threadChunkSizes{makeChunkSizes(eachChunkSize, actualNthreads)};
		std::vector< std::pair<std::vector<uint64_t>::const_iterator, std::vector<uint64_t>::const_iterator> > threadPairs;
		std::vector<uint32_t>::difference_type cumThreadChunkSize{cumChunkSize};
		std::for_each(
			threadChunkSizes.cbegin(),
			threadChunkSizes.cend(),
			[&threadPairs, &cumThreadChunkSize, this](const size_t &chunkSize) {
				std::pair<std::vector<uint64_t>::const_iterator, std::vector<uint64_t>::const_iterator> tmpPair{
					matrix_.cbegin() + cumThreadChunkSize,
					matrix_.cbegin() + cumThreadChunkSize + static_cast<std::vector<uint32_t>::difference_type>(chunkSize)
				};
				threadPairs.emplace_back(tmpPair);
				cumThreadChunkSize += static_cast<std::vector<uint32_t>::difference_type>(chunkSize);
			}
		);
		std::vector<std::string> outStrings(actualNthreads);
		std::vector< std::future<void> > tasks;
		tasks.reserve(actualNthreads);
		size_t iThread{0};
		std::for_each(
			threadPairs.cbegin(),
			threadPairs.cend(),
			[&iThread, &tasks, &outStrings, &locusNames, this](auto pairIt) {
				tasks.emplace_back(
					std::async(
						[iThread, pairIt, &outStrings, &locusNames, this]{
							outStrings.at(iThread) = stringify_(pairIt.first, pairIt.second, locusNames);
						}
					)
				);
				++iThread;
			}
		);
		for (const auto &eachThread : tasks) {
			eachThread.wait();
		}

		for (const auto &eachString : outStrings) {
			outStream.write( eachString.c_str(), static_cast<std::streamsize>( eachString.size() ) );
		}
		cumChunkSize += static_cast<std::vector<uint32_t>::difference_type>(eachChunkSize);
	}
	outStream.close();
}

std::string SimilarityMatrix::stringify_(std::vector<uint64_t>::const_iterator start, std::vector<uint64_t>::const_iterator end,
								const std::vector<std::string> &locusNames) {
	std::string outString;
	if ( locusNames.empty() ) {
		for (auto matIt = start; matIt != end; ++matIt) {
			const RowColIdx currentPair{recoverRCindexes(matIt)};
			const std::string similarityValue = stringLookUp_.at( static_cast<size_t>( (*matIt) & valueMask_ ) );
			outString += std::to_string(currentPair.iRow + 1) + "\t" + std::to_string(currentPair.jCol + 1) + "\t" + similarityValue + "\n";
		}
		return outString;
	}
	for (auto matIt = start; matIt != end; ++matIt) {
		const RowColIdx currentPair{recoverRCindexes(matIt)};
		const std::string similarityValue = stringLookUp_.at( static_cast<size_t>( (*matIt) & valueMask_ ) );
		outString += locusNames[currentPair.iRow] + "\t" + locusNames[currentPair.jCol] + "\t" + similarityValue + "\n";
	}
	return outString;
}

void SimilarityMatrix::insert_(const FullIdxValue &indexWithSimilarity) {
	const uint64_t packedElement{(indexWithSimilarity.fullIdx) << valueSize_ | static_cast<uint64_t>(indexWithSimilarity.quantSimilarity)};
	if ( matrix_.empty() ) {
		matrix_.push_back(packedElement);
		return;
	}
	// redundant element at back
	if ( (matrix_.back() >> valueSize_) == indexWithSimilarity.fullIdx ) {
		return;
	}
	// the additional element is past the current end
	if ( (matrix_.back() >> valueSize_) < indexWithSimilarity.fullIdx ) {
		matrix_.push_back(packedElement);
		return;
	}
	// will be inserting the new element
	if (fullIndexWithSimilarity.fullIdx < lastCumulativeIndex_) {
		const uint64_t newFirstCumIdx{std::min(firstCumulativeIndex_, fullIndexWithSimilarity.fullIdx)};
		uint64_t firstDiff{firstCumulativeIndex_ - newFirstCumIdx}; // 0 or the distance to newVecIndex if it is in front of firstCumulativeIndex_

		// resolution of the possible index bit-field overload by padding with 0 values
		const uint64_t nBFmax = firstDiff / maxIdxBitfield_;
		std::vector<uint32_t> wholeBitField(nBFmax, padding_);
		firstDiff = firstDiff % maxIdxBitfield_;

		uint32_t firstElement = static_cast<uint32_t>(firstDiff) << valueSize_;
		const auto firstValue{matrix_[0] & valueMask_};
		firstElement |= firstValue;
		matrix_[0]    = firstElement;   // if the firstCumulativeIndex_ has not changed, still 0 index
		auto matIt    = matrix_.begin();
		uint64_t currentIdx{recoverFullVIdx(matIt, newFirstCumIdx)};
		while (currentIdx < fullIndexWithSimilarity.fullIdx) {
			++matIt;
			currentIdx = recoverFullVIdx(matIt, currentIdx);
		}
		// ignore if the element already present
		if (currentIdx == fullIndexWithSimilarity.fullIdx) {
			matrix_[0] &= valueMask_;
			return;
		}

		auto nextDiff = static_cast<uint32_t>(currentIdx - fullIndexWithSimilarity.fullIdx);
		auto currDiff = ( (*matIt) >> valueSize_ ) - nextDiff;
		currDiff      = currDiff << valueSize_;
		currDiff     |= static_cast<uint32_t>(fullIndexWithSimilarity.quantSimilarity);
		matIt = std::next( matrix_.insert(matIt, currDiff) );
		const uint32_t nextVal{(*matIt) & valueMask_};
		nextDiff  = nextDiff << valueSize_;
		nextDiff |= nextVal;
		*matIt    = nextDiff;
		matrix_.insert( matIt, wholeBitField.cbegin(), wholeBitField.cend() );
		matrix_[0]           &= valueMask_;
		firstCumulativeIndex_ = newFirstCumIdx;
		return;
	}
	uint64_t newDiff = fullIndexWithSimilarity.fullIdx - lastCumulativeIndex_;

	// if the difference overflows the bit-field, pad with 0-valued entries
	const size_t nBFmax = newDiff / maxIdxBitfield_;
	std::vector<uint32_t> wholeBitField(nBFmax, padding_);
	std::move( wholeBitField.begin(), wholeBitField.end(), std::back_inserter(matrix_) );
	newDiff = newDiff % maxIdxBitfield_;

	uint32_t newEntry{static_cast<uint32_t>(newDiff) << valueSize_};
	newEntry |= static_cast<uint32_t>(fullIndexWithSimilarity.quantSimilarity);
	matrix_.push_back(newEntry);
	lastCumulativeIndex_ = fullIndexWithSimilarity.fullIdx;
}
