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

#include <numeric>
#include <vector>
#include <array>
#include <string>
#include <utility>  // for std::pair
#include <iterator>
#include <cstdint>
#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <future>
#include <thread>

#include "similarityMatrix.hpp"
#include "vashFunctions.hpp"

using namespace BayesicSpace;

void BayesicSpace::chunkedAppend(std::vector<uint32_t> &source, std::vector<uint32_t> &target) {
	const auto sqrtNelementsToMove = static_cast<std::vector<uint32_t>::difference_type>( std::sqrt( static_cast<float>( source.size() ) ) );
	std::vector<uint32_t>::difference_type iBlock{0};
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

uint64_t BayesicSpace::recoverFullVIdx(std::vector<uint32_t>::const_iterator packedElementIt, uint64_t precedingFullIdx) noexcept {
	return precedingFullIdx + static_cast<uint64_t>( (*packedElementIt) >> SimilarityMatrix::valueSize_ );
};

RowColIdx BayesicSpace::recoverRCindexes(const uint64_t &fullIdx) noexcept {
	constexpr double tfiCoeff{8.0};
	RowColIdx result{};

	const auto row = static_cast<uint64_t>((1.0 + sqrt(1.0 + tfiCoeff * static_cast<double>(fullIdx))) / 2.0);
	result.jCol    = static_cast<uint32_t>(fullIdx - row * (row - 1) / 2);
	result.iRow    = static_cast<uint32_t>(row);

	return result;
}

RowColIdx BayesicSpace::recoverRCindexes(std::vector<uint32_t>::const_iterator packedElementIt, uint64_t &precedingFullIdx) noexcept {
	const uint64_t thisFullIdx{recoverFullVIdx(packedElementIt, precedingFullIdx)};

	RowColIdx result{recoverRCindexes(thisFullIdx)};

	precedingFullIdx = thisFullIdx;

	return result;
}

constexpr std::array<float, 255> SimilarityMatrix::floatLookUp_{
	0.0000F, 0.0040F, 0.0079F, 0.0119F, 0.0158F, 0.0198F, 0.0237F, 0.0277F, 0.0316F, 0.0356F, 0.0395F,
	0.0435F, 0.0474F, 0.0514F, 0.0553F, 0.0593F, 0.0632F, 0.0672F, 0.0711F, 0.0751F, 0.0790F, 0.0830F,
	0.0869F, 0.0909F, 0.0948F, 0.0988F, 0.1027F, 0.1067F, 0.1106F, 0.1146F, 0.1185F, 0.1225F, 0.1264F,
	0.1304F, 0.1343F, 0.1383F, 0.1422F, 0.1462F, 0.1501F, 0.1541F, 0.1580F, 0.1620F, 0.1659F, 0.1699F,
	0.1738F, 0.1778F, 0.1817F, 0.1857F, 0.1896F, 0.1936F, 0.1975F, 0.2015F, 0.2054F, 0.2094F, 0.2133F,
	0.2173F, 0.2212F, 0.2252F, 0.2291F, 0.2331F, 0.2370F, 0.2410F, 0.2449F, 0.2489F, 0.2528F, 0.2568F,
	0.2607F, 0.2647F, 0.2686F, 0.2726F, 0.2765F, 0.2805F, 0.2844F, 0.2884F, 0.2923F, 0.2962F, 0.3002F,
	0.3042F, 0.3081F, 0.3121F, 0.3160F, 0.3200F, 0.3239F, 0.3279F, 0.3318F, 0.3358F, 0.3397F, 0.3436F,
	0.3476F, 0.3516F, 0.3555F, 0.3595F, 0.3634F, 0.3674F, 0.3713F, 0.3753F, 0.3792F, 0.3832F, 0.3871F,
	0.3911F, 0.3950F, 0.3990F, 0.4029F, 0.4069F, 0.4108F, 0.4148F, 0.4187F, 0.4227F, 0.4266F, 0.4306F,
	0.4345F, 0.4385F, 0.4424F, 0.4464F, 0.4503F, 0.4543F, 0.4582F, 0.4622F, 0.4661F, 0.4701F, 0.4740F,
	0.4780F, 0.4819F, 0.4859F, 0.4898F, 0.4938F, 0.4977F, 0.5016F, 0.5056F, 0.5096F, 0.5135F, 0.5175F,
	0.5214F, 0.5254F, 0.5293F, 0.5332F, 0.5372F, 0.5412F, 0.5451F, 0.5490F, 0.5530F, 0.5570F, 0.5609F,
	0.5649F, 0.5688F, 0.5728F, 0.5767F, 0.5807F, 0.5846F, 0.5886F, 0.5925F, 0.5964F, 0.6004F, 0.6044F,
	0.6083F, 0.6123F, 0.6162F, 0.6202F, 0.6241F, 0.6281F, 0.6320F, 0.6360F, 0.6399F, 0.6438F, 0.6478F,
	0.6518F, 0.6557F, 0.6597F, 0.6636F, 0.6676F, 0.6715F, 0.6755F, 0.6794F, 0.6834F, 0.6873F, 0.6912F,
	0.6952F, 0.6992F, 0.7031F, 0.7071F, 0.7110F, 0.7150F, 0.7189F, 0.7229F, 0.7268F, 0.7308F, 0.7347F,
	0.7386F, 0.7426F, 0.7466F, 0.7505F, 0.7545F, 0.7584F, 0.7624F, 0.7663F, 0.7703F, 0.7742F, 0.7782F,
	0.7821F, 0.7860F, 0.7900F, 0.7940F, 0.7979F, 0.8019F, 0.8058F, 0.8098F, 0.8137F, 0.8177F, 0.8216F,
	0.8256F, 0.8295F, 0.8335F, 0.8374F, 0.8414F, 0.8453F, 0.8493F, 0.8532F, 0.8572F, 0.8611F, 0.8651F,
	0.8690F, 0.8730F, 0.8769F, 0.8809F, 0.8848F, 0.8888F, 0.8927F, 0.8967F, 0.9006F, 0.9046F, 0.9085F,
	0.9125F, 0.9164F, 0.9204F, 0.9243F, 0.9283F, 0.9322F, 0.9362F, 0.9401F, 0.9441F, 0.9480F, 0.9520F,
	0.9559F, 0.9599F, 0.9638F, 0.9678F, 0.9717F, 0.9757F, 0.9796F, 0.9836F, 0.9875F, 0.9915F, 0.9954F,
	0.9994F, 1.0000F
};

constexpr std::array<const char*, 255> SimilarityMatrix::stringLookUp_{
	"0.0000", "0.0040", "0.0079", "0.0119", "0.0158", "0.0198", "0.0237", "0.0277", "0.0316", "0.0356", "0.0395",
	"0.0435", "0.0474", "0.0514", "0.0553", "0.0593", "0.0632", "0.0672", "0.0711", "0.0751", "0.0790", "0.0830",
	"0.0869", "0.0909", "0.0948", "0.0988", "0.1027", "0.1067", "0.1106", "0.1146", "0.1185", "0.1225", "0.1264",
	"0.1304", "0.1343", "0.1383", "0.1422", "0.1462", "0.1501", "0.1541", "0.1580", "0.1620", "0.1659", "0.1699",
	"0.1738", "0.1778", "0.1817", "0.1857", "0.1896", "0.1936", "0.1975", "0.2015", "0.2054", "0.2094", "0.2133",
	"0.2173", "0.2212", "0.2252", "0.2291", "0.2331", "0.2370", "0.2410", "0.2449", "0.2489", "0.2528", "0.2568",
	"0.2607", "0.2647", "0.2686", "0.2726", "0.2765", "0.2805", "0.2844", "0.2884", "0.2923", "0.2962", "0.3002",
	"0.3042", "0.3081", "0.3121", "0.3160", "0.3200", "0.3239", "0.3279", "0.3318", "0.3358", "0.3397", "0.3436",
	"0.3476", "0.3516", "0.3555", "0.3595", "0.3634", "0.3674", "0.3713", "0.3753", "0.3792", "0.3832", "0.3871",
	"0.3911", "0.3950", "0.3990", "0.4029", "0.4069", "0.4108", "0.4148", "0.4187", "0.4227", "0.4266", "0.4306",
	"0.4345", "0.4385", "0.4424", "0.4464", "0.4503", "0.4543", "0.4582", "0.4622", "0.4661", "0.4701", "0.4740",
	"0.4780", "0.4819", "0.4859", "0.4898", "0.4938", "0.4977", "0.5016", "0.5056", "0.5096", "0.5135", "0.5175",
	"0.5214", "0.5254", "0.5293", "0.5332", "0.5372", "0.5412", "0.5451", "0.5490", "0.5530", "0.5570", "0.5609",
	"0.5649", "0.5688", "0.5728", "0.5767", "0.5807", "0.5846", "0.5886", "0.5925", "0.5964", "0.6004", "0.6044",
	"0.6083", "0.6123", "0.6162", "0.6202", "0.6241", "0.6281", "0.6320", "0.6360", "0.6399", "0.6438", "0.6478",
	"0.6518", "0.6557", "0.6597", "0.6636", "0.6676", "0.6715", "0.6755", "0.6794", "0.6834", "0.6873", "0.6912",
	"0.6952", "0.6992", "0.7031", "0.7071", "0.7110", "0.7150", "0.7189", "0.7229", "0.7268", "0.7308", "0.7347",
	"0.7386", "0.7426", "0.7466", "0.7505", "0.7545", "0.7584", "0.7624", "0.7663", "0.7703", "0.7742", "0.7782",
	"0.7821", "0.7860", "0.7900", "0.7940", "0.7979", "0.8019", "0.8058", "0.8098", "0.8137", "0.8177", "0.8216",
	"0.8256", "0.8295", "0.8335", "0.8374", "0.8414", "0.8453", "0.8493", "0.8532", "0.8572", "0.8611", "0.8651",
	"0.8690", "0.8730", "0.8769", "0.8809", "0.8848", "0.8888", "0.8927", "0.8967", "0.9006", "0.9046", "0.9085",
	"0.9125", "0.9164", "0.9204", "0.9243", "0.9283", "0.9322", "0.9362", "0.9401", "0.9441", "0.9480", "0.9520",
	"0.9559", "0.9599", "0.9638", "0.9678", "0.9717", "0.9757", "0.9796", "0.9836", "0.9875", "0.9915", "0.9954",
	"0.9994", "1.0000"
};

constexpr uint32_t SimilarityMatrix::maxIdxBitfield_{0x00FFFFFF};
constexpr uint32_t SimilarityMatrix::valueMask_{0x000000FF};
constexpr uint32_t SimilarityMatrix::padding_{0xFFFFFFFF};
constexpr uint32_t SimilarityMatrix::valueSize_{8};
constexpr uint32_t SimilarityMatrix::maxValueIdx_{0x000000FE};

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
	const auto jCol{static_cast<uint64_t>(rowColOrdered.first)};
	const auto iRow{static_cast<uint64_t>(rowColOrdered.second)};
	const uint64_t newVecIndex = (iRow - 1UL) * iRow / 2 + jCol;

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
	if (toMerge.matrix_.size() == 1) {
		FullIdxValue tmp{};
		tmp.fullIdx         = toMerge.firstCumulativeIndex_;
		tmp.quantSimilarity = static_cast<uint8_t>( toMerge.matrix_.at(0) );
		this->insert_(tmp);
		return;
	}
	// figure out which matrix goes in front
	if (firstCumulativeIndex_ > toMerge.firstCumulativeIndex_) {
		std::swap(matrix_, toMerge.matrix_);
		std::swap(firstCumulativeIndex_, toMerge.firstCumulativeIndex_);
		std::swap(lastCumulativeIndex_, toMerge.lastCumulativeIndex_);
	}
	if (lastCumulativeIndex_ <= toMerge.firstCumulativeIndex_) { // append in this case
		// fix the differential value for the first element that will be moved
		// it used to count off the previous element in toMerge.matrix_ (if any),
		// now needs to be off the last element in this->matrix_
		// the rest will be automatically correct
		uint64_t firstMoveDiff{toMerge.firstCumulativeIndex_ - lastCumulativeIndex_};

		// resolution of the possible index bit-field overload by padding with 0 values
		const uint64_t nBFmax = firstMoveDiff / maxIdxBitfield_;
		std::vector<uint32_t> wholeBitField(nBFmax, padding_);
		chunkedAppend(wholeBitField, matrix_);
		firstMoveDiff = firstMoveDiff % maxIdxBitfield_;

		uint32_t firstElement = static_cast<uint32_t>(firstMoveDiff) << valueSize_;
		const uint32_t firstValue{toMerge.matrix_.front() & valueMask_};
		firstElement           |= firstValue;
		toMerge.matrix_.front() = firstElement;

		chunkedAppend(toMerge.matrix_, matrix_);

		lastCumulativeIndex_ = toMerge.lastCumulativeIndex_;
		toMerge.matrix_.clear();
		toMerge.firstCumulativeIndex_ = 0;
		toMerge.lastCumulativeIndex_  = 0;
		return;
	}

	const float chunkSizeFloat = std::sqrt( static_cast<float>( matrix_.size() + toMerge.matrix_.size() ) );
	const auto  chunkSizeDiff  = static_cast<std::vector<uint32_t>::difference_type>(chunkSizeFloat);

	// pack the first similarity value into the first full index
	const uint64_t packedThisFirstFullIdx    = (firstCumulativeIndex_ << valueSize_) | static_cast<uint64_t>(matrix_.front() & valueMask_);
	const uint64_t packedToMergeFirstFullIdx = (toMerge.firstCumulativeIndex_ << valueSize_) | static_cast<uint64_t>(toMerge.matrix_.front() & valueMask_);

	uint64_t cumIdx{0};
	std::vector<uint32_t> mergedMatrix;
	std::vector<uint64_t> fullIndexBuffer;
	uint64_t thisCumIdx{packedThisFirstFullIdx};
	uint64_t toMergeCumIdx{packedToMergeFirstFullIdx};
	uint64_t lastCumIdx{firstCumulativeIndex_};        // for re-creating the diff index
	auto convertToFullIdx = [&cumIdx, &fullIndexBuffer](const uint32_t &currDiff) {
		const auto similarityValue{static_cast<uint64_t>(currDiff & valueMask_)};
		const auto idxDiff{static_cast<uint64_t>(currDiff >> valueSize_)};
		cumIdx = (cumIdx >> valueSize_) + idxDiff;
		cumIdx = (cumIdx << valueSize_) | similarityValue;
		fullIndexBuffer.push_back(cumIdx);
	};
	auto convertToDiffs = [&lastCumIdx, &mergedMatrix](uint64_t &currPackedFullIdx) {
		const auto currValue{static_cast<uint32_t>( currPackedFullIdx & static_cast<uint64_t>(valueMask_) )};
		currPackedFullIdx = currPackedFullIdx >> valueSize_;
		auto currDiff     = static_cast<uint32_t>(currPackedFullIdx - lastCumIdx);
		currDiff          = (currDiff << valueSize_) | currValue;
		lastCumIdx        = currPackedFullIdx;
		mergedMatrix.push_back(currDiff);
	};

	while (std::max( matrix_.size(), toMerge.matrix_.size() ) > 0) {
		auto fiMidPosition{std::distance( fullIndexBuffer.begin(), fullIndexBuffer.end() )};
		const auto thisRemainingSize{std::distance( matrix_.cbegin(), matrix_.cend() )};
		const auto thisChunkEndIt{std::next( matrix_.cbegin(), std::min(chunkSizeDiff, thisRemainingSize) )};
		cumIdx = thisCumIdx;

		std::for_each(matrix_.cbegin(), thisChunkEndIt, convertToFullIdx);
		matrix_.erase(matrix_.begin(), thisChunkEndIt);
		thisCumIdx = cumIdx;
		std::inplace_merge(
			fullIndexBuffer.begin(), 
			std::next(fullIndexBuffer.begin(), fiMidPosition), 
			fullIndexBuffer.end(), 
			[](const uint64_t &firstPacked, const uint64_t &secondPacked) {
				return (firstPacked >> valueSize_) < (secondPacked >> valueSize_);
			}
		);
		fiMidPosition = std::distance( fullIndexBuffer.begin(), fullIndexBuffer.end() );
		cumIdx        = toMergeCumIdx;
		const auto toMergeRemainingSize{std::distance( toMerge.matrix_.cbegin(), toMerge.matrix_.cend() )};
		const auto toMergeChunkEndIt{std::next( toMerge.matrix_.cbegin(), std::min(chunkSizeDiff, toMergeRemainingSize) )};
		std::for_each(toMerge.matrix_.cbegin(), toMergeChunkEndIt, convertToFullIdx);
		toMerge.matrix_.erase(toMerge.matrix_.begin(), toMergeChunkEndIt);
		toMergeCumIdx = cumIdx;

		std::inplace_merge(
			fullIndexBuffer.begin(), 
			std::next(fullIndexBuffer.begin(), fiMidPosition), 
			fullIndexBuffer.end(), 
			[](const uint64_t &firstPacked, const uint64_t &secondPacked) {
				return (firstPacked >> valueSize_) < (secondPacked >> valueSize_);
			}
		);
		auto lastUnqIt = std::unique(
			fullIndexBuffer.begin(),
			fullIndexBuffer.end(),
			[](const uint64_t &firstPacked, const uint64_t &secondPacked) {
				return (firstPacked >> valueSize_) == (secondPacked >> valueSize_);
			}
		);
		fullIndexBuffer.erase( lastUnqIt, fullIndexBuffer.end() );

		// no subsequent matrix blocks can have full indexes smaller than min(thisCumIdx, toMergeCumIdx)
		const uint64_t lowerBoundCumIdx{std::min(thisCumIdx >> valueSize_, toMergeCumIdx >> valueSize_)};
		auto fiMoveIdxIt = std::lower_bound(
			fullIndexBuffer.begin(),
			fullIndexBuffer.end(),
			lowerBoundCumIdx,
			[](const uint64_t &currPckIdx, const uint64_t &lbValue) {
				return (currPckIdx >> valueSize_) <= lbValue;
			}
		);
		std::for_each(
			fullIndexBuffer.begin(),
			fiMoveIdxIt,
			convertToDiffs
		);
		fullIndexBuffer.erase(fullIndexBuffer.begin(), fiMoveIdxIt);
	}

	// move over any straggling full index values
	std::for_each(fullIndexBuffer.begin(),
		fullIndexBuffer.end(),
		convertToDiffs
	);
	fullIndexBuffer.clear();
	std::swap(mergedMatrix, matrix_);
	lastCumulativeIndex_ = std::max(toMerge.lastCumulativeIndex_, lastCumulativeIndex_);

	// complete the move by resetting all members of toMerge
	toMerge.firstCumulativeIndex_ = 0;
	toMerge.lastCumulativeIndex_  = 0;
}

void SimilarityMatrix::save(const std::string &outFileName, const size_t &nThreads, const std::string &locusNameFile) const {
	if ( matrix_.empty() ) {
		return;
	}
	// determine how many line strings can fit in RAM
	const RowColIdx lastIdxPair{recoverRCindexes(lastCumulativeIndex_)};
	const std::string lastEntry{std::to_string(lastIdxPair.iRow + 1) + "\t" + std::to_string(lastIdxPair.jCol + 1) + "\t" + stringLookUp_[0] + "\n"};
	const size_t maxInRAM = getAvailableRAM() / ( static_cast<size_t>(2) * lastEntry.size() );      // use half to leave resources for other operations
	const size_t nChunks  = std::max(matrix_.size() / maxInRAM, 1UL);
	std::vector<size_t> chunkSizes{makeChunkSizes(matrix_.size(), nChunks)};

	std::fstream outStream;
	outStream.open(outFileName, std::ios::out | std::ios::binary | std::ios::app);
	uint64_t runningFullIdx{firstCumulativeIndex_};
	std::vector<uint32_t>::difference_type cumChunkSize{0};
	for (const auto &eachChunkSize : chunkSizes) {
		const size_t actualNthreads = std::min(eachChunkSize, nThreads);
		std::vector<size_t> threadChunkSizes{makeChunkSizes(eachChunkSize, actualNthreads)};
		std::vector< std::pair<std::vector<uint32_t>::const_iterator, std::vector<uint32_t>::const_iterator> > threadPairs;

		// pre-calculate cumulative index ranges
		std::vector<uint64_t> cumIndexes;
		std::vector<uint32_t>::difference_type cumThreadChunkSize{cumChunkSize};
		std::for_each(
			threadChunkSizes.cbegin(),
			threadChunkSizes.cend(),
			[&threadPairs, &cumThreadChunkSize, &runningFullIdx, &cumIndexes, this](const size_t &chunkSize) {
				std::pair<std::vector<uint32_t>::const_iterator, std::vector<uint32_t>::const_iterator> tmpPair{
					matrix_.cbegin() + cumThreadChunkSize,
					matrix_.cbegin() + cumThreadChunkSize + static_cast<std::vector<uint32_t>::difference_type>(chunkSize)
				};
				cumIndexes.push_back(runningFullIdx);
				for (auto matIt = tmpPair.first; matIt != tmpPair.second; ++matIt) {
					runningFullIdx = recoverFullVIdx(matIt, runningFullIdx);
				}
				threadPairs.emplace_back(tmpPair);
				cumThreadChunkSize += static_cast<std::vector<uint32_t>::difference_type>(chunkSize);
			}
		);

		std::vector<std::string> outStrings(actualNthreads);
		std::vector<std::string> locusNames;
		if ( !locusNameFile.empty() ) {
			locusNames = getLocusNames(locusNameFile);
			const RowColIdx lastRC{recoverRCindexes(lastCumulativeIndex_)};
			// testing only the row index, since the column index must be smaller
			if ( lastRC.iRow >= locusNames.size() ) {
				throw std::string("ERROR: number of rows exceeds locus name count in ") + std::string( static_cast<const char*>(__PRETTY_FUNCTION__) );
			}
		}
		std::vector< std::future<void> > tasks;
		tasks.reserve(actualNthreads);
		size_t iThread{0};
		std::for_each(
			threadPairs.cbegin(),
			threadPairs.cend(),
			[&iThread, &tasks, &outStrings, &cumIndexes, &locusNames, this](auto pairIt) {
				tasks.emplace_back(
					std::async(
						[iThread, &cumIndexes, pairIt, &outStrings, &locusNames, this]{
							outStrings.at(iThread) = stringify_(pairIt.first, pairIt.second, cumIndexes.at(iThread), locusNames);
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

std::string SimilarityMatrix::stringify_(std::vector<uint32_t>::const_iterator start, std::vector<uint32_t>::const_iterator end,
									const uint64_t &startCumulativeIndex, const std::vector<std::string> &locusNames) {
	std::string outString;
	uint64_t runningFullIdx{startCumulativeIndex};
	if ( locusNames.empty() ) {
		for (auto matIt = start; matIt != end; ++matIt) {
			const RowColIdx currentPair{recoverRCindexes(matIt, runningFullIdx)};
			if ( ( (*matIt) & valueMask_) != valueMask_ ) { // skip the padding value
				const std::string similarityValue = stringLookUp_.at( static_cast<size_t>( (*matIt) & valueMask_ ) );
				outString += std::to_string(currentPair.iRow + 1) + "\t" + std::to_string(currentPair.jCol + 1) + "\t" + similarityValue + "\n";
			}
		}
		return outString;
	}
	for (auto matIt = start; matIt != end; ++matIt) {
		const RowColIdx currentPair{recoverRCindexes(matIt, runningFullIdx)};
		if ( ( (*matIt) & valueMask_) != valueMask_ ) { // skip the padding value
			const std::string similarityValue = stringLookUp_.at( static_cast<size_t>( (*matIt) & valueMask_ ) );
			outString += locusNames[currentPair.iRow] + "\t" + locusNames[currentPair.jCol] + "\t" + similarityValue + "\n";
		}
	}
	return outString;
}

void SimilarityMatrix::insert_(const FullIdxValue &fullIndexWithSimilarity) {
	if ( matrix_.empty() ) {
		firstCumulativeIndex_ = fullIndexWithSimilarity.fullIdx;
		lastCumulativeIndex_  = fullIndexWithSimilarity.fullIdx;
		matrix_.push_back( static_cast<uint32_t>(fullIndexWithSimilarity.quantSimilarity) );
		return;
	}
	// redundant elements at either end
	if ( (fullIndexWithSimilarity.fullIdx == firstCumulativeIndex_) ||
			(fullIndexWithSimilarity.fullIdx == lastCumulativeIndex_) ) {
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
