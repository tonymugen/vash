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
 * Definitions and interface documentation for a compact representation of a (possibly sparse) similarity matrix.
 *
 */

#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <cstddef>

namespace BayesicSpace {
	struct IndexedSimilarity;
	struct MatrixInitializer;
	struct RowColIdx;
	class SimilarityMatrix;

	/** \brief Index-similarity pair 
	 *
	 * The index is a distance from the preceding value in a vectorized similarity matrix.
	 * Similarity values are quantized to fit in a single byte.
	 */
	struct IndexedSimilarity {
		uint32_t diffIndex;
		uint8_t  similarity;
	};

	/** \brief Object with constants to initialize a `SimilarityMatrix` object
	 *
	 * `previousCumulativeIndex` is the overall index of the last element of the previous portion
	 * of the overall matrix.
	 * `nRows` is the number of rows (or columns, since `SimilarityMatrix` is square) of the overall matrix.
	 */
	struct MatrixInitializer {
		uint64_t previousCumulativeIndex;
		uint32_t nRows;
	};

	/** \brief Row and column index pair */
	struct RowColIdx {
		uint32_t iRow;
		uint32_t jCol;
	};

	/** \brief Similarity matrix
	 *
	 * A representation of a symmetric similarity matrix, excluding the diagonal.
	 * Stores the lower triangle by row.
	 * The representation is efficient if the matrix is sparse.
	 * The object can be a section of a larger matrix. In this case,
	 * a start index (that addresses the end of the vectorized preceding chunk of the larger matrix)
	 * can also be provided.
	 */
	class SimilarityMatrix {
	public:
		/** \brief Default constructor */
		SimilarityMatrix() noexcept : previousCumulativeIndex_{0}, lastCumulativeIndex_{0}, nRows_{0} {};
		/** \brief Constructor
		 *
		 * \param[in] initializer initializer object
		 */
		SimilarityMatrix(const MatrixInitializer &initializer) noexcept :
			previousCumulativeIndex_{initializer.previousCumulativeIndex}, lastCumulativeIndex_{initializer.previousCumulativeIndex}, nRows_{initializer.nRows} {};
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		SimilarityMatrix(const SimilarityMatrix &toCopy) = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `SimilarityMatrix` object
		 */
		SimilarityMatrix& operator=(const SimilarityMatrix &toCopy) = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		SimilarityMatrix( SimilarityMatrix &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `SimilarityMatrix` object
		 */
		SimilarityMatrix& operator=( SimilarityMatrix &&toMove) noexcept = default;
		/** \brief Destructor */
		~SimilarityMatrix() = default;

		/** \brief Object size in bytes */
		[[nodiscard]] size_t size() const noexcept { 
			return	sizeof(IndexedSimilarity) * matrix_.size() + // matrix size
					sizeof(uint32_t) +                           // nRows_ size
					2 * sizeof(uint64_t);                        // index sizes
		};
		/** \brief Append a value (updating the index) 
		 *
		 * Appending the data according to the provided row and column indexes.
		 * The indexes should belong the to overall matrix.
		 * Row index must be larger than the column index. If not, the values are swapped.
		 * If the indexes are equal or fall in front of the last included value, throws an exception.
		 *
		 * \param[in] rowColPair row and index pair
		 * \param[in] quantSimilarity quantized similarity value
		 */
		void append(const RowColIdx &rowColPair, uint8_t quantSimilarity);
	private:
		/** \brief Vectorized data representation */
		std::vector<IndexedSimilarity> matrix_;
		/** \brief Last cumulative index of the previous sub-matrix */
		uint64_t previousCumulativeIndex_;
		/** \brief Cumulative (for the whole matrix) value of the last index */
		uint64_t lastCumulativeIndex_;
		/** \brief Number of rows (or columns) in the full matrix */
		uint32_t nRows_;
		/** \brief Floating point look-up table
		 *
		 * Used to substitute floating-point values that correspond to the
		 * quantized eight-bit representation in `IndexedSimilarity`.
		 */
		static const std::array<float, 256> floatLookUp_;
	};
}
