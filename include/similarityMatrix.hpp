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
#include <string>
#include <cstdint>
#include <cstddef>

namespace BayesicSpace {
	struct RowColIdx;
	struct FullIdxValue;
	struct JaccardPair;
	struct DiffElementPair;
	struct FullIdxTrio;
	class SimilarityMatrix;

	/** \brief Row and column index pair */
	struct RowColIdx {
		uint32_t iRow;
		uint32_t jCol;
	};
	/** \brief Full vectorized index and similarity value */
	struct FullIdxValue {
		uint64_t fullIdx;
		uint8_t quantSimilarity;
	};
	/** \brief Pair of integers to calculate Jaccard similarity */
	struct JaccardPair {
		uint64_t nIntersect;
		uint64_t nUnion;
	};

	/** \brief Append one vector to another by chunks
	 * 
	 * Moves the contents of the source vector to the end of the target vector.
	 * Uses \f$ \sqrt{N} \f$, where \f$ N \f$ is the size of the source vector,
	 * extra memory. The source vector is cleared.
	 *
	 * \param[in] source the source vector, is cleared as a result
	 * \param[in,out] target the vector accepting the data from `source`
	 */
	void chunkedAppend(std::vector<uint64_t> &source, std::vector<uint64_t> &target);
	/** \brief Recover row and column indexes 
	 *
	 * Recovers the row and column indexes from the matrix element.
	 *
	 * \param[in] vecIdx index into the vectorized matrix
	 * \return row and column index pair
	 */
	[[gnu::warn_unused_result]] RowColIdx recoverRCindexes(const uint64_t &vecIdx) noexcept;

	/** \brief Similarity matrix
	 *
	 * A representation of a square symmetric similarity matrix, excluding the diagonal.
	 * Stores only the values present in the lower triangle by row.
	 * The representation is memory efficient if the matrix is sparse and attempts
	 * to strike a compromise between memory use and matrix manipulation speed.
	 */
	class SimilarityMatrix {
	public:
		/** \brief Default constructor */
		SimilarityMatrix() noexcept  = default;
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

		/** \brief Matrix element size
		 *
		 * \return matrix element size in bytes
		 */
		[[gnu::warn_unused_result]] static size_t elementSize() noexcept {return sizeof(uint64_t);};
		/** \brief Object size in bytes 
		 *
		 * \return object size in bytes
		 */
		[[gnu::warn_unused_result]] size_t objectSize() const noexcept { 
			return	elementSize() * matrix_.size();
		};
		/** \brief Number of elements in the matrix
		 *
		 * \return number of elements
		 */
		[[gnu::warn_unused_result]] size_t nElements() const noexcept {return matrix_.size();};
		/** \brief Insert a value (updating the index) 
		 *
		 * Inserts a new value into the matrix. Addresses the lower triangle of the similarity matrix,
		 * therefore the row index must be larger than the column index. If not, the values are swapped.
		 * If the indexes are equal or the row index is 0, throws an exception.
		 * Inserts a quantized value of the Jaccard similarity calculated from the intersection and union counts provided.
		 *
		 * \param[in] rowColPair row and index pair
		 * \param[in] jaccardCounts intersection and union counts for Jaccard similarity
		 */
		void insert(const RowColIdx &rowColPair, const JaccardPair &jaccardCounts);
		/** \brief Merge two matrices 
		 *
		 * Merge a matrix with the current object, destroying the donor object.
		 * Duplicated indexes are discarded even if they differ in similarity values.
		 *
		 * \param[in] toMerge object to merge
		 */
		void merge(SimilarityMatrix &toMerge);
		/** \brief Save to file 
		 *
		 * Uses multi-threaded data prep to speed up saving.
		 * If the output file already exists, appends to it.
		 *
		 * \param[in] outFileName output file name
		 * \param[in] nThreads number of threads
		 * \param[in] locusNameFile name of the file with locus names (empty by default)
		 */
		void save(const std::string &outFileName, const size_t &nThreads, const std::string &locusNameFile = "") const;
	private:
		/** \brief Vectorized data representation 
		 *
		 * The index of the vectorized (by row) lower triangle and value are packed into 64-bit integers.
		 * The first byte is the quantized similarity value (indexing the look-up table).
		 * The rest encode the vectorized index of the element.
		 */
		std::vector<uint64_t> matrix_;

		// static members
		/** \brief Floating point look-up table
		 *
		 * Used to substitute floating-point values that correspond to the
		 * quantized representation in the `matrix_`.
		 */
		static const std::array<float, 256> floatLookUp_;
		/** \brief String look-up table
		 *
		 * Used to substitute string representations (for display) of the floating-point values
		 * that correspond to the quantized representation in the `matrix_`.
		 */
		static const std::array<const char*, 256> stringLookUp_;
		/** \brief Maximal index bit-field value */
		static const uint64_t maxIdxBitfield_;
		/** \brief Maximal row and column value 
		 *
		 * Depends on `maxIdxBitfield_`
		 */
		static const uint32_t maxRowColValue_;
		/** \brief Mask that isolates the value bit-field */
		static const uint64_t valueMask_;
		/** \brief Maximal 8-bit index into the `float` value table */
		static const uint64_t maxValueIdx_;
		/** \brief Size of the value bit-field in bits */
		static const uint64_t valueSize_;
		/** \brief Convert matrix data to string with locus names
		 *
		 * Construct a string from a portion of the matrix for saving.
		 * Add locus names if the `locusNames` vector is not empty.
		 * Enables multi-threaded saving to file, since conversion to string is the bottleneck for `fstream`.
		 *
		 * \param[in] start start iterator for the matrix
		 * \param[in] end end iterator for the matrix
		 * \param[in] locusNames locus name vector
		 *  \return output string
		 */
		[[gnu::warn_unused_result]] static std::string stringify_(std::vector<uint64_t>::const_iterator start, std::vector<uint64_t>::const_iterator end,
								const std::vector<std::string> &locusNames);
		/** \brief Insert a value (updating the index) 
		 *
		 * Inserts a new value into the matrix according to the full vectorized matrix index.
		 *
		 * \param[in] indexWithSimilarity full index and the corresponding quantized similarity
		 */
		void insert_(const FullIdxValue &indexWithSimilarity);
	};
}
