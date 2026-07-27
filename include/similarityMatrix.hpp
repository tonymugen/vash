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
 * \version 0.3
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
#include <functional>

namespace BayesicSpace {
	struct RowColIdx;
	struct FullIdxValue;
	struct JaccardPair;
	struct DiffElementPair;
	struct FullIdxTrio;
	struct InOutFileNames;
	class RowColCursor;
	class SimilarityMatrix;
	class SimilarityMatrixSink;

	/** \brief Row and column index pair */
	struct RowColIdx {
		/** \brief Row index */
		uint32_t iRow;
		/** \brief Column index */
		uint32_t jCol;
	};
	/** \brief Full vectorized index and similarity value */
	struct FullIdxValue {
		/** \brief Full index of a vectorized triangular matrix */
		uint64_t fullIdx;
		/** \brief Quantized similarity value */
		uint8_t quantSimilarity;
	};
	/** \brief Pair of integers to calculate Jaccard similarity */
	struct JaccardPair {
		/** \brief Intersection size */
		uint64_t nIntersect;
		/** \brief Union size */
		uint64_t nUnion;
	};

	/** \brief Counts parametrizing a parallel or memory-bounded operation
	 *
	 * Groups a driving count with an associated upper bound, so the two are passed together and
	 * cannot be transposed at a call site. The meaning of each field is operation-specific: e.g.
	 * `parallelBuild` reads `count` as the number of blocks and `ceiling` as the thread cap, while
	 * `SimilarityMatrixSink` reads `count` as the number of save threads and `ceiling` as the
	 * element budget.
	 */
	struct WorkloadLimits {
		/** \brief Driving count (e.g. number of blocks or threads) */
		size_t count;
		/** \brief Associated upper bound (e.g. thread ceiling or element budget) */
		size_t ceiling;
	};

	/** \brief Input and output file names
	 *
	 * Groups input and output file names.
	 */
	struct InOutFileNames {
		/** \brief Input file name */
		std::string inputFileName;
		/** \brief Output file name */
		std::string outputFileName;
	};

	/** \brief Recover row and column indexes
	 *
	 * Recovers the row and column indexes from the matrix element.
	 *
	 * \param[in] vecIdx index into the vectorized matrix
	 * \return row and column index pair
	 */
	[[nodiscard]] RowColIdx recoverRCindexes(const uint64_t &vecIdx) noexcept;

	/** \brief Ascending-order row and column cursor
	 *
	 * Recovers row and column indexes from vectorized matrix indexes, like `recoverRCindexes()`,
	 * but for a range of indexes visited in ascending order. Indexes must be presented to
	 * `advanceTo()` in non-decreasing order, none of them smaller than the index the cursor was
	 * constructed with; the results are otherwise undefined. Amortized over such a range this is
	 * cheaper per element than `recoverRCindexes()`, which remains the way to resolve an isolated
	 * index. A cursor holds the position of one traversal, so concurrently visited ranges (for
	 * example the slices of a threaded save) each need their own.
	 */
	class RowColCursor {
	public:
		/** \brief Default constructor (deleted) */
		RowColCursor() = delete;
		/** \brief Constructor with the first index of a range
		 *
		 * \param[in] firstVecIdx smallest vectorized index the cursor will be presented with
		 */
		explicit RowColCursor(const uint64_t &firstVecIdx) noexcept;
		/** \brief Copy constructor
		 *
		 * \param[in] toCopy object to copy
		 */
		RowColCursor(const RowColCursor &toCopy) noexcept = default;
		/** \brief Copy assignment operator
		 *
		 * \param[in] toCopy object to copy
		 * \return `RowColCursor` object
		 */
		RowColCursor& operator=(const RowColCursor &toCopy) noexcept = default;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		RowColCursor(RowColCursor &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `RowColCursor` object
		 */
		RowColCursor& operator=(RowColCursor &&toMove) noexcept = default;
		/** \brief Destructor */
		~RowColCursor() = default;

		/** \brief Row and column indexes of the next element
		 *
		 * Advances the cursor to `vecIdx` and returns the corresponding index pair.
		 * The index must be no smaller than the one passed to the previous call.
		 *
		 * \param[in] vecIdx index into the vectorized matrix
		 * \return row and column index pair
		 */
		[[nodiscard]] RowColIdx advanceTo(const uint64_t &vecIdx) noexcept;
	private:
		/** \brief Current row index */
		uint64_t row_;
		/** \brief Vectorized index of the first element of the current row */
		uint64_t rowStart_;
		/** \brief Vectorized index of the first element of the next row */
		uint64_t nextRowStart_;
	};

	/** \brief Build a similarity matrix from independent blocks in parallel
	 *
	 * Computes `blockCounts.count` matrix blocks concurrently and consolidates them into a single object.
	 *
	 * \param[in] blockCounts block count (`count`) and thread ceiling (`ceiling`)
	 * \param[in] blockToMatrix callable mapping a block index to its `SimilarityMatrix`
	 * \return consolidated `SimilarityMatrix`
	 */
	[[nodiscard]] SimilarityMatrix parallelBuild(const WorkloadLimits &blockCounts, const std::function<SimilarityMatrix(size_t)> &blockToMatrix);

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
		SimilarityMatrix() noexcept = default;
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
		[[nodiscard]] static size_t elementSize() noexcept { return sizeof(uint64_t); };
		/** \brief Object size in bytes 
		 *
		 * \return object size in bytes
		 */
		[[nodiscard]] size_t objectSize() const noexcept {
			size_t stringBytes{0};
			for (const auto &eachBuffer : saveBuffers_) {
				stringBytes += eachBuffer.capacity();
			}
			return	( elementSize() * matrix_.size() ) + stringBytes;
		};
		/** \brief Number of elements in the matrix
		 *
		 * \return number of elements
		 */
		[[nodiscard]] size_t nElements() const noexcept { return matrix_.size(); };
		/** \brief Reserve element capacity
		 *
		 * Pre-allocates storage for at least `nElements` matrix elements so that subsequent
		 * insertions up to that count do not reallocate. Used to claim a memory budget up front,
		 * before other allocations reduce the available contiguous space. `nElements` is the 3/4
		 * (matrix) share of the memory split; the derived string-scratch budget for `save()` (the
		 * remaining 1/4) is set here as well.
		 *
		 * \param[in] nElements number of matrix elements to reserve capacity for
		 */
		void reserve(const size_t &nElements);
		/** \brief Remove all elements
		 *
		 * Drops every stored element but retains the allocated capacity, so a reserved buffer
		 * can be refilled without reallocating.
		 */
		void clear() noexcept { matrix_.clear(); };
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
		/** \brief Append another matrix without ordering
		 *
		 * Moves the packed elements of `toAppend` onto the end of this matrix and clears `toAppend`,
		 * without restoring the sorted, de-duplicated invariant. Use to accumulate many pre-sorted
		 * blocks cheaply, avoiding the repeated full-vector rebuilds of `merge()`; the invariant must
		 * then be restored with a single `sortAndDeduplicate()` call before any `merge()` or `save()`.
		 *
		 * \param[in,out] toAppend matrix whose elements are moved in and then cleared
		 */
		void append(SimilarityMatrix &toAppend);
		/** \brief Restore the sorted, de-duplicated invariant
		 *
		 * Sorts the packed elements by vectorized index and drops duplicate indexes.
		 * Needed only to finalize a sequence of `append()` calls; all other mutators
		 * keep the invariant.
		 */
		void sortAndDeduplicate();
		/** \brief Merge pre-sorted, de-duplicated runs into one matrix
		 *
		 * Merges `runs`, each already sorted and de-duplicated (every block built through `insert()`
		 * is), into a single sorted, de-duplicated matrix; an index shared across runs is collapsed to
		 * one element. There is no full re-sort: the merge cost is \f$ O(N \log k) \f$ in the total
		 * element count \f$ N \f$ and run count \f$ k \f$, versus \f$ O(N \log N) \f$ for a re-sort.
		 *
		 * \param[in,out] runs individually sorted, de-duplicated matrices; all are cleared and freed
		 * \param[in] maxThreads worker-thread cap used to size the partitioning; 0 uses the hardware concurrency
		 * \return merged sorted, de-duplicated matrix
		 */
		[[nodiscard]] static SimilarityMatrix mergeSortedRuns(std::vector<SimilarityMatrix> &runs, size_t maxThreads = 0);
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
		 * Sizes its string scratch from currently-available RAM; prefer the reusable-buffer overload
		 * when the memory budget is managed externally.
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
		/** \brief Per-thread string scratch buffers for saving to file (reused across chunks)
		 *
		 * Scratch state for the logically-const `save()`, hence `mutable`.
		 */
		mutable std::vector<std::string> saveBuffers_;
		/** \brief Combined byte budget for the save string buffers (0 until `reserve()` is called) */
		size_t saveBufferBudget_{0};

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
		/** \brief Length of every `stringLookUp_` entry
		 *
		 * All quantized values render to the same fixed width, so the entries can be appended
		 * without measuring them and output lines have a constant-width value field.
		 */
		static const size_t valueStringLength_;
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
		 * Enables multi-threaded saving to file.
		 * The `target` string is cleared first (its capacity is retained, enabling buffer reuse).
		 * The elements in `[start, end)` must be in ascending index order.
		 *
		 * \param[in] start start iterator for the matrix
		 * \param[in] end end iterator for the matrix
		 * \param[in] locusNames locus name vector
		 * \param[out] target string the output is written into
		 */
		static void stringify_(std::vector<uint64_t>::const_iterator start, std::vector<uint64_t>::const_iterator end,
								const std::vector<std::string> &locusNames, std::string &target);
		/** \brief Insert a value (updating the index) 
		 *
		 * Inserts a new value into the matrix according to the full vectorized matrix index.
		 *
		 * \param[in] indexWithSimilarity full index and the corresponding quantized similarity
		 */
		void insert_(const FullIdxValue &indexWithSimilarity);
	};

	/** \brief Memory-bounded sink for similarity matrices
	 *
	 * Accumulates `SimilarityMatrix` blocks into an internal buffer and streams them to a file,
	 * keeping resident memory below a fixed element budget set at construction. The buffer's
	 * capacity is reserved up front, so the budget is claimed while memory is still available
	 * and the buffer never reallocates while filling. When adding a block would exceed the
	 * budget, the buffer is de-duplicated and, if still over budget, saved to file and cleared.
	 *
	 * Because saved elements leave memory, duplicate index pairs that recur in blocks added
	 * after a flush cannot be de-duplicated against the already-saved contents. Blocks are
	 * therefore flushed as late as possible (only when genuinely over budget) to minimize such
	 * cross-flush duplication in the output file.
	 */
	class SimilarityMatrixSink {
	public:
		/** \brief Default constructor (deleted) */
		SimilarityMatrixSink() = delete;
		/** \brief Constructor
		 *
		 * Reserves capacity for `saveLimits.ceiling` matrix elements immediately. Each block passed to
		 * `add()` must contain no more than that many elements so the reserved buffer is never exceeded.
		 *
		 * \param[in] fileNames output file name and (optional) locus-name input file name
		 * \param[in] saveLimits number of save threads (`count`) and element budget (`ceiling`); the
		 *            budget is also the reserved buffer capacity
		 */
		SimilarityMatrixSink(const InOutFileNames &fileNames, const WorkloadLimits &saveLimits);
		/** \brief Copy constructor (deleted) */
		SimilarityMatrixSink(const SimilarityMatrixSink &toCopy) = delete;
		/** \brief Copy assignment operator (deleted) */
		SimilarityMatrixSink& operator=(const SimilarityMatrixSink &toCopy) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		SimilarityMatrixSink(SimilarityMatrixSink &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to move
		 * \return `SimilarityMatrixSink` object
		 */
		SimilarityMatrixSink& operator=(SimilarityMatrixSink &&toMove) noexcept = default;
		/** \brief Destructor */
		~SimilarityMatrixSink() = default;

		/** \brief Add a block
		 *
		 * Appends the elements of `block` to the buffer, first flushing to file if the combined
		 * size would exceed the budget. Clears `block`.
		 *
		 * \param[in,out] block matrix whose elements are moved in and then cleared
		 */
		void add(SimilarityMatrix &block);
		/** \brief Flush remaining buffered elements
		 *
		 * De-duplicates and saves whatever remains in the buffer, then clears it. Call once after
		 * the last `add()`.
		 */
		void finalize();
		/** \brief Number of buffered elements
		 *
		 * \return element count currently held in the buffer (may include not-yet-de-duplicated pairs)
		 */
		[[nodiscard]] size_t bufferedElements() const noexcept { return buffer_.nElements(); };
	private:
		/** \brief Accumulation buffer
		 *
		 * Holds the matrix data (reserved to the 3/4 matrix share of the memory budget) and owns the
		 * per-thread save string scratch (the remaining 1/4, sized inside `SimilarityMatrix::reserve`).
		 */
		SimilarityMatrix buffer_;
		/** \brief Output file name */
		std::string outFileName_;
		/** \brief Locus-name input file name (empty if unused) */
		std::string locusNameFile_;
		/** \brief Number of threads passed to the buffer's `save()` */
		size_t nThreads_;
		/** \brief Element budget and reserved matrix-buffer capacity (3/4 of the memory budget) */
		size_t maxElements_;
		/** \brief De-duplicate and save the buffer, then clear it (capacity retained) */
		void flush_();
	};
}
