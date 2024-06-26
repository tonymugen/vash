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
 * \version 0.5
 *
 * Definitions and interface documentation for classes that take binary variant files and generate lossy summaries with hashing.
 *
 */

#pragma once

#include <cstddef>
#include <fstream>
#include <vector>
#include <utility>        // for std::pair
#include <string>
#include <thread>

#include "similarityMatrix.hpp"

namespace BayesicSpace {
	struct LocationWithLength;
	struct CountAndSize;
	struct IndividualAndSketchCounts;
	struct BedDataStats;
	struct InOutFileNames;
	struct SparsityParameters;
	struct HashGroup;
	struct HashGroupItPairCount;
	class GenoTableBin;
	class GenoTableHash;

	/** \brief Window location and extent
	 *
	 * Groups the start index and length (in number of elements) of a window spanning container elements.
	 */
	struct LocationWithLength {
		size_t start;
		size_t length;
	};

	/** \brief Number of items and their size
	 *
	 * Groups the number of items with size of each.
	 */
	struct CountAndSize {
		size_t count;
		size_t size;
	};

	/** \brief Number of individuals and sketches
	 *
	 * Groups the number of individuals with sketch numbers for hashing.
	 */
	struct IndividualAndSketchCounts {
		uint32_t nIndividuals;
		uint16_t kSketches;
	};

	/** \brief Attributes of _.bed_ format loci
	 *
	 * Data attributes of locus groups for reading _.bed_ files.
	 */
	struct BedDataStats {
		size_t firstLocusIdx;
		size_t nLociPerThread;
		size_t nBytesPerLocus;
		size_t nBytesToRead;
		size_t nLociToRead;
		size_t nMemChunks;     // number of chunks to read into memory
	};
	/** \brief Input and output file names
	 *
	 * Groups input and output file names.
	 */
	struct InOutFileNames {
		std::string inputFileName;
		std::string outputFileName;
	};
	/** \brief LD matrix sparsity parameters
	 *
	 * The number of rows per band and similarity value cut-off control
	 * sparsity of the similarity matrix. Greater values for either result
	 * in fewer elements included in the matrix.
	 */
	struct SparsityParameters{
		size_t nRowsPerBand;
		float similarityCutOff;
	};
	/** \brief Hash-derived group
	 *
	 * Includes the cumulative number of pairs from all previous
	 * (including the current) groups.
	 */
	struct HashGroup {
		uint64_t cumulativeNpairs;
		std::vector<uint32_t> locusIndexes;
	};

	/** \brief Hash group vector iterator and element number
	 *
	 * Used to delimit ranges of locus pairs for chunked processing.
	 */
	struct HashGroupItPairCount {
		// how many pairs already processed
		size_t pairCount{0};
		std::vector<HashGroup>::const_iterator hgIterator;
	};

	/** \brief Class to store binary compressed genotype tables
	 *
	 * Converts genotype data to a lossy compressed binary code.
	 * Genotypes are stored in memory in a one-bit format: bit set for the minor allele, unset for the major.
	 * Bits corresponding to missing data are unset (this is the same as mean imputation), heterozygotes are set with a 50% probability.
	 */
	class GenoTableBin {
	public:
		/** \brief Default constructor */
		GenoTableBin() : nIndividuals_{0}, nLoci_{0}, binLocusSize_{0}, nThreads_{1} {};
		/** \brief Constructor with input file name
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed).
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * If necessary, alleles are re-coded so that the set bit is always the minor allele.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] logFileName name of the log file
		 */
		GenoTableBin(const std::string &inputFileName, const uint32_t &nIndividuals, std::string logFileName) : 
						GenoTableBin( inputFileName, nIndividuals, std::move(logFileName), std::thread::hardware_concurrency() ) {};
		/** \brief Constructor with input file name and thread count
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed).
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * If necessary, alleles are re-coded so that the set bit is always the minor allele.
		 * The number of threads requested is maximum to be used, depending on available system resources.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] logFileName name of the log file
		 * \param[in] nThreads maximal number of threads to use
		 */
		GenoTableBin(const std::string &inputFileName, const uint32_t &nIndividuals, std::string logFileName, const size_t &nThreads);
		/** \brief Constructor with count vector
		 *
		 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * The counts are checked and re-coded if necessary so that set bits represent the minor allele. This function should run faster if the 0 is the major allele homozygote.
		 * While the above values are the norm, any negative number will be interpreted as missing, any odd number as 1, and any (non-0) even number as 2.
		 * The input is a vectorized matrix of genotypes. The original matrix has individuals on rows, and is vectorized by row.
		 *
		 * \param[in] maCounts vector of minor allele numbers
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] logFileName name of the log file
		 */
		GenoTableBin(const std::vector<int> &maCounts, const uint32_t &nIndividuals, std::string logFileName) :
						GenoTableBin( maCounts, nIndividuals, std::move(logFileName), std::thread::hardware_concurrency() ) {};
		/** \brief Constructor with count vector and thread count
		 *
		 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * The counts are checked and re-coded if necessary so that set bits represent the minor allele. This function should run faster if the 0 is the major allele homozygote.
		 * While the above values are the norm, any negative number will be interpreted as missing, any odd number as 1, and any (non-0) even number as 2.
		 * The input is a vectorized matrix of genotypes. The original matrix has individuals on rows, and is vectorized by row.
		 * The number of threads requested is maximum to be used, depending on available system resources.
		 *
		 * \param[in] maCounts vector of minor allele numbers
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] logFileName name of the log file
		 * \param[in] nThreads maximal number of threads to use
		 */
		GenoTableBin(const std::vector<int> &maCounts, const uint32_t &nIndividuals, std::string logFileName, const size_t &nThreads);

		/** \brief Copy constructor (deleted) */
		GenoTableBin(const GenoTableBin &toCopy) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTableBin& operator=(const GenoTableBin &toCopy) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		GenoTableBin(GenoTableBin &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to be moved
		 * \return `GenoTableBin` object
		 */
		GenoTableBin& operator=(GenoTableBin &&toMove) noexcept = default;
		/** \brief Destructor */
		~GenoTableBin() = default;

		/** \brief Save the binary genotype file
		 *
		 * Saves the binary approximate genotype data to a binary file.
		 *
		 * \param[in] outFileName output file name
		 */
		void saveGenoBinary(const std::string &outFileName) const;
		/** \brief All by all Jaccard similarity LD with locus names
		 *
		 * Calculates linkage disequilibrium among all loci using Jaccard similarity and \f$r^2\f$ as the statistics.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * Row and column locus names are also included in the tab-delimited output file.
		 * The lower triangle is vectorized by column (i.e. all correlations of the first locus, then all remaining correlations of the second, etc.).
		 * If the result does not fit in RAM, calculates in blocks and saves to disk periodically.
		 *
		 * \param[in] bimAndLDnames name of the input _.bim_ file that has locus names and the output LD value file name
		 */
		void allJaccardLD( const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks = static_cast<size_t>(1) ) const;
		/** \brief Save the log to a file
		 *
		 * Log file name provided at construction.
		 */
		void saveLogFile() const;
	private:
		/** \brief Binarized genotype table
		 *
		 * Stores one bit per genotype. Heterozygotes are randomly assigned, missing data are assigned 0.
		 */
		std::vector<uint8_t> binGenotypes_;
		/** \brief Number of individuals */
		uint32_t nIndividuals_;
		/** \brief Number of loci */
		uint32_t nLoci_;
		/** \brief Binarized locus size in bytes */
		size_t binLocusSize_;
		/** \brief Maximal number of threads to use */
		size_t nThreads_;
		/** \brief Leading bytes for .bed files */
		static const size_t nMagicBytes_;
		/** \brief One set bit for masking */
		static const uint8_t oneBit_;
		/** \brief Size of one byte in bits */
		static const uint8_t byteSize_;
		/** \brief Number of .bed genotypes per byte */
		static const uint8_t bedGenoPerByte_;
		/** \brief 64 bit word size in bytes */
		static const uint8_t llWordSize_;
		/** \brief Log messages */
		mutable std::string logMessages_;
		/** \brief Log file name */
		std::string logFileName_;
		/** \brief Binarize a range of loci from _.bed_ file input
		 *
		 * Binarizes a range of loci from a vector of input from a _.bed_ file.
		 *
		 * \param[in] bedData vector of _.bed_ file input
		 * \param[in] bedLocusIndRange indexes of the first and last locus in the _.bed_ vector
		 * \param[in] locusSpan window with the overall index of the first locus in the range and locus length
		 */
		void bed2binBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const LocationWithLength &locusSpan);
		/** \brief Multi-threaded binarization of _.bed_ file input 
		 *
		 * Binarizes loci from a _.bed_ format byte-vector using multiple threads.
		 *
		 * \param[in] bedData vector of _.bed_ file input
		 * \param[in] threadRanges vector of locus index ranges, one per thread
		 * \param[in] locusSpan window with the overall index of the first locus in the range and locus length
		 * \return new value of `firstLocusInd`
		 */
		size_t bed2binThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const LocationWithLength &locusSpan);
		/** \brief Wraps _.bed_ file binarization 
		 *
		 * Wraps _.bed_ format lossy conversion to binary (with one bit per locus).
		 *
		 * \param[in] locusGroupStats _.bed_ locus group attributes
		 * \param[in,out] bedStream _.bed_ file to be converted
		 * \return new start individual index
		 */
		size_t bed2bin_(const BedDataStats &locusGroupStats, std::fstream &bedStream);
		/** \brief Binarize minor allele counts in a locus block
		 *
		 * Binarizes a portion of a vector of per-individual minor allele counts (0, 1, or 2; see the count vector constructor documentation for details).
		 *
		 * \param[in] macData vector of minor allele counts
		 * \param[in] locusIndRange locus index range
		 */
		void mac2binBlk_(const std::vector<int> &macData, const std::pair<size_t, size_t> &locusIndRange);
		/** \brief Jaccard similarity in a block of loci
		 *
		 * \param[in] blockRange row/column index pair range
		 * \return `SimilarityMatrix` object with compressed indexed similarity values
		 */
		[[gnu::warn_unused_result]] SimilarityMatrix jaccardBlock_(const std::pair<RowColIdx, RowColIdx> &blockRange) const;
		/** \brief Jaccard similarity between locus pairs using multiple threads
		 *
		 * \param[in] indexPairs vector of row/column index ranges, one per thread
		 * \return `SimilarityMatrix` object with compressed indexed similarity values
		 */
		[[gnu::warn_unused_result]] SimilarityMatrix jaccardThreaded_(const std::vector< std::pair<RowColIdx, RowColIdx> > &indexPairs) const;
		/** \brief Calculate the union and intersection Jaccard similarity pair
		 *
		 * \param[in] rowColumn indexes of the locus pair to compare
		 * \return Intersection and union for Jaccard similarity calculation between two binarized loci
		 */
		[[gnu::warn_unused_result]] JaccardPair makeJaccardPair_(const RowColIdx &rowColumn) const;
	};
	/** \brief Class to store compressed genotype tables
	 *
	 * Provides facilities to store and manipulate compressed genotype tables.
	 * Genotypes are stored in a one-bit format: bit set for the minor allele, unset for the major.
	 * Bits corresponding to missing data are unset (this is the same as mean imputation), heterozygotes are set with a 50% probability.
	 */
	class GenoTableHash {
	public:
		/** \brief Default constructor */
		GenoTableHash() : nIndividuals_{0}, kSketches_{0}, sketchSize_{0}, nLoci_{0}, locusSize_{0}, nFullWordBytes_{0}, nThreads_{1}, emptyBinIdxSeed_{0} {};
		/** \brief Constructor with input file name and thread number
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed) format.
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * If necessary, alleles are re-coded so that the set bit is always the minor allele.
		 * The binary stream is then hashed using a one-permutation hash (OPH; one sketch per locus).
		 * Bits are permuted using the Fisher-Yates-Durstenfeld algorithm.
		 * Filling in empty bins using the Mai _et al._ (2020) algorithm.
		 * The number of threads specified is the maximal that will be used. Actual number depends on system resources.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] indivSketchCounts number of individuals and sketches
		 * \param[in] nThreds maximal number of threads to use
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::string &inputFileName, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, std::string logFileName);
		/** \brief Constructor with input file name
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed) format.
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * If necessary, alleles are re-coded so that the set bit is always the minor allele.
		 * The input is a vectorized matrix of genotypes. The original matrix has individuals on rows, and is vectorized by row.
		 * The binary stream is then hashed using a one-permutation hash (OPH; one sketch per locus).
		 * Bits are permuted using the Fisher-Yates-Durstenfeld algorithm.
		 * Filling in empty bins using the Mai _et al._ (2020) algorithm.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] indivSketchCounts number of individuals and sketches
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::string &inputFileName, const IndividualAndSketchCounts &indivSketchCounts, std::string logFileName) :
							GenoTableHash( inputFileName, indivSketchCounts, std::thread::hardware_concurrency(), std::move(logFileName) ) {};
		/** \brief Constructor with count vector and thread number
		 *
		 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * The counts are checked and re-coded if necessary so that set bits represent the minor allele. This function should run faster if the 0 is the major allele homozygote.
		 * While the above values are the norm, any negative number will be interpreted as missing, any odd number as 1, and any (non-0) even number as 2.
		 * The input is a vectorized matrix of genotypes. The original matrix has individuals on rows, and is vectorized by row.
		 * The binary stream is then hashed using a one-permutation hash (OPH; one sketch per locus).
		 * Bits are permuted using the Fisher-Yates-Durstenfeld algorithm.
		 * Filling in empty bins using the Mai _et al._ (2020) algorithm.
		 * The number of threads specified is the maximal that will be used. Actual number depends on system resources.
		 *
		 * \param[in] maCounts vector of minor allele numbers
		 * \param[in] indivSketchCounts number of individuals and sketches
		 * \param[in] nThreds maximal number of threads to use
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::vector<int> &maCounts, const IndividualAndSketchCounts &indivSketchCounts, const size_t &nThreads, std::string logFileName);
		/** \brief Constructor with count vector
		 *
		 * Input is a vector of minor allele counts (0, 1, or 2) or -9 for missing data.
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * The counts are checked and re-coded if necessary so that set bits represent the minor allele. This function should run faster if the 0 is the major allele homozygote.
		 * While the above values are the norm, any negative number will be interpreted as missing, any odd number as 1, and any (non-0) even number as 2.
		 * The binary stream is then hashed using a one-permutation hash (OPH; one sketch per locus).
		 * Bits are permuted using the Fisher-Yates-Durstenfeld algorithm.
		 * Filling in empty bins using the Mai _et al._ (2020) algorithm.
		 *
		 * \param[in] maCounts vector of minor allele numbers
		 * \param[in] indivSketchCounts number of individuals and sketches
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::vector<int> &maCounts, const IndividualAndSketchCounts &indivSketchCounts, std::string logFileName) :
				GenoTableHash( maCounts, indivSketchCounts, std::thread::hardware_concurrency(), std::move(logFileName) ) {};

		/** \brief Copy constructor (deleted) */
		GenoTableHash(const GenoTableHash &toCopy) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTableHash& operator=(const GenoTableHash &toCopy) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		GenoTableHash(GenoTableHash &&toMove) noexcept = default;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to be moved
		 * \return `GenoTableHash` object
		 */
		GenoTableHash& operator=(GenoTableHash &&toMove) noexcept = default;
		/** \brief Destructor */
		~GenoTableHash() = default;

		/** \brief All by all LD from hashes
		 *
		 * Calculates linkage disequilibrium among all loci using a modified OPH.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * All values belong to the same group. Row and column locus names are also included in the tab-delimited output file.
		 * The lower triangle is vectorized by column (i.e. all correlations of the first locus, then all remaining correlations of the second, etc.).
		 * If `suggestNchunks` is set, processing the data at least in the given number of chunks even if everything fits in RAM.
		 * If the resulting chunks are still too big to fit in RAM, the number is adjusted up.
		 * Otherwise, set the number of chunks automatically.
		 * If the .bim file name is left blank or the file does not exist, base-1 locus indexes are used instead of locus names.
		 *
		 * \param[in] similarityCutOff only save pairs with at least this similarity
		 * \param[in] bimAndLDnames name of the _.bim_ file with locus names and the output LD results file
		 * \param[in] suggestNchunks force processing in chunks
		 */
		void allHashLD( const float &similarityCutOff, const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks = static_cast<size_t>(1) ) const;
		/** \brief Assign groups by linkage disequilibrium (LD)
		 *
		 * The sketch matrix is divided into bands, `nRowsPerBand` rows per band (must be 1 or greater).
		 * Locus pairs are included in the pair hash table if all rows in at least one band match.
		 * The resulting hash table has groups with at least two loci per group (indexed by a hash of the index vector in the group).
		 * Locus indexes are in increasing order within each group. Groups are sorted by first and second locus indexes.
		 * Some locus pairs may end up in more than one group, but no groups are completely identical in locus composition.
		 *
		 * \param[in] nRowsPerBand number of rows per sketch matrix band
		 * \return locus index hash table
		 */
		[[gnu::warn_unused_result]] std::vector<HashGroup> makeLDgroups(const size_t &nRowsPerBand) const;
		/** \brief Assign groups by LD and save to a file with locus names
		 *
		 * Assign groups as above and save locus names with their group IDs to a file.
		 * If the .bim file name is left blank or the file does not exist, base-1 locus indexes are used instead of locus names.
		 *
		 * \param[in] nRowsPerBand number of rows per sketch matrix band
		 * \param[in] bimAndGroupNames _.bim_ and output group file name
		 */
		void makeLDgroups(const size_t &nRowsPerBand, const InOutFileNames &bimAndGroupNames) const;
		/** \brief LD in groups
		 *
		 * Group loci according to LD using the algorithm for `makeLDgroups` and calculate similarity within groups.
		 * Output LD (Jaccard similarity) estimates with group IDs and locus names.
		 * If `suggestNchunks` is set, processing the data at least in the given number of chunks even if everything fits in RAM.
		 * If the resulting chunks are still too big to fit in RAM, the number is adjusted up.
		 * Otherwise, set the number of chunks automatically.
		 * If the .bim file name is left blank or the file does not exist, base-1 locus indexes are used instead of locus names.
		 *
		 * \param[in] sparsityValues `SparsityParameters` object that controls output matrix sparsity
		 * \param[in] similarityCutOff only save pairs with at least this similarity
		 * \param[in] bimAndLDnames _.bim_ and output LD file names
		 * \param[in] suggestNchunks force processing in chunks
		 */
		void ldInGroups(const SparsityParameters &sparsityValues, const InOutFileNames &bimAndLDnames, const size_t &suggestNchunks = static_cast<size_t>(1) ) const;
		/** \brief Save the log to a file
		 *
		 * Log file name provided at construction.
		 */
		void saveLogFile() const;
	private:
		/** \brief Vector of sketches
		 *
		 * A sketch is the position of the first set bit in a bin of permuted bits.
		 */
		std::vector<uint16_t> sketches_;
		/** \brief Number of individuals, possibly padded */
		uint32_t nIndividuals_;
		/** \brief Number of sketches */
		uint16_t kSketches_;
		/** \brief Sketch size */
		uint32_t sketchSize_;
		/** \brief Number of loci */
		uint32_t nLoci_;
		/** \brief Locus size in bytes */
		uint32_t locusSize_;
		/** \brief Number of bytes in divisible by `llWordSize_` */
		size_t nFullWordBytes_;
		/** \brief Maximal number of threads to use */
		size_t nThreads_;
		/** \brief Random index progression seed 
		 *
		 * Seeds the random index progression used to fill empty bins in OPH sketches.
		 * The index set must be the same across loci (although not necessarily the same number is actually used).
		 */
		uint64_t emptyBinIdxSeed_;
		/** \brief Log messages */
		mutable std::string logMessages_;
		/** \brief Log file name */
		std::string logFileName_;
		/** \brief Leading bytes for .bed files */
		static const size_t nMagicBytes_;
		/** \brief One set bit for masking */
		static const uint8_t oneBit_;
		/** \brief Size of one byte in bits */
		static const uint8_t byteSize_;
		/** \brief Number of .bed genotypes per byte */
		static const uint8_t bedGenoPerByte_;
		/** \brief 64 bit word size in bytes */
		static const uint8_t llWordSize_;
		/** \brief Mask to round to the nearest whole-byte count */
		static const uint32_t roundMask_;
		/** 64-bit word with all bits set */
		static const uint64_t allBitsSet_;
		/** \brief 64-bit word size in bits */
		static const size_t wordSizeInBits_;
		/** \brief Value corresponding to an empty token */
		static const uint16_t emptyBinToken_;
		/** \brief Permute bits 
		 *
		 * Permutes individual bits in a vector of bytes according to the provided index vector.
		 * The index vector stores swap addresses of the corresponding positions.
		 *
		 * \param[in] permutationIdx permutation index vector
		 * \param[in,out] binLocus vector of bytes to be permuted
		 *
		 */
		void permuteBits_(const std::vector<size_t> &permutationIdx, std::vector<uint8_t> &binLocus) const;
		/** \brief Single-locus one-permutation hash
		 *
		 * Generates an OPH of a binarized locus. The locus data are modified by the function.
		 *
		 * \param[in] locusInd locus index
		 * \param[in] permutation permutation to be applied to each locus
		 * \param[in,out] binLocus vector of genotypes for a locus
		 */
		void locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint8_t> &binLocus);
		/** \brief OPH from _.bed_ file input
		 *
		 * Hashes a portion of a vector of input from a _.bed_ file that corresponds to a range of loci.
		 *
		 * \param[in] bedData _.bed_ file input
		 * \param[in] bedLocusIndRange range of locus indexes in the block
		 * \param[in] bedLocusSpan position and of the first _.bed_ locus with its size
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in] padIndiv additional individuals, `first` is the placement index, `second` is the index of the individual to add
		 */
		void bed2ophBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const LocationWithLength &bedLocusSpan,
									const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv);
		/** \brief OPH from _.bed_ file input using multiple threads
		 *
		 * Hashes input from a _.bed_ file using multiple threads.
		 *
		 * \param[in] bedData _.bed_ file input
		 * \param[in] threadRanges vector of locus ranges, one per thread
		 * \param[in] bedLocusSpan position of the first _.bed_ locus with its size
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in] padIndiv additional individuals, `first` is the placement index, `second` is the index of the individual to add
		 * \return new value of `firstLocusInd`
		 */
		size_t bed2ophThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const LocationWithLength &bedLocusSpan,
							const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv);
		/** \brief Wraps _.bed_ file to binarization 
		 *
		 * Wraps _.bed_ format hashing.
		 *
		 * \param[in] locusGroupStats _.bed_ locus group attributes
		 * \param[in,out] bedStream _.bed_ file to be converted
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in] padIndiv additional individuals, `first` is the placement index, `second` is the index of the individual to add
		 * \return new start individual index
		 */
		size_t bed2oph_(const BedDataStats &locusGroupStats, std::fstream &bedStream, const std::vector<size_t> &permutation,
							const std::vector< std::pair<size_t, size_t> > &padIndiv);
		/** \brief OPH from minor allele counts
		 *
		 * Hashes a portion of a vector of per-individual minor allele counts (0, 1, or 2; see the count vector constructor documentation for details).
		 * The vector portion corresponds to a block of loci.
		 *
		 * \param[in] macData vector of minor allele counts
		 * \param[in] blockRange range of loci in the block
		 * \param[in] permutation permutation to be applied to each locus 
		 */
		void mac2ophBlk_(const std::vector<int> &macData, const std::pair<size_t, size_t> &blockRange,
				const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv);
		/** \brief OPH from minor allele counts using multiple threads
		 *
		 * Hashes input from a vector of minor allele counts using multiple threads.
		 *
		 * \param[in] macData vector of minor allele counts
		 * \param[in] threadRanges vector of locus ranges, one per thread
		 * \param[in] locusBlock locus block start and size
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in] padIndiv additional individuals, `first` is the placement index, `second` is the index of the individual to add
		 * \return new value of `firstLocusInd`
		 */
		size_t mac2ophThreaded_(const std::vector<int> &macData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const LocationWithLength &locusBlock,
							const std::vector<size_t> &permutation, const std::vector< std::pair<size_t, size_t> > &padIndiv);
		/** \brief Hash-based similarity in a block of loci
		 *
		 * The provided row and column values index the `locusIndexes` vector that translates them to the actual locus indexes.
		 * This is necessary for hash group processing.
		 *
		 * \param[in] blockRange row/column index pair range
		 * \param[in] locusIndexes vector of locus indexes
		 * \param[in] similarityCutOff only save pairs with at least this similarity
		 * \return `SimilarityMatrix` object with compressed indexed similarity values
		 */
		[[gnu::warn_unused_result]] SimilarityMatrix hashJacBlock_(const std::pair<RowColIdx, RowColIdx> &blockRange, const std::vector<uint32_t> &locusIndexes, const float &similarityCutOff) const;
		/** \brief Hash-based similarity in a range of locus groups
		 *
		 * The range points to a hash table of locus indexes.
		 * Start and end can be within a group.
		 *
		 * \param[in] blockRange range of locus hash table groups
		 * \param[in] similarityCutOff only save pairs with at least this similarity
		 * \return a `SimilarityMatrix` object
		 */
		[[gnu::warn_unused_result]] SimilarityMatrix hashJacBlock_(const std::pair<HashGroupItPairCount, HashGroupItPairCount> &blockRange, const float &similarityCutOff) const;
		/** \brief Hash-based similarity between locus pairs using multiple threads
		 *
		 * The provided row and column values index the `locusIndexes` vector that translates them to the actual locus indexes.
		 * This is necessary for hash group processing.
		 *
		 * \param[in] indexPairs vector of row/column index ranges, one per thread
		 * \param[in] locusIndexes vector of locus indexes
		 * \param[in] similarityCutOff only save pairs with at least this similarity
		 * \return `SimilarityMatrix` object with compressed indexed similarity values
		 */
		[[gnu::warn_unused_result]] SimilarityMatrix hashJacThreaded_(const std::vector< std::pair<RowColIdx, RowColIdx> > &indexPairs,
				const std::vector<uint32_t> &locusIndexes, const float &similarityCutOff) const;
		/** \brief Threaded hash-based similarity in ranges of locus groups
		 *
		 * The ranges point to a hash table of locus indexes.
		 * Starts and ends can be within a group.
		 *
		 * \param[in] blockRanges ranges of locus hash table groups
		 * \param[in] similarityCutOff only save pairs with at least this similarity
		 * \return a `SimilarityMatrix` object
		 */
		[[gnu::warn_unused_result]] SimilarityMatrix hashJacThreaded_(const std::vector< std::pair<HashGroupItPairCount, HashGroupItPairCount> > &blockRanges, const float &similarityCutOff) const;
		/** \brief Calculate the union and intersection Jaccard similarity pair
		 *
		 * \param[in] rowColumn indexes of the locus pair to compare
		 * \return Intersection and union for Jaccard similarity calculation between two loci using OPH
		 */
		[[gnu::warn_unused_result]] JaccardPair makeJaccardPair_(const RowColIdx &rowColumn) const noexcept;
	};
}
