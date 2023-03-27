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
#include <array>
#include <unordered_map>
#include <utility>        // for std::pair
#include <string>
#include <thread>
#include <mutex>

#include "random.hpp"

namespace BayesicSpace {
	struct IndexedPairSimilarity; 
	class GenoTableBin;
	class GenoTableHash;

	/** \brief Jaccard value with indexes
	 *
	 * Groups a Jaccard similarity value of two elements (e.g., loci or individuals) with their indexes and the ID of the group they belong to.
	 */
	struct IndexedPairSimilarity {
		float similarityValue;
		size_t element1ind;
		size_t element2ind;
		uint32_t groupID;
	};

	/** \brief LD value with indexes
	 *
	 * Groups a Jaccard similarity value and the D linkage disequilibrium statistic of two loci with their indexes.
	 */
	struct IndexedPairLD {
		float jaccard;
		float d;
		size_t element1ind;
		size_t element2ind;
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
		GenoTableBin(const std::string &inputFileName, const size_t &nIndividuals, std::string logFileName) : GenoTableBin( inputFileName, nIndividuals, std::move(logFileName), std::thread::hardware_concurrency() ){};
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
		GenoTableBin(const std::string &inputFileName, const size_t &nIndividuals, std::string logFileName, const size_t &nThreads);
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
		GenoTableBin(const std::vector<int> &maCounts, const size_t &nIndividuals, std::string logFileName) : GenoTableBin( maCounts, nIndividuals, std::move(logFileName), std::thread::hardware_concurrency() ){};
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
		GenoTableBin(const std::vector<int> &maCounts, const size_t &nIndividuals, std::string logFileName, const size_t &nThreads);

		/** \brief Copy constructor (deleted) */
		GenoTableBin(const GenoTableBin &toCopy) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTableBin operator=(const GenoTableBin &toCopy) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		GenoTableBin(GenoTableBin &&toMove) noexcept;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to be moved
		 * \return `GenoTableBin` object
		 */
		GenoTableBin& operator=(GenoTableBin &&toMove) noexcept;
		/** \brief Destructor */
		~GenoTableBin() = default;

		/** \brief Save the binary genotype file
		 *
		 * Saves the binary approximate genotype data to a binary file.
		 *
		 * \param[in] outFileName output file name
		 */
		void saveGenoBinary(const std::string &outFileName) const;
		/** \brief All by all Jaccard similarity LD
		 *
		 * Calculates linkage disequilibrium among all loci using Jaccard similarity as the statistic.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * All values belong to the same group. Row and column (1-base) indexes of the similarity matrix are also included in the tab-delimited output file.
		 * The lower triangle is vectorized by column (i.e. all correlations of the first locus, then all remaining correlations of the second, etc.).
		 *
		 * \param[in] ldFileName name of the output file
		 * \return lower triangle of the LD matrix
		 */
		void allJaccardLD(const std::string &ldFileName) const;
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
		size_t nIndividuals_;
		/** \brief Number of loci */
		size_t nLoci_;
		/** \brief Binarized locus size in bytes */
		size_t binLocusSize_;
		/** \brief Maximal number of threads to use */
		size_t nThreads_;
		/** \brief Random number generator */
		RanDraw rng_;
		/** \brief The mutex */
		mutable std::mutex mtx_;
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
		/** \brief Maximum number of loci for all by all LD */
		static const size_t maxNlocusPairs_;
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
		 * \param[in] firstLocusInd overall index of the first locus in the range
		 * \param[in] bedLocusLength number of bytes in each locus
		 */
		void bed2binBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const size_t &firstLocusInd, const size_t &bedLocusLength);
		/** \brief Multi-threaded binarization of _.bed_ file input 
		 *
		 * Binarizes loci from a _.bed_ file using multiple threads.
		 *
		 * \param[in] bedData vector of _.bed_ file input
		 * \param[in] threadRanges vector of locus index ranges, one per thread
		 * \param[in] firstLocusInd overall index of the first locus in the range
		 * \param[in] bedLocusLength number of bytes in each locus
		 * \return new value of `firstLocusInd`
		 */
		size_t bed2binThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const size_t &firstLocusInd, const size_t &bedLocusLength);
		/** \brief Binarize minor allele counts in a locus block
		 *
		 * Binarizes a portion of a vector of per-individual minor allele counts (0, 1, or 2; see the count vector constructor documentation for details).
		 *
		 * \param[in] macData vector of minor allele counts
		 * \param[in] locusIndRange locus index range
		 * \param[in] randVecLen length of the random bit vector (for heterozygote resolution)
		 */
		void mac2binBlk_(const std::vector<int> &macData, const std::pair<size_t, size_t> &locusIndRange, const size_t &randVecLen);
		/** \brief Jaccard similarity in a block of loci
		 *
		 * \param[in] blockVecRange block index range in `ldVec`
		 * \param[in] blockStartAll index of the block start in the overall vectorized LD matrix
		 * \param[out] ldVec vectorized lower triangle of the Jaccard and D similarity matrix with locus pair indexes
		 */
		void jaccardBlock_(const std::pair<size_t, size_t> &blockVecRange, const size_t &blockStartAll, std::vector<IndexedPairLD> &ldVec) const;
		/** \brief Jaccard similarity between locus pairs using multiple threads
		 *
		 * \param[in] pairIndRanges vector of pair ranges, one range per thread
		 * \param[in] blockStartAll index of the block start in the overall vectorized LD matrix
		 * \param[out] ldVec vectorized lower triangle of the Jaccard and D similarity matrix with locus pair indexes
		 * \return new block start index
		 */
		size_t jaccardThreaded_(const std::vector< std::pair<size_t, size_t> > &pairIndRanges, const size_t &blockStartAll, std::vector<IndexedPairLD> &ldVec) const;
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
		GenoTableHash() : nIndividuals_{0}, kSketches_{0}, fSketches_{0.0}, sketchSize_{0}, nLoci_{0}, locusSize_{0}, nFullWordBytes_{0}, nThreads_{1} {};
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
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] kSketches number of sketches per locus
		 * \param[in] nThreds maximal number of threads to use
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::string &inputFileName, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, std::string logFileName);
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
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] kSketches the number of sketches per locus
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::string &inputFileName, const size_t &nIndividuals, const size_t &kSketches, std::string logFileName) :
							GenoTableHash( inputFileName, nIndividuals, kSketches, std::thread::hardware_concurrency(), std::move(logFileName) ) {};
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
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] kSketches the number of sketches per locus
		 * \param[in] nThreds maximal number of threads to use
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, std::string logFileName);
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
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] kSketches the number of sketches per locus
		 * \param[in] logFileName name of the log file
		 */
		GenoTableHash(const std::vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches, std::string logFileName) :
				GenoTableHash( maCounts, nIndividuals, kSketches, std::thread::hardware_concurrency(), std::move(logFileName) ) {};

		/** \brief Copy constructor (deleted) */
		GenoTableHash(const GenoTableHash &toCopy) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTableHash operator=(const GenoTableHash &toCopy) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] toMove object to move
		 */
		GenoTableHash(GenoTableHash &&toMove) noexcept;
		/** \brief Move assignment operator
		 *
		 * \param[in] toMove object to be moved
		 * \return `GenoTableHash` object
		 */
		GenoTableHash& operator=(GenoTableHash &&toMove) noexcept;
		/** \brief Destructor */
		~GenoTableHash() = default;

		/** \brief All by all LD from hashes
		 *
		 * Calculates linkage disequilibrium among all loci using a modified OPH.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * All values belong to the same group. Row and column indexes (1-base) of the similarity matrix are also included in the tab-delimited output file.
		 * The lower triangle is vectorized by column (i.e. all correlations of the first locus, then all remaining correlations of the second, etc.).
		 *
		 * \param[in] ldFileName name of file to save the results
		 */
		void allHashLD(const std::string &ldFileName) const;
		/** \brief Assign groups by linkage disequilibrium (LD)
		 *
		 * The sketch matrix is divided into bands, `nRowsPerBand` rows per band (must be 1 or greater).
		 * Locus pairs are included in the pair hash table if all rows in at least one band match.
		 * The resulting hash table has groups with at least two loci per group (indexed by a hash of OPH sketches within bands).
		 * Locus indexes are in increasing order within each group.
		 * Some locus pairs may end up in more than one group, but no groups are completely identical in locus composition.
		 *
		 * \param[in] nRowsPerBand number of rows per sketch matrix band
		 *
		 * \return locus index hash table
		 */
		std::unordered_map< uint32_t, std::vector<size_t> > makeLDgroups(const size_t &nRowsPerBand) const;
		/** \brief Linkage disequilibrium (LD) in groups
		 *
		 * Group loci according to LD using the algorithm for `makeLDgroups` and calculate similarity within  groups.
		 *
		 * \param[in] nRowsPerBand number of rows per sketch matrix band
		 * \param[in] outFileName name of the output file
		 */
		void ldInGroups(const size_t &nRowsPerBand, const std::string &outFileName) const;
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
		size_t nIndividuals_;
		/** \brief Number of sketches */
		size_t kSketches_;
		/** \brief Number of sketches, float representation */
		float fSketches_;
		/** \brief Sketch size */
		size_t sketchSize_;
		/** \brief Number of loci */
		size_t nLoci_;
		/** \brief Locus size in bytes */
		size_t locusSize_;
		/** \brief Number of bytes in divisible by `llWordSize_` */
		size_t nFullWordBytes_;
		/** \brief Maximal number of threads to use */
		size_t nThreads_;
		/** \brief Random number generator */
		mutable RanDraw rng_;
		/** \brief The mutex */
		mutable std::mutex mtx_;
		/** \brief Log messages */
		mutable std::string logMessages_;
		/** \brief Log file name */
		std::string logFileName_;
		/** \brief Maximum number that does not overflow a triangle of an all by all comparison matrix */
		static const size_t maxPairs_;
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
		static const uint64_t roundMask_;
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
		 * The `seeds` vector may be appended by the function if additional seeds are required.
		 *
		 * \param[in] locusInd locus index
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in,out] seeds random number seeds for empty bin filling
		 * \param[in,out] binLocus vector of genotypes for a locus
		 */
		void locusOPH_(const size_t &locusInd, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds, std::vector<uint8_t> &binLocus);
		/** \brief OPH from _.bed_ file input
		 *
		 * Hashes a portion of a vector of input from a _.bed_ file that corresponds to a range of loci.
		 *
		 * \param[in] bedData _.bed_ file input
		 * \param[in] bedLocusIndRange range of locus indexes in the block
		 * \param[in] firstLocusInd overall index of the first locus in the range
		 * \param[in] locusLength number of bytes in each locus
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in] padIndiv additional individuals, `first` is the placement index, `second` is the index of the individual to add
		 * \param[in,out] seeds random number seeds for empty bin filling
		 */
		void bed2ophBlk_(const std::vector<char> &bedData, const std::pair<size_t, size_t> &bedLocusIndRange, const size_t &firstLocusInd,
							const size_t &bedLocusLength, const std::vector<size_t> &permutation, std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds);
		/** \brief OPH from _.bed_ file input using multiple threads
		 *
		 * Hashes input from a _.bed_ file using multiple threads.
		 *
		 * \param[in] bedData _.bed_ file input
		 * \param[in] threadRanges vector of locus ranges, one per thread
		 * \param[in] firstLocusInd overall index of the first locus in the range
		 * \param[in] locusLength number of bytes in each locus
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in] padIndiv additional individuals, `first` is the placement index, `second` is the index of the individual to add
		 * \param[in,out] seeds random number seeds for empty bin filling
		 * \return new value of `firstLocusInd`
		 */
		size_t bed2ophThreaded_(const std::vector<char> &bedData, const std::vector< std::pair<size_t, size_t> > &threadRanges, const size_t &firstLocusInd,
							const size_t &bedLocusLength, const std::vector<size_t> &permutation, std::vector< std::pair<size_t, size_t> > &padIndiv, std::vector<uint32_t> &seeds);
		/** \brief OPH from minor allele counts
		 *
		 * Hashes a portion of a vector of per-individual minor allele counts (0, 1, or 2; see the count vector constructor documentation for details).
		 * The vector portion corresponds to a block of loci.
		 *
		 * \param[in] macData vector of minor allele counts
		 * \param[in] startLocusInd index of the first locus in block
		 * \param[in] endLocusInd index of one past the last locus in block
		 * \param[in] randVecLen length of the random bit vector (for heterozygote resolution)
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in,out] seeds random number seeds for empty bin filling
		 */
		void mac2ophBlk_(const std::vector<int> &macData, const size_t &startLocusInd, const size_t &endLocusInd, const size_t &randVecLen, const std::vector<size_t> &permutation, std::vector<uint32_t> &seeds);
		/** \brief Hash-based similarity in a block of loci
		 *
		 * Pairwise hash-estimated Jaccard similarity among loci in a block continuous in a vectorized lower triangle of similarity values.
		 * The range of indexes refers to a vectorized by column lower triangle of a similarity matrix.
		 *
		 * \param[in] blockRange index range in a block in `hashJacVec`
		 * \param[in] blockStartAll index of the block start in the overall vectorized LD matrix
		 * \param[out] hashJacVec vector of similarity values with associated locus indexes and group IDs
		 */
		void hashJacBlock_(const std::pair<size_t, size_t> &blockRange, const size_t &blockStartAll, std::vector<IndexedPairSimilarity> &hashJacVec) const;
		/** \brief Hash-based similarity among indexed loci
		 *
		 * Pairwise hash-estimated Jaccard similarities among loci indexed by the provided vector. This is for blocked estimates.
		 * The index range refers to the portion of the vectorized by column lower triangle of the resulting block similarity matrix.
		 * The index vector contains indexes of locus pairs included in LD calculations.
		 *
		 * \param[in] blockStartVec index of the block start
		 * \param[in] blockEndVec index of the block end
		 * \param[in, out] indexedJacc vector of similarity values with associated locus indexes and group IDs
		 */
		void hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, std::vector< IndexedPairSimilarity > &indexedJacc) const;
		/** \brief Hash-based similarity using multiple threads
		 *
		 * Pairwise hash-estimated Jaccard similarity among loci in a block continuous in a vectorized lower triangle of similarity values using multiple threads.
		 * The ranges of indexes refer to a vectorized by column lower triangle of a similarity matrix.
		 *
		 * \param[in] threadRanges vector of block ranges, one per tread, in `hashJacVec`
		 * \param[in] blockStartAll index of the block start in the overall vectorized LD matrix
		 * \param[out] hashJacVec vector of similarity values with associated locus indexes and group IDs
		 * \return new block start index
		 */
		size_t hashJacThreaded_(const std::vector< std::pair<size_t, size_t> > &threadRanges, const size_t &blockStartAll, std::vector<IndexedPairSimilarity> &hashJacVec) const;
	};
}

