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

#ifndef gvhash_hpp
#define gvhash_hpp

#include <cstddef>
#include <vector>
#include <array>
#include <utility>  // for std::pair
#include <string>
#include <thread>
#include <mutex>

#include "bayesicUtilities/random.hpp"

using std::vector;
using std::array;
using std::pair;
using std::string;
using std::thread;
using std::mutex;

namespace BayesicSpace {
	class GenoTableBinCPP;
	class GenoTableBin;
	class GenoTableHash;

	/** \brief Count set bits in a 16-bit word
	 *
	 * Counting the set bits using Karnigan's method. Passing by value to modify the copy and also because the address is much bigger than 16 bits.
	 *
	 * \param[in] inVal input value
	 * \return number of bits set
	 */
	uint16_t countSetBits(uint16_t inVal);
	/** \brief Count set bits in a vector
	 *
	 * Counting the set bits in a vector of bytes using Karnigan's method.
	 *
	 * \param[in] inVec input vector
	 * \return number of bits set
	 */
	uint32_t countSetBits(const vector<uint8_t> &inVec);
	/** \brief Count set bits in a range within a vector
	 *
	 * Counting the set bits in a range within a vector of bytes using Karnigan's method.
	 *
	 * \param[in] inVec input vector
	 * \param[in] start staring index
	 * \param[in] length number of bytes to process
	 * \return number of bits set
	 */
	uint32_t countSetBits(const vector<uint8_t> &inVec, const size_t &start, const size_t &length);
	/** \brief Get available RAM
	 *
	 * Estimates available RAM. If `procfs` is mounted, uses information from there. Otherwise, sets available RAM to 2 GiB.
	 *
	 * \return estimated available RAM in bytes
	 */
	size_t getAvailableRAM();

	/** \brief Class to store binary compressed genotype tables
	 *
	 * Converts genotype data to a lossy compressed binary code.
	 * Genotypes are stored in memory in a one-bit format: bit set for the minor allele, unset for the major.
	 * Bits corresponding to missing data are unset (this is the same as mean imputation), heterozygotes are set with a 50% probability.
	 */
	class GenoTableBin {
	public:
		/** \brief Default constructor */
		GenoTableBin(){};
		/** \brief Constructor with input file name
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed).
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * If necessary, alleles are re-coded so that the set bit is always the minor allele.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] nIndividuals number of genotyped individuals
		 */
		GenoTableBin(const string &inputFileName, const size_t &nIndividuals) : GenoTableBin( inputFileName, nIndividuals, thread::hardware_concurrency() ){};
		/** \brief Constructor with input file name and thread count
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed).
		 * Heterozygotes are assigned the major or minor allele at random, missing genotypes are assigned the major allele.
		 * If necessary, alleles are re-coded so that the set bit is always the minor allele.
		 * The number of threads requested is maximum to be used, depending on available system resources.
		 *
		 * \param[in] inputFileName input file name
		 * \param[in] nIndividuals number of genotyped individuals
		 * \param[in] nThreads maximal number of threads to use
		 */
		GenoTableBin(const string &inputFileName, const size_t &nIndividuals, const size_t &nThreads);
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
		 */
		GenoTableBin(const vector<int> &maCounts, const size_t &nIndividuals) : GenoTableBin( maCounts, nIndividuals, thread::hardware_concurrency() ){};
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
		 * \param[in] nThreads maximal number of threads to use
		 */
		GenoTableBin(const vector<int> &maCounts, const size_t &nIndividuals, const size_t &nThreads);

		/** \brief Copy constructor (deleted) */
		GenoTableBin(const GenoTableBin &in) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTableBin operator=(const GenoTableBin &in) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] in object to move
		 */
		GenoTableBin(GenoTableBin &&in) noexcept;
		/** \brief Move assignment operator
		 *
		 * \param[in] in object to be moved
		 * \return `GenoTableBin` object
		 */
		GenoTableBin& operator=(GenoTableBin &&in) noexcept;

		/** \brief Save the binary genotype file
		 *
		 * Saves the binary approximate genotype data to a binary file.
		 *
		 * \param[in] outFileName output file name
		 */
		void saveGenoBinary(const string &outFileName) const;
		/** \brief All by all Jaccad similarity LD
		 *
		 * Calculates linkage disequilibrium among all loci using a corrected Jaccard similarity as the statistic.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * The lower triangle is vectorized by column (i.e. all correlations of the first locus, then all remaining correlations of the second, etc.).
		 *
		 * \param[in] ldFileName name of the output file
		 * \return lower triangle of the LD matrix
		 */
		void allJaccardLD(const string &ldFileName) const;
	protected:
		/** \brief Binarized genotype table
		 *
		 * Stores one bit per genotype. Heterozygotes are randomly assigned, missing data are assigned 0.
		 */
		vector<uint8_t> binGenotypes_;
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
		mutable mutex mtx_;
		/** \brief Leading bytes for .bed files */
		static const array<char, 3> magicBytes_;
		/** \brief One set bit for masking */
		static const uint8_t oneBit_;
		/** \brief Size of one byte in bits */
		static const uint8_t byteSize_;
		/** \brief 64 bit word size in bytes */
		static const uint8_t llWordSize_;
		/** \brief Maximum number of loci for all by all LD */
		static const size_t maxNlocusPairs_;
		/** \brief Binarize a range of loci from _.bed_ file input
		 *
		 * Binarizes a range of loci from a vector of input from a _.bed_ file.
		 *
		 * \param[in] bedData vector of _.bed_ file input
		 * \param[in] firstBedLocusInd index of the first locus in the _.bed_ vector
		 * \param[in] lastBedLocusInd index of one past the last locus in the _.bed_ vector
		 * \param[in] firstLocusInd overall index of the first locus in the range
		 * \param[in] bedLocusLength number of bytes in each locus
		 * \param[in] randVecLen length of the random bit vector (for heterozygote resolution)
		 */
		void bed2binBlk_(const vector<char> &bedData, const size_t &firstBedLocusInd, const size_t &lastBedLocusInd, const size_t &firstLocusInd, const size_t &bedLocusLength, const size_t &randVecLen);
		/** \brief Binarize minor allele counts in a locus block
		 *
		 * Binarizes a portion of a vector of per-individual minor allele counts (0, 1, or 2; see the count vector constructor documentation for details).
		 *
		 * \param[in] macData vector of minor allele counts
		 * \param[in] startLocusInd first locus index
		 * \param[in] endLocusInd one past the last locus index
		 * \param[in] randVecLen length of the random bit vector (for heterozygote resolution)
		 */
		void mac2binBlk_(const vector<int> &macData, const size_t &startLocusInd, const size_t &endLocusInd, const size_t &randVecLen);
		/** \brief Jaccard similarity in a block of loci
		 *
		 * \param[in] blockStartVec index of the block start in `jaccardVec`
		 * \param[in] blockEndVec index of one past the block end in `jaccardVec`
		 * \param[in] blockStartAll index of the block start in the overall vectorized LD matrix
		 * \param[out] jaccardVec vectorized lower triangle of the Jaccard similarity matrix
		 */
		void jaccardBlock_(const size_t &blockStartVec, const size_t &blockEndVec, const size_t &blockStartAll, vector<float> &jaccardVec) const;
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
		GenoTableHash() : kSketches_{0} {};
		/** \brief Constructor with input file name and thread number
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed).
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
		GenoTableHash(const string &inputFileName, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, const string &logFileName);
		/** \brief Constructor with input file name
		 *
		 * The file should be in the `plink` [.bed format](https://www.cog-genomics.org/plink/1.9/formats#bed).
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
		GenoTableHash(const string &inputFileName, const size_t &nIndividuals, const size_t &kSketches, const string &logFileName) : GenoTableHash(inputFileName, nIndividuals, kSketches, thread::hardware_concurrency(), logFileName) {};
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
		GenoTableHash(const vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches, const size_t &nThreads, const string &logFileName);
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
		GenoTableHash(const vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches, const string &logFileName) : GenoTableHash(maCounts, nIndividuals, kSketches, thread::hardware_concurrency(), logFileName) {};

		/** \brief Copy constructor (deleted) */
		GenoTableHash(const GenoTableHash &in) = delete;
		/** \brief Copy assignment operator (deleted) */
		GenoTableHash operator=(const GenoTableHash &in) = delete;
		/** \brief Move constructor
		 *
		 * \param[in] in object to move
		 */
		GenoTableHash(GenoTableHash &&in) noexcept;
		/** \brief Move assignment operator
		 *
		 * \param[in] in object to be moved
		 * \return `GenoTableHash object
		 */
		GenoTableHash& operator=(GenoTableHash &&in) noexcept;

		/** \brief All by all LD from hashes
		 *
		 * Calculates linkage disequilibrium among all loci using a modified OPH.
		 * Result is a vectorized lower triangle of the symmetric \f$N \times N\f$ similarity matrix, where \f$N\f$ is the number of loci.
		 * The lower triangle is vectorized by column (i.e. all correlations of the first locus, then all remaining correlations of the second, etc.).
		 *
		 * \param[in] ldFileName name of file to save the results
		 */
		void allHashLD(const string &ldFileName) const;
		/** \brief Assign groups by local linkage disequilibrium (LD)
		 *
		 * Group loci by LD along the genome. The algorithm is
		 * Start by using simHash on the first `kSketchSubset` of the first locus OPH. Proceed along the genome, for each locus
		 * - simHash `kSketchSubset` of the OPH
		 * - compare to the latest group simHash
		 * - if the Hamming distance from the latest group is less than `hammingCutoff`, add the locus index to the group
		 * - if not, compare to up to `lookBackNumber` of groups back along the genome, adding to the first group that meets the cut-off
		 * - if none of the previous groups are close enough, start a new group, labeling it with the current simHash.
		 *
		 * \param[in] hammingCutoff the maximum Hamming distance for group inclusion
		 * \param[in] kSketchSubset number of OPH sketches to use for simHash
		 * \param[in] lookBackNumber number of previous groups to consider
		 *
		 * \return group IDs for each locus
		 */
		vector< vector<size_t> > makeLDgroups(const uint16_t &hammingCutoff, const size_t &kSketchSubset, const size_t &lookBackNumber) const;
		/** \brief Calculates linkage disequilibrium (LD) in local groups
		 *
		 * Group loci according to LD using the algorithm for `makeLDgroups` and calculate similarity within  groups.
		 * All hash values are used for similarity calculations, even if only a subset is considered for similarity grouping.
		 *
		 * \param[in] hammingCutoff the maximum Hamming distance for group inclusion
		 * \param[in] kSketchSubset number of OPH sketches to use for simHash
		 * \param[in] lookBackNumber number of previous groups to consider
		 * \param[in] smallestGrpSize groups with fewer loci than this will be discarded from LD calculations
		 * \param[in] outFileName name of the output file
		 */
		void ldInGroups(const uint16_t &hammingCutoff, const size_t &kSketchSubset, const size_t &lookBackNumber, const size_t &smallestGrpSize, const string &outFileName) const;
		/** \brief Save the log to a file
		 *
		 * Log file name provided at construction.
		 */
		void saveLogFile() const;
	protected:
		/** \brief Vector of sketches
		 *
		 * A sketch is the position of the first set bit in a bin of permuted bits.
		 */
		vector<uint16_t> sketches_;
		/** \brief Number of individuals */
		size_t nIndividuals_;
		/** \brief Number of sketches */
		size_t kSketches_;
		/** \brief Sketch size */
		size_t sketchSize_;
		/** \brief Number of loci */
		size_t nLoci_;
		/** \brief Locus size in bytes */
		size_t locusSize_;
		/** \brief Maximal number of threads to use */
		size_t nThreads_;
		/** \brief Random number generator */
		RanDraw rng_;
		/** \brief The mutex */
		mutable mutex mtx_;
		/** \brief Log messages */
		mutable string logMessages_;
		/** \brief Log file name */
		string logFileName_;
		/** \brief Maximum number that does not overflow a triangle of an all by all comparison matrix */
		static const size_t maxPairs_;
		/** \brief Leading bytes for .bed files */
		static const array<char, 3> magicBytes_;
		/** \brief One set bit for masking */
		static const uint8_t oneBit_;
		/** \brief Size of one byte in bits */
		static const uint8_t byteSize_;
		/** \brief 64 bit word size in bytes */
		static const uint8_t llWordSize_;
		/** \brief MurMurHash number of blocks */
		static const size_t nblocks_;
		/** \brief MurMurHash key length */
		static const uint32_t mmhKeyLen_;
		/** \brief Value corresponding to an empty token */
		static const uint16_t emptyBinToken_;
		/** \brief MurMurHash c1 constant */
		static const uint32_t c1_;
		/** \brief MurMurHash c2 constant */
		static const uint32_t c2_;
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
		void locusOPH_(const size_t &locusInd, const vector<size_t> &permutation, vector<uint32_t> &seeds, vector<uint8_t> &binLocus);
		/** \brief OPH from _.bed_ file input
		 *
		 * Hashes a portion of a vector of input from a _.bed_ file that corresponds to a range of loci.
		 *
		 * \param[in] bedData _.bed_ file input
		 * \param[in] firstBedLocusInd index of the first locus in the _.bed_ vector
		 * \param[in] lastBedLocusInd index of one past the last locus in the _.bed_ vector
		 * \param[in] firstLocusInd overall index of the first locus in the range
		 * \param[in] locusLength number of bytes in each locus
		 * \param[in] randVecLen length of the random bit vector (for heterozygote resolution)
		 * \param[in] permutation permutation to be applied to each locus 
		 * \param[in,out] seeds random number seeds for empty bin filling
		 */
		void bed2ophBlk_(const vector<char> &bedData, const size_t &firstBedLocusInd, const size_t &lastBedLocusInd, const size_t &firstLocusInd,
							const size_t &bedLocusLength, const size_t &randVecLen, const vector<size_t> &permutation, vector<uint32_t> &seeds);
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
		void mac2ophBlk_(const vector<int> &macData, const size_t &startLocusInd, const size_t &endLocusInd, const size_t &randVecLen, const vector<size_t> &permutation, vector<uint32_t> &seeds);
		/** \brief MurMurHash to fill in empty bins
		 *
		 * Generates a 32-bit hash of an index value using the MurMurHash3 algorithm.
		 *
		 * \param[in] key the key to be hashed
		 * \param[in] seed the seed
		 *
		 * \return the hash value
		 */
		uint32_t murMurHash_(const size_t &key, const uint32_t &seed) const;
		/** \brief 16 bit MurMurHash
		 *
		 * Generates a 16-bit hash of a 16 bit value using the MurMurHash3 algorithm.
		 *
		 * \param[in] key the key to be hashed
		 * \param[in] seed the seed
		 *
		 * \return the hash value
		 */
		uint16_t murMurHash_(const uint16_t &key, const uint32_t &seed) const;
		/** \brief MurMurHash of a sketch portion
		 *
		 * Generates a 32-bit hash of a portion of an OPH sketch using the MurMurHash3 algorithm.
		 * If the number of elements provided is odd, it is rounded down to the next even number.
		 *
		 * \param[in] startInd index of the first sketch element (to the `sketches_` vector)
		 * \param[in] nElements number of elements to hash
		 * \param[in] seed the seed
		 *
		 * \return the hash value
		 */
		uint32_t murMurHash_(const size_t &startInd, const size_t &nElements, const uint32_t &seed) const;
		/** \brief SimHash of OPH sketches
		 *
		 * Takes the collection of OPH sketches for each locus and returns a simHash.
		 *
		 * \param[in] startInd index of the first sketch in the `sketches_` vector
		 * \param[in] kSketches number of sketches to hash
		 * \param[in] seed seed value to use with murMurHash on each sketch element
		 * \return hash value
		 */
		uint16_t simHash_(const size_t &startInd, const size_t &kSketches, const uint32_t &seed) const;
		/** \brief Hash-based similarity in a block of loci
		 *
		 * Pairwise hash-estimated Jaccard similarity among loci in a block continuous in a vectorized lower triangle of similarity values.
		 * The range of indexes refers to a vectorized by column lower triangle of a similarity matrix.
		 *
		 * \param[in] blockStartVec index of the block start in `hashJacVec`
		 * \param[in] blockEndVec index of one past the block end in `hashJacVec`
		 * \param[in] blockStartAll index of the block start in the overall vectorized LD matrix
		 * \param[out] hashJacVec vectorized lower triangle of the hash-estimated Jaccard similarity matrix
		 */
		void hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, const size_t &blockStartAll, vector<float> &hashJacVec) const;
		/** \brief Hash-based similarity among indexed loci
		 *
		 * Pairwise hash-estimated Jaccard similarities among loci indexed by the provided vector. This is for blocked estimates.
		 * The index range refers to the portion of the vectorized by column lower triangle of the resulting block similarity matrix.
		 * The index vector contains indexes of locus pairs included in LD calculations.
		 *
		 * \param[in] blockStartVec index of the block start in `idxVector` and `hashJacVec`
		 * \param[in] blockEndVec index of the block end in `idxVector` and `hashJacVec`
		 * \param[in] idxVector vector of locus pair indexes
		 * \param[out] hashJacVec vectorized blocked lower triangle of the hash-estimated Jaccard similarity matrix
		 */
		void hashJacBlock_(const size_t &blockStartVec, const size_t &blockEndVec, const vector< pair<size_t, size_t> > &idxVector, vector<float> &hashJacVec) const;
		/** \brief Hamming distance
		 *
		 * Calculates the bit-wise Hamming distance between two 16-bit variables. Passing the variables by value since they are much smaller than addresses.
		 *
		 * \param[in] first first value
		 * \param[in] second second value
		 * \return Hamming distance
		 *
		 */
		uint16_t inline hammingDistance_(uint16_t first, uint16_t second) const {return countSetBits(first ^ second); }
	};
}

#endif // gvhash_hpp
