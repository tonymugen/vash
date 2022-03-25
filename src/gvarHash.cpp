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
 * \version 0.1
 *
 * Implementation of classes that take binary variant files and generate lossy summaries with hashing.
 *
 */

#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <iterator>
#include <algorithm>
#include <limits>
#include <fstream>
#include <future>
#include <mutex>

#include <chrono>
#include <cmath>

#include "gvarHash.hpp"

using std::vector;
using std::array;
using std::distance;
using std::string;
using std::to_string;
using std::move;
using std::numeric_limits;
using std::fstream;
using std::ios;
using std::streampos;
using std::future;
using std::async;
using std::mutex;
using std::lock_guard;

using std::chrono::high_resolution_clock;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::milli;

using namespace BayesicSpace;
// External functions
uint16_t BayesicSpace::countSetBits(uint16_t inVal) {
	uint16_t totSet = 0;
	for (; inVal; ++totSet) {
		inVal &= inVal - 1;
	}
	return totSet;
}

uint32_t BayesicSpace::countSetBits(const vector<uint8_t> &inVec) {
	uint32_t totSet = 0;
	for (const auto &in : inVec){
		uint8_t v = in;
		for (; v; ++totSet) {
			v &= v - 1;
		}
	}
	return totSet;
}

uint32_t BayesicSpace::countSetBits(const vector<uint8_t> &inVec, const size_t &start, const size_t &length) {
	uint32_t totSet = 0;
	for (size_t i = start; i < start + length; ++i){
		uint8_t v = inVec[i];
		for (; v; ++totSet) {
			v &= v - 1;
		}
	}
	return totSet;
}

// GenoTableBin methods
constexpr array<char, 3> GenoTableBin::magicBytes_ = {0x6c, 0x1b, 0x01};   // Leading bytes for .bed files 
constexpr uint8_t  GenoTableBin::oneBit_           = 0b00000001;           // One set bit for masking 
constexpr uint8_t  GenoTableBin::byteSize_         = 8;                    // Size of one byte in bits 
constexpr uint8_t  GenoTableBin::llWordSize_       = 8;                    // 64 bit word size in bytes 

// Constructors
GenoTableBin::GenoTableBin(const string &inputFileName, const size_t &nIndividuals) : nIndividuals_{nIndividuals}, nLoci_{0} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in ") + string(__FUNCTION__);
	} else if (nIndividuals > numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
		throw string("ERROR: the number of individuals (") + to_string(nIndividuals) + string( ") is too big to make a square relationship matrix in ") + string(__FUNCTION__);
	}
	size_t nBedBytes = nIndividuals_ / 4 + static_cast<bool>(nIndividuals_ % 4);
	fstream inStr;
	// Start by measuring file size
	inStr.open(inputFileName.c_str(), ios::in | ios::binary | ios::ate);
	const size_t N = static_cast<uint64_t>( inStr.tellg() ) - 3UL;       // first three bytes are majic
	inStr.close();
	nLoci_ = N / nBedBytes;
	inStr.open(inputFileName.c_str(), ios::in | ios::binary);
	char magicBuf[magicBytes_.size()]{};
	inStr.read( magicBuf, magicBytes_.size() );
	if ( inStr.eof() ){
		throw string("ERROR: No loci in the input .bed file ") + inputFileName + string(" in ") + string(__FUNCTION__);
	} else if (magicBuf[0] != magicBytes_[0]){
		throw string("ERROR: first magic byte in input .bed file is not the expected value in ") + string(__FUNCTION__);
	} else if (magicBuf[1] != magicBytes_[1]){
		throw string("ERROR: second magic byte in input .bed file is not the expected value in ") + string(__FUNCTION__);
	} else if (magicBuf[2] != magicBytes_[2]){
		throw string("ERROR: third magic byte in input .bed file is not the expected value in ") + string(__FUNCTION__);
	}
	// Generate the binary genotype table while reading the .bed file
	locusSize_ = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	binGenotypes_.resize(locusSize_ * nIndividuals_, 0);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	vector<char> bedLocus(nBedBytes, 0);
	// Create a vector to store random bytes for stochastic heterozygote resolution
	vector<uint64_t> rand( nBedBytes / llWordSize_ + static_cast<bool>(nBedBytes % llWordSize_) );
	uint8_t *randBytes   = reinterpret_cast<uint8_t*>( rand.data() );
	const float fNind    = static_cast<float>(nIndividuals_);
	const size_t endBed  = nBedBytes - 2UL + (nBedBytes & 1UL);
	const size_t addIndv = nIndividuals_ - endBed * 4;
	while ( inStr.read(bedLocus.data(), nBedBytes) ) {
		// Fill the random byte vector
		for (auto &rv : rand){
			rv = rng_.ranInt();
		}
		vector<uint8_t> binLocus(locusSize_, 0);
		vector<uint8_t> missMasks(locusSize_, 0);
		size_t iBinGeno = 0;                       // binGenotypes_ vector index
		// Two bytes of .bed code go into one byte of my binary representation
		// Therefore, work on two consecutive bytes of .bed code in the loop
		for (size_t iBed = 0; iBed < endBed ; iBed += 2){                     // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
			bedLocus[iBed]      = ~bedLocus[iBed];                            // flip so that homozygous second allele (usually minor) is set to 11
			uint8_t offsetToBin = 0;                                          // move the .bed mask by this much to align with the binarized byte
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedLocus[iBed] & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedLocus[iBed] & ( oneBit_ << (iInByteG + 1) );
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks[iBinGeno]      |= curMissMask >> offsetToBin;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask      |= randBytes[iBed] & (firstBitMask << 1);
				firstBitMask       &= secondBitMask >> 1;
				binLocus[iBinGeno] |= firstBitMask >> offsetToBin;
				++offsetToBin;
			}
			const size_t nextIbed = iBed + 1;
			bedLocus[nextIbed]    = ~bedLocus[nextIbed];
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedLocus[nextIbed] & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedLocus[nextIbed] & ( oneBit_ << (iInByteG + 1) );
				// Keep track of missing genotypes to revert them if I have to flip bits later on
				const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
				missMasks[iBinGeno]      |= curMissMask << offsetToBin;
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask      |= randBytes[nextIbed] & (firstBitMask << 1);
				firstBitMask       &= secondBitMask >> 1;
				binLocus[iBinGeno] |= firstBitMask << offsetToBin; // keep adding to the current binarized byte, so switch the direction of shift
				--offsetToBin;
			}
			++iBinGeno;
		}
		for (size_t iBed = endBed; iBed < nBedBytes; ++iBed){
			bedLocus[iBed] = ~bedLocus[iBed];
		}
		uint8_t inBedByteOffset = 0;
		for (size_t iInd = 0; iInd < addIndv; ++iInd){
			const size_t curBedByte = endBed + iInd / 4;
			uint8_t firstBitMask    = bedLocus[curBedByte] & (oneBit_ << inBedByteOffset);
			const uint8_t secondBBO = inBedByteOffset + 1;
			uint8_t secondBitMask   = bedLocus[curBedByte] & (oneBit_ << secondBBO);
			// Keep track of missing genotypes to revert them if I have to flip bits later on
			const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
			missMasks.back()         |= (curMissMask >> inBedByteOffset) << iInd;
			// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
			secondBitMask   |= randBytes[curBedByte] & (firstBitMask << 1);
			firstBitMask    &= secondBitMask >> 1;
			firstBitMask     = firstBitMask >> inBedByteOffset;
			firstBitMask     = firstBitMask << iInd;
			binLocus.back() |= firstBitMask;
			inBedByteOffset += 2;
			inBedByteOffset  = inBedByteOffset % 8;
		}
		float aaCount = static_cast<float>( countSetBits(binLocus) ) / fNind;
		if (aaCount > 0.5){ // always want the alternative to be the minor allele
			for (size_t iBL = 0; iBL < locusSize_; ++iBL){
				binGenotypes_.push_back( (~binLocus[iBL]) & (~missMasks[iBL]) );
			}
			binGenotypes_.back() &= lastByteMask; // unset the remainder bits
			aaCount = 1.0 - aaCount;
		} else {
			for (const auto &bg : binLocus){
				binGenotypes_.push_back(bg);
			}
		}
		aaf_.push_back(aaCount);
		++nLoci_;
	}
	inStr.close();
}

GenoTableBin::GenoTableBin(const vector<int> &maCounts, const size_t &nIndividuals) : nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in ") + string(__FUNCTION__);
	}
	if (maCounts.size() % nIndividuals){
		throw string("ERROR: length of allele count vector (") + to_string( maCounts.size() ) + string(" is not divisible by the provided number of individuals (") +
			to_string(nIndividuals) + string( ") in ") + string(__FUNCTION__);
	}
	if ( maCounts.empty() ){
		throw string("ERROR: empty vector of minor allele counts in ") + string(__FUNCTION__);
	}
	locusSize_                 = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	uint8_t remainderInd       = static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
	remainderInd               = byteSize_ - remainderInd;

	binGenotypes_.resize(nLoci_ * locusSize_, 0);

	// Create a vector to store random bytes for stochastic heterozygote resolution
	vector<uint64_t> rand( binGenotypes_.size() / llWordSize_ + static_cast<bool>(binGenotypes_.size() % llWordSize_) );
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	// Fill the random byte vector
	for (auto &rv : rand){
		rv = rng_.ranInt();
	}
	const float fNind = static_cast<float>(nIndividuals_);

	for (size_t jLoc = 0; jLoc < nLoci_; ++jLoc) {
		size_t iIndiv         = 0;
		const size_t begByte  = jLoc * locusSize_;
		const size_t begIndiv = jLoc * nIndividuals_;
		vector<uint8_t> missMasks(locusSize_, 0);
		size_t i0Byte = 0;                                                                         // to index the missMasks vector
		for (size_t iByte = begByte; iByte < begByte + locusSize_ - 1; ++iByte){                   // treat the last byte separately
			for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte){
				uint8_t curIndiv          = static_cast<uint8_t>(maCounts[begIndiv + iIndiv]);     // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= 0b10000011;                                            // mask everything in the middle
				const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
				missMasks[i0Byte]        |= (missingMask << iInByte);
				curIndiv                 &= 0b00000011;
				const uint8_t randMask    = (randBytes[iByte] >> iInByte) & oneBit_;               // 0b00000000 or 0b00000001 with equal chance
				uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
				curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
				binGenotypes_[iByte]     |= curBitMask << iInByte;
				++iIndiv;
			}
			++i0Byte;
		}
		// now deal with the last byte in the individual
		for (uint8_t iRem = 0; iRem < remainderInd; ++iRem){
			uint8_t curIndiv          = static_cast<uint8_t>(maCounts[begIndiv + iIndiv]);        // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= 0b10000011;                                               // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                            // 0b00000001 iff is missing (negative value)
			missMasks.back()         |= (missingMask << iRem);
			curIndiv                 &= 0b00000011;
			const uint8_t randMask    = (randBytes[begByte + locusSize_ - 1] >> iRem) & oneBit_;  // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                  // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                             // zero it out if missing value is set

			binGenotypes_[begByte + locusSize_ - 1] |= curBitMask << iRem;
			++iIndiv;
		}
		float maf = static_cast<float>( countSetBits(binGenotypes_, begByte, locusSize_) ) / fNind;
		if (maf > 0.5){ // always want the alternative to be the minor allele
			i0Byte = 0;
			for (size_t i = begByte; i < begByte + locusSize_; ++i){
				binGenotypes_[i] = (~binGenotypes_[i]) & (~missMasks[i0Byte]);
				++i0Byte;
			}
			binGenotypes_[begByte + locusSize_ - 1] &= lastByteMask; // unset the remainder bits
			maf = 1.0 - maf;
		}
		aaf_.push_back(maf);
	}
}

GenoTableBin::GenoTableBin(GenoTableBin &&in) noexcept {
	if (this != &in){
		binGenotypes_ = move(in.binGenotypes_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;
	}
}

GenoTableBin& GenoTableBin::operator=(GenoTableBin &&in) noexcept {
	if (this != &in){
		*this = move(in);
	}
	return *this;
}

void GenoTableBin::saveGenoBinary(const string &outFileName) const {
	fstream out;
	out.open(outFileName.c_str(), ios::out | ios::binary | ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), binGenotypes_.size() );
	out.close();
}

vector<float> GenoTableBin::allJaccardLD() const {
	// TODO: a try/catch block that writes to file directly if allocation fails
	vector<float> LDmat(nLoci_ * (nLoci_ - 1) / 2, 0.0);
	vector< future<void> > tasks;
	tasks.reserve(nLoci_);
	size_t blockBeg = 0; // index in LDmat of the first element of the block of similarities
	for (size_t iRow = 0; iRow < nLoci_; ++iRow) {
		tasks.emplace_back( async([this, iRow, blockBeg, &LDmat]{jaccardBlock_(iRow, blockBeg, LDmat);}) );
		blockBeg += nLoci_ - iRow - 1;
	}
	return LDmat;
}

void GenoTableBin::bed2bin_(const vector<char> &bedData, const size_t &bedBegInd, const size_t &binBegInd, const size_t &locusLength, const size_t &randVecLen){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	const size_t endWholeBed   = bedBegInd + locusLength - 2UL + (locusLength & 1UL);
	const size_t addIndv       = nIndividuals_ - endWholeBed * 4UL;
	const size_t lastBinInd    = binBegInd + locusSize_ - 1;
	// Fill the random byte vector
	vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	for (auto &rv : rand){
		rv = rng_.ranInt();
	}
	vector<uint8_t> missMasks(locusSize_, 0);
	size_t iBinGeno = binBegInd;                       // binLocus vector index
	uint8_t bedByte = 0;
	// Two bytes of .bed code go into one byte of my binary representation
	// Therefore, work on two consecutive bytes of .bed code in the loop
	for (size_t iBed = bedBegInd; iBed < endWholeBed ; iBed += 2){                // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
		bedByte             = ~bedData[iBed];                                  // flip so that homozygous second allele (usually minor) is set to 11
		uint8_t offsetToBin = 0;                                               // move the .bed mask by this much to align with the binarized byte
		for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
			uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
			uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
			// Keep track of missing genotypes to revert them if I have to flip bits later on
			const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
			missMasks[iBinGeno]      |= curMissMask >> offsetToBin;
			// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
			secondBitMask           |= randBytes[iBed] & (firstBitMask << 1);
			firstBitMask            &= secondBitMask >> 1;
			binGenotypes_[iBinGeno] |= firstBitMask >> offsetToBin;
			++offsetToBin;
		}
		const size_t nextIbed = iBed + 1;
		bedByte               = ~bedData[nextIbed];
		for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
			uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
			uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
			// Keep track of missing genotypes to revert them if I have to flip bits later on
			const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
			missMasks[iBinGeno]      |= curMissMask << offsetToBin;
			// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
			secondBitMask           |= randBytes[nextIbed] & (firstBitMask << 1);
			firstBitMask            &= secondBitMask >> 1;
			binGenotypes_[iBinGeno] |= firstBitMask << offsetToBin; // keep adding to the current binarized byte, so switch the direction of shift
			--offsetToBin;
		}
		++iBinGeno;
	}
	uint8_t inBedByteOffset = 0;
	for (size_t iInd = 0; iInd < addIndv; ++iInd){
		const size_t curBedByte = endWholeBed + iInd / 4;
		uint8_t firstBitMask    = (~bedData[curBedByte]) & (oneBit_ << inBedByteOffset);
		const uint8_t secondBBO = inBedByteOffset + 1;
		uint8_t secondBitMask   = (~bedData[curBedByte]) & (oneBit_ << secondBBO);
		// Keep track of missing genotypes to revert them if I have to flip bits later on
		const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
		missMasks.back()         |= (curMissMask >> inBedByteOffset) << iInd;
		// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
		secondBitMask             |= randBytes[curBedByte] & (firstBitMask << 1);
		firstBitMask              &= secondBitMask >> 1;
		firstBitMask               = firstBitMask >> inBedByteOffset;
		firstBitMask               = firstBitMask << iInd;
		binGenotypes_[lastBinInd] |= firstBitMask;
		inBedByteOffset           += 2;
		inBedByteOffset            = inBedByteOffset % 8;
	}
	float aaCount = static_cast<float>( countSetBits(binGenotypes_, binBegInd, locusSize_) ) / static_cast<float>(nIndividuals_);
	if (aaCount > 0.5){ // always want the alternative to be the minor allele
		for (size_t iBL = binBegInd; iBL < binBegInd + locusSize_; ++iBL){
			binGenotypes_[iBL] = (~binGenotypes_[iBL]) & (~missMasks[iBL]);
		}
		const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
		binGenotypes_[lastBinInd] &= lastByteMask; // unset the remainder bits
		aaCount = 1.0 - aaCount;
	}
	aaf_.push_back(aaCount);
}

void GenoTableBin::mac2bin_(const vector<int> &macData, const size_t &bedBegInd, const size_t &binBegInd, const size_t &randVecLen){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	const uint8_t remainderInd = byteSize_ - static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	// Create a vector to store random bytes for stochastic heterozygote resolution
	vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	// Fill the random byte vector
	for (auto &rv : rand){
		rv = rng_.ranInt();
	}
	size_t iIndiv         = 0;
	const size_t begIndiv = bedBegInd * nIndividuals_;
	vector<uint8_t> missMasks(locusSize_, 0);
	size_t i0Byte = 0;                                                                         // to index the missMasks vector
	for (size_t iByte = binBegInd; iByte < binBegInd + locusSize_ - 1; ++iByte){                                   // treat the last byte separately
		for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte){
			uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);     // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= 0b10000011;                                            // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
			missMasks[i0Byte]        |= (missingMask << iInByte);
			curIndiv                 &= 0b00000011;
			const uint8_t randMask    = (randBytes[iByte] >> iInByte) & oneBit_;               // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
			binGenotypes_[iByte]     |= curBitMask << iInByte;
			++iIndiv;
		}
		++i0Byte;
	}
	// now deal with the last byte in the individual
	for (uint8_t iRem = 0; iRem < remainderInd; ++iRem){
		uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);         // cramming down to one byte because I do not care what the actual value is
		curIndiv                 &= 0b10000011;                                                // mask everything in the middle
		const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
		missMasks.back()         |= (missingMask << iRem);
		curIndiv                 &= 0b00000011;                                                
		const uint8_t randMask    = (randBytes[locusSize_ - 1] >> iRem) & oneBit_;             // 0b00000000 or 0b00000001 with equal chance
		uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                   // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
		curBitMask               &= ~missingMask;                                              // zero it out if missing value is set

		binGenotypes_[binBegInd + locusSize_ - 1] |= curBitMask << iRem;
		++iIndiv;
	}
	float maf = static_cast<float>( countSetBits(binGenotypes_, binBegInd, locusSize_) ) / static_cast<float>(nIndividuals_);
	if (maf > 0.5){ // always want the alternative to be the minor allele
		i0Byte = 0;
		for (size_t i = binBegInd; i < binBegInd + locusSize_; ++i){
			binGenotypes_[i] = (~binGenotypes_[i]) & (~missMasks[i0Byte]);
			++i0Byte;
		}
		const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
		binGenotypes_[binBegInd + locusSize_ - 1] &= lastByteMask; // unset the remainder bits
		maf = 1.0 - maf;
	}
	aaf_.push_back(maf);
}

void GenoTableBin::jaccardBlock_(const size_t &iLocus, const size_t &blockInd, vector<float> &jaccardVec) const {
	vector<uint8_t> locus(locusSize_);
	size_t ind = blockInd;
	for (size_t jCol = iLocus + 1; jCol < nLoci_; ++jCol){
		size_t rowInd = iLocus * locusSize_;
		size_t colInd = jCol * locusSize_;
		for (size_t iBinLoc = 0; iBinLoc < locusSize_; ++iBinLoc){
			locus[iBinLoc] = binGenotypes_[rowInd + iBinLoc] & binGenotypes_[colInd + iBinLoc];
		}
		const uint32_t uni = countSetBits(locus);
		for (size_t iBinLoc = 0; iBinLoc < locusSize_; ++iBinLoc){
			locus[iBinLoc] = binGenotypes_[rowInd + iBinLoc] | binGenotypes_[colInd + iBinLoc];
		}
		const uint32_t isect = countSetBits(locus);
		jaccardVec[ind]      = static_cast<float>(uni) / static_cast<float>(isect) - aaf_[iLocus] * aaf_[jCol];
		++ind;
	}
}

// GenoTableHash methods
constexpr array<char, 3> GenoTableHash::magicBytes_ = {0x6c, 0x1b, 0x01};              // Leading bytes for .bed files 
constexpr uint8_t  GenoTableHash::oneBit_           = 0b00000001;                      // One set bit for masking 
constexpr uint8_t  GenoTableHash::byteSize_         = 8;                               // Size of one byte in bits 
constexpr uint8_t  GenoTableHash::llWordSize_       = 8;                               // 64 bit word size in bytes 
constexpr size_t   GenoTableHash::nblocks_          = sizeof(size_t) / 4;              // MurMurHash number of blocks 
constexpr uint32_t GenoTableHash::mmhKeyLen_        = sizeof(size_t);                  // MurMurHash key length 
constexpr uint16_t GenoTableHash::emptyBinToken_    = numeric_limits<uint16_t>::max(); // Value corresponding to an empty token 
constexpr uint32_t GenoTableHash::c1_               = 0xcc9e2d51;                      // MurMurHash c1 constant 
constexpr uint32_t GenoTableHash::c2_               = 0x1b873593;                      // MurMurHash c2 constant 

// Constructors
GenoTableHash::GenoTableHash(const string &inputFileName, const size_t &nIndividuals, const size_t &kSketches) : nIndividuals_{nIndividuals}, kSketches_{kSketches}, nLoci_{0} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in ") + string(__FUNCTION__);
	} else if (nIndividuals > numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
		throw string("ERROR: the number of individuals (") + to_string(nIndividuals) + string( ") is too big to make a square relationship matrix in ") + string(__FUNCTION__);
	}
	sketchSize_ = nIndividuals_ / kSketches + static_cast<bool>(nIndividuals_ % kSketches);
	if (sketchSize_ >= emptyBinToken_){
		throw string("ERROR: Number of sketches (") + to_string(kSketches_) + string(") implies sketch size (") +
			to_string(sketchSize_) + string(") that is larger than ") + to_string(emptyBinToken_) +
			string( ", the largest allowed value in ") + string(__FUNCTION__);
	}
	const size_t nBedBytes = nIndividuals_ / 4 + static_cast<bool>(nIndividuals_ % 4);
	fstream inStr;
	// Start by measuring file size
	inStr.open(inputFileName.c_str(), ios::in | ios::binary | ios::ate);
	const size_t N = static_cast<uint64_t>( inStr.tellg() ) - 3UL;       // first three bytes are majic
	inStr.close();
	nLoci_ = N / nBedBytes;
	// Calculate the actual sketch number based on the realized sketch size
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	inStr.open(inputFileName.c_str(), ios::in | ios::binary);
	char magicBuf[magicBytes_.size()]{};
	inStr.read( magicBuf, magicBytes_.size() );
	if ( inStr.eof() ){
		throw string("ERROR: No loci in the input .bed file ") + inputFileName + string(" in ") + string(__FUNCTION__);
	} else if (magicBuf[0] != magicBytes_[0]){
		throw string("ERROR: first magic byte in input .bed file is not the expected value in ") + string(__FUNCTION__);
	} else if (magicBuf[1] != magicBytes_[1]){
		throw string("ERROR: second magic byte in input .bed file is not the expected value in ") + string(__FUNCTION__);
	} else if (magicBuf[2] != magicBytes_[2]){
		throw string("ERROR: third magic byte in input .bed file is not the expected value in ") + string(__FUNCTION__);
	}
	// Generate the binary genotype table while reading the .bed file
	locusSize_              = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	const size_t ranVecSize = nBedBytes / llWordSize_ + static_cast<bool>(nBedBytes % llWordSize_);
	vector<char> bedLocus(nBedBytes, 0);

	// generate the sequence of random integers; each column must be permuted the same
	vector<size_t> ranInts;
	size_t i = nIndividuals_;
	while (i >= 2UL){
		ranInts.push_back( rng_.ranInt() % i ); // need 0 <= j <= i, so i is actually i+1 (compared to the Wikipedia description)
		--i;
	}
	vector<uint32_t> seeds;
	seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
	size_t locusBegin = 0;
	while ( inStr.read(bedLocus.data(), nBedBytes) ) {
		bed2oph_(bedLocus, locusBegin, nBedBytes, ranVecSize, ranInts, seeds);
	}
	inStr.close();
}

GenoTableHash::GenoTableHash(const vector<int> &maCounts, const size_t &nIndividuals, const size_t &kSketches) : nIndividuals_{nIndividuals}, kSketches_{kSketches}, nLoci_{maCounts.size() / nIndividuals} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in ") + string(__FUNCTION__);
	}
	if (maCounts.size() % nIndividuals){
		throw string("ERROR: length of allele count vector (") + to_string( maCounts.size() ) + string(" is not divisible by the provided number of individuals (") +
			to_string(nIndividuals) + string( ") in ") + string(__FUNCTION__);
	}
	if ( maCounts.empty() ){
		throw string("ERROR: empty vector of minor allele counts in ") + string(__FUNCTION__);
	}
	sketchSize_ = nIndividuals_ / kSketches_ + static_cast<bool>(nIndividuals_ % kSketches_);
	if (sketchSize_ >= emptyBinToken_){
		throw string("ERROR: Number of sketches (") + to_string(kSketches_) + string(") implies sketch size (") +
			to_string(sketchSize_) + string(") that is larger than ") + to_string(emptyBinToken_) +
			string( ", the largest allowed value in ") + string(__FUNCTION__);
	}
	locusSize_              = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	const size_t ranVecSize = locusSize_ / llWordSize_ + static_cast<bool>(locusSize_ % llWordSize_);
	// Calculate the actual sketch number based on the realized sketch size
	sketches_.resize(kSketches_ * nLoci_, emptyBinToken_);
	// generate the sequence of random integers; each column must be permuted the same
	vector<size_t> ranInts;
	size_t i = nIndividuals_;
	while (i >= 2UL){
		ranInts.push_back( rng_.ranInt() % i ); // need 0 <= j <= i, so i is actually i+1 (compared to the Wikipedia description)
		--i;
	}
	vector<uint32_t> seeds;
	seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
	for (size_t jLoc = 0; jLoc < nLoci_; ++jLoc) {
		mac2oph_(maCounts, jLoc, ranVecSize, ranInts, seeds);
	}
}

GenoTableHash::GenoTableHash(GenoTableHash &&in) noexcept : kSketches_{in.kSketches_} {
	if (this != &in){
		sketches_     = move(in.sketches_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;
		sketchSize_   = in.sketchSize_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;
		in.sketchSize_   = 0;
	}
}

GenoTableHash& GenoTableHash::operator=(GenoTableHash &&in) noexcept {
	if (this != &in){
		*this = move(in);
	}
	return *this;
}

vector<float> GenoTableHash::allHashLD() const {
	if ( sketches_.empty() ){
		throw string("ERROR: Cannot calculate hash-based LD on empty sketches in ") + string(__FUNCTION__);
	}
	vector<float> LDmat(nLoci_ * (nLoci_ - 1) / 2, 0.0);
	const size_t kSketches = sketches_.size() / nLoci_;
	const float fNind = 1.0 / static_cast<float>(kSketches);
	vector< future<void> > tasks;
	tasks.reserve(nLoci_);
	size_t blockBeg = 0; // index in LDmat of the first element of the block of similarities
	for (size_t iRow = 0; iRow < nLoci_; ++iRow){
		tasks.emplace_back( async([this, iRow, blockBeg, kSketches, fNind, &LDmat]{hashJacBlock_(iRow, blockBeg, kSketches, fNind, LDmat); }) );
		blockBeg += nLoci_ - iRow - 1;
	}
	return LDmat;
}

vector<uint16_t> GenoTableHash::assignGroups(const size_t &nElements) const {
	const size_t kSketches = sketches_.size() / nLoci_;
	if (kSketches == 0){
		throw string("ERROR: Number of sketches must be non-zero in ") + string(__FUNCTION__);
	}
	if (nElements > kSketches){
		throw string("ERROR: Number of elements to consider (") + to_string(nElements) + string(") is greater than the number of sketches (") 
						+ to_string(kSketches) + string( ") in ") + string(__FUNCTION__);
	}
	vector<uint16_t> grpID;
	grpID.reserve(nLoci_);
	const uint32_t seed = static_cast<uint32_t>( rng_.ranInt() );
	for (size_t locusBeg = 0; locusBeg < sketches_.size(); locusBeg += kSketches){
		grpID.push_back( murMurHash_(locusBeg, nElements, seed) );
	}
	return grpID;
}

vector<uint16_t> GenoTableHash::assignGroups() const {
	const size_t kSketches = sketches_.size() / nLoci_;
	if (kSketches == 0){
		throw string("ERROR: Number of sketches must be non-zero in ") + string(__FUNCTION__);
	}
	vector<uint16_t> grpID;
	grpID.reserve(nLoci_);
	const uint32_t seed = static_cast<uint32_t>( rng_.ranInt() );
	for (size_t locusBeg = 0; locusBeg < sketches_.size(); locusBeg += kSketches){
		grpID.push_back( simHash_(locusBeg, kSketches, seed) );
	}
	return grpID;
}

void GenoTableHash::groupByLD(const uint16_t &hammingCutoff, const size_t &kSketches, const size_t &lookBackNumber, const string &outFileName) const {
	if ( sketches_.empty() ){
		throw string("ERROR: no OPH sketches generated before calling in ") + string(__FUNCTION__);
	}
	const size_t totSketches = sketches_.size() / nLoci_;
	if (kSketches == 0){
		throw string("ERROR: number of OPH sketches to simHash cannot be 0 in ") + string(__FUNCTION__);
	} else if (kSketches > totSketches){
		throw string("ERROR: number of OPH sketches to simHash cannot exceed the total number of sketches in ") + string(__FUNCTION__);
	}
	const uint32_t seed = static_cast<uint32_t>( rng_.ranInt() );
	vector< vector<size_t> > ldGroup;                              // each element is a vector of indexes into sketches_
	vector<uint16_t> activeHashes(lookBackNumber, 0);              // hashes that are under consideration in a simplified ring buffer

	auto time1 = high_resolution_clock::now();
	// Start by filling the active hash buffer
	//
	// initialize with the first locus
	activeHashes.back() = simHash_(0, kSketches, seed);
	auto latestHashIt   = activeHashes.rbegin();                   // points to the latest hash under consideration
	size_t locusInd     = totSketches;                             // current locus index (starts from the second locus)
	ldGroup.emplace_back(vector<size_t>{0});
	// loop through loci until buffer is full or no more loci left
	while ( ( latestHashIt != activeHashes.rend() ) && ( locusInd < sketches_.size() ) ){
		const uint16_t curHash = simHash_(locusInd, kSketches, seed);
		auto fwdHashIt         = activeHashes.cbegin() + distance( latestHashIt, activeHashes.rend() ) - 1;
		auto ldgIt             = ldGroup.rbegin();   // latest group is at the end
		bool noCollision       = true;
		while ( fwdHashIt != activeHashes.cend() ){
			if ( ldgIt == ldGroup.rend() ){
				throw string("ERROR: got past the first group in ") + string(__FUNCTION__) + string(" while filling in the buffer") +
					string("\n This should never happen. Please report this as a bug to github.com/tonymugen");
			}
			if (hammingDistance_(*fwdHashIt, curHash) <= hammingCutoff){
				noCollision = false;
				break;
			}
			++fwdHashIt;
			++ldgIt;
		}
		if (noCollision){
			++latestHashIt;
			if ( latestHashIt == activeHashes.rend() ){
				ldGroup.emplace_back(vector<size_t>{locusInd});
				locusInd       += totSketches;
				activeHashes[0] = curHash;
				break;
			} else {
				*latestHashIt = curHash;
				ldGroup.emplace_back(vector<size_t>{locusInd});
				locusInd += totSketches;
			}
		} else {
			ldgIt->push_back(locusInd);
			locusInd += totSketches;
		}
	}
	// If there are loci left after filling the buffer, continue.
	while ( locusInd < sketches_.size() ){
		const uint16_t curHash = simHash_(locusInd, kSketches, seed);
		size_t latestHashInd   = 0;
		auto ldgIt             = ldGroup.rbegin();   // latest group is at the end
		bool noCollision       = true;
		for (size_t iTestHash = 0; iTestHash < lookBackNumber; ++iTestHash){
			if ( ldgIt == ldGroup.rend() ){
				throw string("ERROR: got past the first group in ") + string(__FUNCTION__) +
					string("\n This should never happen. Please report this as a bug to github.com/tonymugen");
			}
			if (hammingDistance_(activeHashes[(latestHashInd + iTestHash) % lookBackNumber], curHash) <= hammingCutoff){
				noCollision = false;
				break;
			}
			++ldgIt;
		}
		if (noCollision){
			++latestHashInd;
			latestHashInd %= lookBackNumber;
			ldGroup.emplace_back(vector<size_t>{locusInd});
			locusInd += totSketches;
		} else {
			ldgIt->push_back(locusInd);
			locusInd += totSketches;
		}
	}
	// estimate Jaccard similarities within groups
	//vector< future<void> > tasks;
	//tasks.reserve( ldGroup.size() );
	vector< vector<float> > hashJacGroups;
	const float fNind = 1.0 / static_cast<float>(totSketches);
	for (const auto &ldg : ldGroup){
		if (ldg.size() >= 2){
			hashJacGroups.emplace_back( hashJacBlock_(ldg, totSketches, fNind) );
		}
	}
	auto time2 = high_resolution_clock::now();
	duration<float, milli> execTime = time2 - time1;
	// Save the results
	fstream out;
	out.open(outFileName.c_str(), ios::out | ios::trunc);
	out << "# " << execTime.count() << "\n";
	size_t iValidLDG = 0;
	for (const auto &ldg : ldGroup){
		if (ldg.size() >= 2){
			size_t ldValInd = 0;
			for (size_t iRow = 0; iRow < ldg.size() - 1; ++iRow) {
				for (size_t jCol = iRow + 1; jCol < ldg.size(); ++jCol){
					out << ldg[iRow] / totSketches + 1 << "\t" << ldg[jCol] / totSketches + 1 << "\t" << hashJacGroups[iValidLDG][ldValInd] << "\n";
					++ldValInd;
				}
			}
			++iValidLDG;
		}
	}
	out.close();
}

void GenoTableHash::locusOPH_(const vector<size_t> &permutation, vector<uint32_t> &seeds, vector<uint8_t> &binLocus){
	// Start with a permutation to make OPH
	size_t iIndiv = nIndividuals_ - 1UL; // safe b/c nIndividuals_ > 1 is checked at construction
	for (const auto &ri : permutation){
		uint16_t firstIdx  = iIndiv % byteSize_;
		size_t firstByte   = (iIndiv / byteSize_);
		uint16_t secondIdx = ri % byteSize_;
		size_t secondByte  = (ri / byteSize_);
		uint16_t diff      = byteSize_ * (firstByte != secondByte); // will be 0 if the same byte is being accessed; then need to swap bits within byte

		// swapping bits within a two-byte variable
		// using the method in https://graphics.stanford.edu/~seander/bithacks.html#SwappingBitsXOR
		// if the same byte is being accessed, secondIdx is not shifted to the second byte
		// This may be affected by endianness (did not test)
		uint16_t twoBytes  = (static_cast<uint16_t>(binLocus[secondByte]) << 8) | ( static_cast<uint16_t>(binLocus[firstByte]) );
		secondIdx         += diff;
		uint16_t x         = ( (twoBytes >> firstIdx) ^ (twoBytes >> secondIdx) ) & 1;
		twoBytes          ^= ( (x << firstIdx) | (x << secondIdx) );

		memcpy( binLocus.data() + firstByte, &twoBytes, sizeof(uint8_t) );
		twoBytes = twoBytes >> diff;
		memcpy( binLocus.data() + secondByte, &twoBytes, sizeof(uint8_t) );
		--iIndiv;
	}
	// Now make the sketches
	vector<size_t> filledIndexes;                     // indexes of the non-empty sketches
	size_t iByte     = 0;
	size_t colEnd    = iByte + locusSize_;
	size_t sketchBeg = nLoci_ * kSketches_;            // nLoci_ is incremented at the end of this loop
	size_t iSeed     = 0;                             // index into the seed vector
	uint8_t iInByte  = 0;
	// A possible optimization is to test a whole byte for 0
	// Will test later
	for (size_t iSketch = 0; iSketch < kSketches_; ++iSketch){
		uint16_t firstSetBitPos = 0;
		while ( ( ( (oneBit_ << iInByte) & binLocus[iByte] ) == 0 ) &&
				(firstSetBitPos < sketchSize_) && (iByte != colEnd) ){
			++iInByte;
			// these are instead of an if statement
			iByte  += iInByte == byteSize_;
			iInByte = iInByte % byteSize_;
			++firstSetBitPos;
		}
		if ( (iByte < colEnd) && (firstSetBitPos < sketchSize_) ){
			filledIndexes.push_back(iSketch);
			sketches_[sketchBeg + iSketch] = firstSetBitPos;

			uint16_t remainder = sketchSize_ - firstSetBitPos;
			uint16_t inByteMod = remainder % byteSize_;
			uint16_t inByteSum = iInByte + inByteMod;

			iByte  += remainder / byteSize_ + inByteSum / byteSize_;
			iInByte = inByteSum % byteSize_;
		}
	}
	if (filledIndexes.size() == 1){
		for (size_t iSk = 0; iSk < kSketches_; ++iSk){ // this will overwrite the one assigned sketch, but the wasted operation should be swamped by the rest
			sketches_[sketchBeg + iSk] = sketches_[filledIndexes[0] + sketchBeg];
		}
	} else if (filledIndexes.size() != kSketches_){
		size_t emptyCount = kSketches_ - filledIndexes.size();
		while (emptyCount){
			for (const auto &f : filledIndexes){
				uint32_t newIdx = murMurHash_(f, seeds[iSeed]) % kSketches_ + sketchBeg;
				if ( sketches_[newIdx] == emptyBinToken_ ){
					sketches_[newIdx] = sketches_[f + sketchBeg];
					--emptyCount;
					break;
				}
			}
			++iSeed;
			mutex mtx;
			lock_guard<mutex> lk(mtx);      // lock before measuring to ensure that the size is valid
			if ( iSeed == seeds.size() ){
				seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
			}
		}
	}
}

void GenoTableHash::bed2oph_(const vector<char> &bedData, const size_t &begInd, const size_t &locusLength, const size_t &randVecLen, const vector<size_t> &permutation, vector<uint32_t> &seeds){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	const size_t endWholeBed   = begInd + locusLength - 2UL + (locusLength & 1UL);
	const size_t addIndv       = nIndividuals_ - endWholeBed * 4UL;
	// Fill the random byte vector
	vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	for (auto &rv : rand){
		rv = rng_.ranInt();
	}
	vector<uint8_t> binLocus(locusSize_, 0);
	vector<uint8_t> missMasks(locusSize_, 0);
	size_t iBinGeno = 0;                       // binLocus vector index
	uint8_t bedByte = 0;
	// Two bytes of .bed code go into one byte of my binary representation
	// Therefore, work on two consecutive bytes of .bed code in the loop
	for (size_t iBed = begInd; iBed < endWholeBed ; iBed += 2){                // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
		bedByte             = ~bedData[iBed];                                  // flip so that homozygous second allele (usually minor) is set to 11
		uint8_t offsetToBin = 0;                                               // move the .bed mask by this much to align with the binarized byte
		for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
			uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
			uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
			// Keep track of missing genotypes to revert them if I have to flip bits later on
			const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
			missMasks[iBinGeno]      |= curMissMask >> offsetToBin;
			// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
			secondBitMask      |= randBytes[iBed] & (firstBitMask << 1);
			firstBitMask       &= secondBitMask >> 1;
			binLocus[iBinGeno] |= firstBitMask >> offsetToBin;
			++offsetToBin;
		}
		const size_t nextIbed = iBed + 1;
		bedByte               = ~bedData[nextIbed];
		for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
			uint8_t firstBitMask  = bedByte & (oneBit_ << iInByteG);
			uint8_t secondBitMask = bedByte & ( oneBit_ << (iInByteG + 1) );
			// Keep track of missing genotypes to revert them if I have to flip bits later on
			const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
			missMasks[iBinGeno]      |= curMissMask << offsetToBin;
			// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
			secondBitMask      |= randBytes[nextIbed] & (firstBitMask << 1);
			firstBitMask       &= secondBitMask >> 1;
			binLocus[iBinGeno] |= firstBitMask << offsetToBin; // keep adding to the current binarized byte, so switch the direction of shift
			--offsetToBin;
		}
		++iBinGeno;
	}
	uint8_t inBedByteOffset = 0;
	for (size_t iInd = 0; iInd < addIndv; ++iInd){
		const size_t curBedByte = endWholeBed + iInd / 4;
		uint8_t firstBitMask    = (~bedData[curBedByte]) & (oneBit_ << inBedByteOffset);
		const uint8_t secondBBO = inBedByteOffset + 1;
		uint8_t secondBitMask   = (~bedData[curBedByte]) & (oneBit_ << secondBBO);
		// Keep track of missing genotypes to revert them if I have to flip bits later on
		const uint8_t curMissMask = ( ( secondBitMask ^ (firstBitMask << 1) ) & secondBitMask ) >> 1;  // 2nd different from 1st, and 2nd set => missing
		missMasks.back()         |= (curMissMask >> inBedByteOffset) << iInd;
		// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
		secondBitMask   |= randBytes[curBedByte] & (firstBitMask << 1);
		firstBitMask    &= secondBitMask >> 1;
		firstBitMask     = firstBitMask >> inBedByteOffset;
		firstBitMask     = firstBitMask << iInd;
		binLocus.back() |= firstBitMask;
		inBedByteOffset += 2;
		inBedByteOffset  = inBedByteOffset % 8;
	}
	float aaCount = static_cast<float>( countSetBits(binLocus) ) / static_cast<float>(nIndividuals_);
	if (aaCount > 0.5){ // always want the alternative to be the minor allele
		for (size_t iBL = 0; iBL < locusSize_; ++iBL){
			binLocus[iBL] = (~binLocus[iBL]) & (~missMasks[iBL]);
		}
		const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
		binLocus.back() &= lastByteMask; // unset the remainder bits
		aaCount = 1.0 - aaCount;
	}
	aaf_.push_back(aaCount);
	locusOPH_(permutation, seeds, binLocus);
}

void GenoTableHash::mac2oph_(const vector<int> &macData, const size_t &begInd, const size_t &randVecLen, const vector<size_t> &permutation, vector<uint32_t> &seeds){
	// Define constants. Some can be taken outside of the function as an optimization
	// Opting for more encapsulation for now unless I find significant performance penalties
	const uint8_t remainderInd = byteSize_ - static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	// Create a vector to store random bytes for stochastic heterozygote resolution
	vector<uint64_t> rand(randVecLen);
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	// Fill the random byte vector
	for (auto &rv : rand){
		rv = rng_.ranInt();
	}
	size_t iIndiv         = 0;
	const size_t begIndiv = begInd * nIndividuals_;
	vector<uint8_t> missMasks(locusSize_, 0);
	vector<uint8_t> binLocus(locusSize_, 0);
	size_t i0Byte = 0;                                                                         // to index the missMasks vector
	for (size_t iByte = 0; iByte < locusSize_ - 1; ++iByte){                                   // treat the last byte separately
		for (uint8_t iInByte = 0; iInByte < byteSize_; ++iInByte){
			uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);     // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= 0b10000011;                                            // mask everything in the middle
			const uint8_t missingMask = curIndiv >> 7;                                         // 0b00000001 iff is missing (negative value)
			missMasks[i0Byte]        |= (missingMask << iInByte);
			curIndiv                 &= 0b00000011;
			const uint8_t randMask    = (randBytes[iByte] >> iInByte) & oneBit_;               // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= ~missingMask;                                          // zero it out if missing value is set
			binLocus[iByte]          |= curBitMask << iInByte;
			++iIndiv;
		}
		++i0Byte;
	}
	// now deal with the last byte in the individual
	for (uint8_t iRem = 0; iRem < remainderInd; ++iRem){
		uint8_t curIndiv          = static_cast<uint8_t>(macData[begIndiv + iIndiv]);         // cramming down to one byte because I do not care what the actual value is
		curIndiv                 &= 0b10000011;                                                // mask everything in the middle
		const uint8_t missingMask = curIndiv >> 7;                                             // 0b00000001 iff is missing (negative value)
		missMasks.back()         |= (missingMask << iRem);
		curIndiv                 &= 0b00000011;                                                
		const uint8_t randMask    = (randBytes[locusSize_ - 1] >> iRem) & oneBit_;             // 0b00000000 or 0b00000001 with equal chance
		uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);                   // if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
		curBitMask               &= ~missingMask;                                              // zero it out if missing value is set

		binLocus[locusSize_ - 1] |= curBitMask << iRem;
		++iIndiv;
	}
	float maf = static_cast<float>( countSetBits(binLocus) ) / static_cast<float>(nIndividuals_);
	if (maf > 0.5){ // always want the alternative to be the minor allele
		i0Byte = 0;
		for (size_t i = 0; i < locusSize_; ++i){
			binLocus[i] = (~binLocus[i]) & (~missMasks[i0Byte]);
			++i0Byte;
		}
		const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;
		binLocus[locusSize_ - 1] &= lastByteMask; // unset the remainder bits
		maf = 1.0 - maf;
	}
	aaf_.push_back(maf);
	throw string("got here");
	locusOPH_(permutation, seeds, binLocus);
}

uint32_t GenoTableHash::murMurHash_(const size_t &key, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	auto blocks = reinterpret_cast<const uint32_t *>(&key);

	for(size_t i = 0; i < nblocks_; ++i){
		uint32_t k1 = blocks[i];

		k1 *= c1_;
		k1  = (k1 << 15) | (k1 >> 17);
		k1 *= c2_;

		hash ^= k1;
		hash  = (hash << 13) | (hash >> 19);
		hash  = hash * 5 + 0xe6546b64;
	}

	// there is no tail since the input is fixed (at eight bytes typically)
	// finalize
	hash ^= mmhKeyLen_;
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;

	return hash;
}

uint16_t GenoTableHash::murMurHash_(const uint16_t &key, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	uint32_t block = static_cast<uint32_t>(key);

	uint32_t k1 = block;

	k1 *= c1_;
	k1 = (k1 << 15) | (k1 >> 17);
	k1 *= c2_;

	hash ^= k1;
	hash  = (hash << 13) | (hash >> 19);
	hash  = hash * 5 + 0xe6546b64;

	// there is no tail since the input is fixed (at eight bytes typically)
	// finalize
	hash ^= mmhKeyLen_;
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;
	
	auto hashV = reinterpret_cast<const uint16_t *>(&hash);

	return hashV[0] ^ hashV[1];
}

uint32_t GenoTableHash::murMurHash_(const size_t &startInd, const size_t &nElements, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	auto blocks    = reinterpret_cast<const uint32_t *>(sketches_.data() + startInd);
	size_t nBlocks = nElements / 2; // each sketch is 16 bits; this implies that only even number of sketches is considered

	for(size_t i = 0; i < nBlocks; ++i){
		uint32_t k1 = blocks[i];

		k1 *= c1_;
		k1 = (k1 << 15) | (k1 >> 17);
		k1 *= c2_;

		hash ^= k1;
		hash  = (hash << 13) | (hash >> 19);
		hash  = hash * 5 + 0xe6546b64;
	}

	// there is no tail since the input is fixed (at eight bytes typically)
	// finalize
	hash ^= mmhKeyLen_;
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;

	return hash;
}

uint16_t GenoTableHash::simHash_(const size_t &startInd, const size_t &kSketches, const uint32_t &seed) const {
	uint16_t hash        = 0;
	const uint16_t one   = 1;
	const size_t twoByte = 2 * byteSize_;
	array<int16_t, twoByte> v{};
	for (size_t iSketch = startInd; iSketch < startInd + kSketches; ++iSketch){
		const uint16_t skHash = murMurHash_(sketches_[iSketch], seed);
		for (uint16_t j = 0; j < twoByte; ++j){
			v[j] += -1 + 2 * static_cast<int16_t>( one & (skHash >> j) );
		}
	}
	for (uint16_t i = 0; i < twoByte; ++i){
		hash |= (static_cast<uint16_t>(v[i] > 0) << i);
	}
	return hash;
}

void GenoTableHash::hashJacBlock_(const size_t &iLocus, const size_t &blockInd, const size_t &kSketches, const float &invK, vector<float> &hashJacVec) const {
	size_t ind = blockInd;
	for (size_t jCol = iLocus + 1; jCol < nLoci_; ++jCol){
		float simVal = 0.0;
		size_t sketchRowInd = iLocus * kSketches;
		size_t sketchColInd = jCol * kSketches;
		for (size_t iSk = 0; iSk < kSketches; ++iSk){
			if (sketches_[sketchRowInd + iSk] == sketches_[sketchColInd + iSk]){
				simVal += 1.0;
			}
		}
		simVal *= invK;
		simVal -= aaf_[iLocus] * aaf_[jCol]; // subtracting expected similarity
		hashJacVec[ind] = simVal;
		++ind;
	}
}

void GenoTableHash::hashJacBlock_(const size_t &iLocus, const vector<size_t> &jLocus, const size_t &kSketches, const float &invK, vector<float> &hashJacVec) const {
	for (size_t jCol = iLocus + 1; jCol < jLocus.size(); ++jCol){
		float simVal = 0.0;
		for (size_t iSk = 0; iSk < kSketches; ++iSk){
			if (sketches_[jLocus[iLocus] + iSk] == sketches_[jLocus[jCol] + iSk]){
				simVal += 1.0;
			}
		}
		simVal *= invK;
		simVal -= aaf_[jLocus[iLocus] / kSketches] * aaf_[jLocus[jCol] / kSketches]; // subtracting expected similarity
		hashJacVec.push_back(simVal);
	}
}

vector<float> GenoTableHash::hashJacBlock_(const vector<size_t> &locusIndexes, const size_t &kSketches, const float &invK) const {
	//vector< future<void> > tasks;
	vector<float> hashJacVec;
	hashJacVec.reserve(locusIndexes.size() * (locusIndexes.size() - 1) / 2);
	//tasks.reserve(locusIndexes.size() - 1); // size checked outside of the function, safe to subtract
	for (size_t iLocus = 0; iLocus < locusIndexes.size() - 1; ++iLocus) {
		hashJacBlock_(iLocus, locusIndexes, kSketches, invK, hashJacVec);
		// it would be cool to multi-thread it, but must take care not to scramble the indexes
		//tasks.emplace_back( async([this, iLocus, &locusIndexes, kSketches, invK, &hashJacVec]{hashJacBlock_(iLocus, locusIndexes, kSketches, invK, hashJacVec);}) );
	}
	return hashJacVec;
}
