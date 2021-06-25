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

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <fstream>

#include <iostream>
#include <bitset>

#include "gvarHash.hpp"

using std::vector;
using std::array;
using std::string;
using std::to_string;
using std::move;
using std::numeric_limits;
using std::fstream;
using std::ios;
using std::streampos;

using namespace BayesicSpace;

const array<char, 3> GenoTable::magicBytes_ = {0x6c, 0x1b, 0x01};
const uint8_t GenoTable::oneBit_            = 0b00000001;
const uint8_t GenoTable::byteSize_          = 8;
const uint8_t GenoTable::llWordSize_        = 8;
const size_t GenoTable::nblocks_            = sizeof(size_t) / 4;
const uint32_t GenoTable::mmhKeyLen_        = sizeof(size_t);
const uint16_t GenoTable::emptyBinToken_    = numeric_limits<uint16_t>::max();
const uint32_t GenoTable::c1_               = 0xcc9e2d51;
const uint32_t GenoTable::c2_               = 0x1b873593;

// Constructors
GenoTable::GenoTable(const string &inputFileName, const size_t &nIndividuals) : nIndividuals_{nIndividuals}, nLoci_{0} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in the GenoTable(const string &, const size_t &) constructor");
	} else if (nIndividuals > numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
		throw string("ERROR: the number of individuals (") + to_string(nIndividuals) + string(") is too big to make a square relationship matrix in the GenoTable(const string &, const size_t &) constructor");
	}
	size_t nBedBytes = nIndividuals_ / 4 + static_cast<bool>(nIndividuals_ % 4);
	fstream inStr;
	inStr.open(inputFileName.c_str(), ios::in | ios::binary);
	char magicBuf[magicBytes_.size()]{};
	inStr.read( magicBuf, magicBytes_.size() );
	if ( inStr.eof() ){
		throw string("ERROR: No loci in the input .bed file ") + inputFileName + string("in the GenoTable(const string &, const size_t &) constructor");
	} else if (magicBuf[0] != magicBytes_[0]){
		throw string("ERROR: first magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	} else if (magicBuf[1] != magicBytes_[1]){
		throw string("ERROR: second magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	} else if (magicBuf[2] != magicBytes_[2]){
		throw string("ERROR: third magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	}
	// Generate the binary genotype table while reading the .bed file
	locusSize_                 = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
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
		size_t iBinGeno = 0;                       // binGenotypes_ vector index
		// Two bytes of .bed code go into one byte of my binary representation
		// Therefore, work on two consecutive bytes of .bed code in the loop
		for (size_t iBed = 0; iBed < endBed ; iBed += 2){                     // the last byte has the padding; will deal with it separately (plus the penultimate byte if nBedBytes is even)
			bedLocus[iBed] = ~bedLocus[iBed];                                 // flip so that homozygous second allele (usually minor) is set to 11
			uint8_t offsetToBin = 0;                                          // move the .bed mask by this much to align with the binarized byte
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedLocus[iBed] & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedLocus[iBed] & ( oneBit_ << (iInByteG + 1) );
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask |= randBytes[iBed] & (firstBitMask << 1);
				firstBitMask  &= secondBitMask >> 1;
				binLocus[iBinGeno] |= firstBitMask >> offsetToBin;
				offsetToBin++;
			}
			const size_t nextIbed = iBed + 1;
			bedLocus[nextIbed]    = ~bedLocus[nextIbed];
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				uint8_t firstBitMask  = bedLocus[nextIbed] & (oneBit_ << iInByteG);
				uint8_t secondBitMask = bedLocus[nextIbed] & ( oneBit_ << (iInByteG + 1) );
				// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
				secondBitMask |= randBytes[nextIbed] & (firstBitMask << 1);
				firstBitMask  &= secondBitMask >> 1;
				binLocus[iBinGeno] |= firstBitMask << offsetToBin; // keep adding to the current binarized byte, so switch the direction of shift
				offsetToBin--;
			}
			iBinGeno++;
		}
		for (size_t iBed = endBed; iBed < nBedBytes; iBed++){
			bedLocus[iBed] = ~bedLocus[iBed];
		}
		uint8_t inBedByteOffset = 0;
		for (size_t iInd = 0; iInd < addIndv; iInd++){
			const size_t curBedByte = endBed + iInd / 4;
			uint8_t firstBitMask    = bedLocus[curBedByte] & (oneBit_ << inBedByteOffset);
			inBedByteOffset++;
			uint8_t secondBitMask = bedLocus[curBedByte] & (oneBit_ << inBedByteOffset);
			// If 1st is set and 2nd is not, we have a heterozygote. In this case, set the 1st with a 50/50 chance
			secondBitMask   |= randBytes[curBedByte] & (firstBitMask << 1);
			firstBitMask    &= secondBitMask >> 1;
			firstBitMask     = firstBitMask >> (inBedByteOffset - 1);
			firstBitMask     = firstBitMask << iInd;
			binLocus.back() |= firstBitMask;
			inBedByteOffset++;
			inBedByteOffset = inBedByteOffset % 8;
		}
		float aaCount = static_cast<float>( countSetBits_(binLocus) ) / fNind;
		if (aaCount > 0.5){ // always want the alternative to be the minor allele
			for (auto &bg : binLocus){
				binGenotypes_.push_back(~bg);
			}
			binGenotypes_.back() &= lastByteMask; // unset the remainder bits
			aaCount = 1.0 - aaCount;
		} else {
			for (const auto &bg : binLocus){
				binGenotypes_.push_back(bg);
			}
		}
		aaf_.push_back(aaCount);
		nLoci_++;
	}
	inStr.close();
}

GenoTable::GenoTable(const vector<int> &maCounts, const size_t &nIndividuals) : nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	if (maCounts.size() % nIndividuals){
		throw string("ERROR: length of allele count vector (") + to_string( maCounts.size() ) + string(" is not divisible by the provided number of individuals (") +
			to_string(nIndividuals) + string(") in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	if ( maCounts.empty() ){
		throw string("ERROR: empty vector of minor allele counts in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	locusSize_                 = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	const uint8_t remainderInd = static_cast<uint8_t>(locusSize_ * byteSize_ - nIndividuals_);
	const uint8_t lastByteMask = static_cast<uint8_t>(0b11111111) >> remainderInd;

	// Create a vector to store random bytes for stochastic heterozygote resolution
	vector<uint64_t> rand( locusSize_ / llWordSize_ + static_cast<bool>(locusSize_ % llWordSize_) );
	uint8_t *randBytes = reinterpret_cast<uint8_t*>( rand.data() );
	const float fNind  = static_cast<float>(nIndividuals_);

	binGenotypes_.resize(nLoci_ * locusSize_, 0);

	for (size_t jLoc = 0; jLoc < nLoci_; jLoc++) {
		// Fill the random byte vector
		for (auto &rv : rand){
			rv = rng_.ranInt();
		}
		size_t iIndiv         = 0;
		const size_t begByte  = jLoc * locusSize_;
		const size_t begIndiv = jLoc * nIndividuals_;
		for (size_t iByte = begByte; iByte < begByte + locusSize_ - 1; iByte++){                   // treat the last byte separately
			for (uint8_t iInByte = 0; iInByte < byteSize_; iInByte++){
				uint8_t curIndiv          = static_cast<uint8_t>(maCounts[begIndiv + iIndiv]);     // cramming down to one byte because I do not care what the actual value is
				curIndiv                 &= 0b10000011;                                            // mask everything in the middle
				const uint8_t missingMask = (~curIndiv) >> 7;                                      // 0b00000000 iff is missing (negative value)
				const uint8_t randMask    = (randBytes[iByte] >> iInByte) & oneBit_;               // 0b00000000 or 0b00000001 with equal chance
				uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // curIndiv == 0b00000010 curBitMask = 0b00000001; if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
				curBitMask               &= missingMask;                                           // zero it out if missing value is set
				binGenotypes_[iByte]     |= curBitMask << iInByte;
				iIndiv++;
			}
		}
		// now deal with the last byte in the individual
		for (uint8_t iRem = 0; iRem < remainderInd; iRem++){
			uint8_t curIndiv          = static_cast<uint8_t>(maCounts[begIndiv + iIndiv]);     // cramming down to one byte because I do not care what the actual value is
			curIndiv                 &= 0b10000011;                                            // mask everything in the middle
			const uint8_t missingMask = (~curIndiv) >> 7;                                      // 0b00000000 iff is missing (negative value)
			const uint8_t randMask    = (randBytes[locusSize_ - 1] >> iRem) & oneBit_;         // 0b00000000 or 0b00000001 with equal chance
			uint8_t curBitMask        = (curIndiv >> 1) ^ (curIndiv & randMask);               // curIndiv == 0b00000010 curBitMask = 0b00000001; if curIndiv == 0b00000001 or 0b00000011 (i.e. het) can be 1 or 0 with equal chance
			curBitMask               &= missingMask;                                           // zero it out if missing value is set
			binGenotypes_.back()     |= curBitMask << iRem;
			iIndiv++;
		}
		float aaCount = static_cast<float>( countSetBits_(binGenotypes_, begByte, locusSize_) ) / fNind;
		if (aaCount > 0.5){ // always want the alternative to be the minor allele
			for (size_t i = begByte; i < begByte + locusSize_; i++){
				binGenotypes_[i] = ~binGenotypes_[i];
			}
			binGenotypes_[begByte + locusSize_ - 1] &= lastByteMask; // unset the remainder bits
			aaCount = 1.0 - aaCount;
		} else {
			for (size_t i = begByte; i < begByte + locusSize_; i++){
				binGenotypes_[i] = ~binGenotypes_[i];
			}
		}
		aaf_.push_back(aaCount);
	}
}

GenoTable::GenoTable(GenoTable &&in){
	if (this != &in){
		binGenotypes_ = move(in.binGenotypes_);
		sketches_     = move(in.sketches_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;
	}
}

GenoTable& GenoTable::operator=(GenoTable &&in){
	if (this != &in){
		binGenotypes_ = move(in.binGenotypes_);
		sketches_     = move(in.sketches_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;

	}
	return *this;
}

void GenoTable::saveGenoBinary(const string &outFileName) const {
	fstream out;
	out.open(outFileName.c_str(), ios::out | ios::binary | ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), binGenotypes_.size() );
	out.close();
}

void GenoTable::makeIndividualOPH(const size_t &kSketches){
	const size_t sketchSize = nIndividuals_ / kSketches + static_cast<bool>(nIndividuals_ % kSketches);
	if (sketchSize >= emptyBinToken_){
		throw string("ERROR: Number of sketches (") + to_string(kSketches) + string(") implies sketch size (") +
			to_string(sketchSize) + string(") that is larger than ") + to_string(emptyBinToken_) +
			string(", the largest allowed value in GenoTable::makeIndividualOPH(const size_t &kSketches)");
	}
	// Calculate the actual sketch number based on the realized sketch size
	sketches_.resize(kSketches * nLoci_, emptyBinToken_);

	// generate the sequence of random integers; each column must be permuted the same
	vector<size_t> ranInts;
	size_t i = nIndividuals_;
	while (i >= 2UL){
		ranInts.push_back( rng_.ranInt() % i ); // need 0 <= j <= i, so i is actually i+1 (compared to the Wikipedia description)
		i--;
	}

	vector<uint32_t> seeds;
	seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
	for (size_t iLoc = 0; iLoc < nLoci_; iLoc++){
		// Start with a permutation
		size_t colInd = iLoc * locusSize_;
		size_t iIndiv = nIndividuals_ - 1UL; // safe b/c nIndividuals_ > 1 is checked at construction
		for (const auto &ri : ranInts){
			uint16_t firstIdx  = iIndiv % byteSize_;
			size_t firstByte   = (iIndiv / byteSize_) + colInd;
			uint16_t secondIdx = ri % byteSize_;
			size_t secondByte  = (ri / byteSize_) + colInd;
			uint16_t diff      = byteSize_ * (firstByte != secondByte); // will be 0 if the same byte is being accessed; then need to swap bits within byte

			// swapping bits within a two-byte variable
			// using the method in https://graphics.stanford.edu/~seander/bithacks.html#SwappingBitsXOR
			// if the same byte is being accessed, secondIdx is not shifted to the second byte
			// This may be affected by endianness (did not test)
			uint16_t twoBytes  = (static_cast<uint16_t>(binGenotypes_[secondByte]) << 8) | (static_cast<uint16_t>(binGenotypes_[firstByte]));
			secondIdx         += diff;
			uint16_t x         = ((twoBytes >> firstIdx) ^ (twoBytes >> secondIdx)) & 1;
			twoBytes          ^= ((x << firstIdx) | (x << secondIdx));

			memcpy(binGenotypes_.data() + firstByte, &twoBytes, sizeof(uint8_t));
			twoBytes = twoBytes >> diff;
			memcpy(binGenotypes_.data() + secondByte, &twoBytes, sizeof(uint8_t));
			iIndiv--;
		}
		// Now make the sketches
		vector<size_t> filledIndexes;                     // indexes of the non-empty sketches
		size_t iByte     = iLoc * locusSize_;
		size_t colEnd    = iByte + locusSize_;
		size_t sketchBeg = iLoc * kSketches;
		size_t iSeed     = 0;                             // index into the seed vector
		uint8_t iInByte  = 0;
		// A possible optimization is to test a whole byte for 0
		// Will test later
		for (size_t iSketch = 0; iSketch < kSketches; iSketch++){
			uint16_t firstSetBitPos = 0;
			while ( ( ( (oneBit_ << iInByte) & binGenotypes_[iByte] ) == 0 ) &&
					(firstSetBitPos < sketchSize) && (iByte != colEnd) ){
				iInByte++;
				// these are instead of an if statement
				iByte  += iInByte == byteSize_;
				iInByte = iInByte % byteSize_;
				firstSetBitPos++;
			}
			if ( (iByte < colEnd) && (firstSetBitPos < sketchSize) ){
				filledIndexes.push_back(iSketch);
				sketches_[sketchBeg + iSketch] = firstSetBitPos;

				uint16_t remainder = sketchSize - firstSetBitPos;
				uint16_t inByteMod = remainder % byteSize_;
				uint16_t inByteSum = iInByte + inByteMod;

				iByte  += remainder / byteSize_ + inByteSum / byteSize_;
				iInByte = inByteSum % byteSize_;
			}
		}
		if (filledIndexes.size() == 1){
			for (size_t iSk = 0; iSk < kSketches; iSk++){ // this will overwrite the one assigned sketch, but the wasted operation should be swamped by the rest
				sketches_[sketchBeg + iSk] = sketches_[filledIndexes[0] + sketchBeg];
			}
		} else if (filledIndexes.size() != kSketches){
			size_t emptyCount = kSketches - filledIndexes.size();
			while (emptyCount){
				for (const auto &f : filledIndexes){
					uint32_t newIdx = murMurHash_(f, seeds[iSeed]) % kSketches + sketchBeg;
					if ( sketches_[newIdx] == emptyBinToken_ ){
						sketches_[newIdx] = sketches_[f + sketchBeg];
						emptyCount--;
						break;
					}
				}
				iSeed++;
				if ( iSeed == seeds.size() ){
					seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
				}
			}
		}
	}
}

void GenoTable::allHashLD(vector<float> &LDmat) const {
	if ( (nLoci_ / 2) > ( numeric_limits<size_t>::max() / (nLoci_ - 1) ) ){ // too many loci to fit in the upper triangle
		throw string("ERROR: Number of loci (") + to_string(nLoci_) + string(") too large to calculate all by all LD in GenoTable::allHashLD(vector<float> &)");
	}
	if ( sketches_.empty() ){
		throw string("ERROR: Cannot has-based calculate LD on empty sketches");
	}
	LDmat.resize(nLoci_ * (nLoci_ - 1) / 2, 0.0);
	const size_t kSketches = sketches_.size() / nLoci_;
	const float fNind = 1.0 / static_cast<float>(kSketches);
	size_t resInd = 0;
	for (size_t iRow = 0; iRow < nLoci_; iRow++) {
		for (size_t jCol = iRow + 1; jCol < nLoci_; jCol++){
			float simVal = 0.0;
			size_t sketchRowInd = iRow * kSketches;
			size_t sketchColInd = jCol * kSketches;
			for (size_t iSk = 0; iSk < kSketches; iSk++){
				if (sketches_[sketchRowInd + iSk] == sketches_[sketchColInd + iSk]){
					simVal += 1.0;
				}
			}
			simVal *= fNind;
			//simVal -= aaf_[iRow] * aaf_[jCol]; // subtracting expected similarity
			LDmat[resInd] = simVal;
			resInd++;
		}
	}
}

void GenoTable::allJaccardLD(vector<float> &LDmat) const {
	if ( (nLoci_ / 2) > ( numeric_limits<size_t>::max() / (nLoci_ - 1) ) ){ // too many loci to fit in the upper triangle
		throw string("ERROR: Number of loci (") + to_string(nLoci_) + string(") too large to calculate all by all LD in GenoTable::allJaccardLD(vector<float> &)");
	}
	LDmat.resize(nLoci_ * (nLoci_ - 1) / 2, 0.0);
	vector<uint8_t> locus(locusSize_);
	size_t resInd = 0;
	for (size_t iRow = 0; iRow < nLoci_; iRow++) {
		for (size_t jCol = iRow + 1; jCol < nLoci_; jCol++){
			size_t rowInd = iRow * locusSize_;
			size_t colInd = jCol * locusSize_;
			for (size_t iBinLoc = 0; iBinLoc < locusSize_; iBinLoc++){
				locus[iBinLoc] = binGenotypes_[rowInd + iBinLoc] & binGenotypes_[colInd + iBinLoc];
			}
			const uint32_t uni = countSetBits_(locus);
			for (size_t iBinLoc = 0; iBinLoc < locusSize_; iBinLoc++){
				locus[iBinLoc] = binGenotypes_[rowInd + iBinLoc] | binGenotypes_[colInd + iBinLoc];
			}
			const uint32_t isect = countSetBits_(locus);
			LDmat[resInd] = static_cast<float>(uni) / static_cast<float>(isect);
			resInd++;
		}
	}
}

uint32_t GenoTable::murMurHash_(const size_t &key, const uint32_t &seed) const {
	uint32_t hash = seed;

	// body
	auto blocks = reinterpret_cast<const uint32_t *>(&key);

	for(size_t i = 0; i < nblocks_; i++){
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

uint32_t GenoTable::countSetBits_(const vector<uint8_t> &inVec) const {
	uint32_t totSet = 0;
	for (const auto &in : inVec){
		uint8_t v = in;
		for (; v; totSet++) {
			v &= v - 1;
		}
	}
	return totSet;
}

uint32_t GenoTable::countSetBits_(const vector<uint8_t> &inVec, const size_t &start, const size_t &length) const {
	uint32_t totSet = 0;
	for (size_t i = start; i < start + length; i++){
		uint8_t v = inVec[i];
		for (; v; totSet++) {
			v &= v - 1;
		}
	}
	return totSet;
}
