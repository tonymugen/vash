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

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>
#include <fstream>

#include <iostream>

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
const uint8_t GenoTable::oneBit_            = 1;
const uint8_t GenoTable::byteSize_          = 8;
const uint8_t GenoTable::llWordSize_        = 8;
const size_t GenoTable::nblocks_            = sizeof(size_t) / 4;
const uint32_t GenoTable::mmhKeyLen_        = sizeof(size_t);
const uint16_t GenoTable::emptyBinToken_    = numeric_limits<uint16_t>::max();
const uint32_t GenoTable::c1_               = 0xcc9e2d51;
const uint32_t GenoTable::c2_               = 0x1b873593;

// Constructors
GenoTable::GenoTable(const string &inputFileName, const size_t &nIndividuals, const size_t &kSketches) : nIndividuals_{nIndividuals}, nLoci_{0} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in the GenoTable(const string &, const size_t &) constructor");
	} else if (nIndividuals > numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
		throw string("ERROR: the number of individuals (") + to_string(nIndividuals) + string(") is too big to make a square relationship matrix in the GenoTable(const string &, const size_t &) constructor");
	}
	size_t nBedBytes = nIndividuals_ / 4 + static_cast<bool>(nIndividuals % 4);
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
	const uint8_t locusRemainder = nIndividuals_ % byteSize_;
	locusSize_ = nIndividuals_ / byteSize_ + static_cast<bool>(locusRemainder);
	vector<char> bedLocus(nBedBytes, 0);
	uint8_t iInByteB  = 0; // index within the current binary byte
	size_t  iBinGeno  = 0; // binGenotypes_ vector index
	array<uint8_t, llWordSize_> rand;
	uint64_t randW = rng_.ranInt();
	memcpy(rand.data(), &randW, llWordSize_);
	uint8_t iInRand   = 0;
	uint8_t iRandByte = 0;
	const float fNind = 2.0 * static_cast<float>(nIndividuals_); // double for diploids
	while (inStr) {
		inStr.read(bedLocus.data(), nBedBytes);
		vector<uint8_t> binLocus(locusSize_, 0);
		float aaCount = 0.0;
		for (const auto &b : bedLocus){
			for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
				if ( b & (oneBit_ << iInByteG) ){
					if ( b & ( oneBit_ << (iInByteG + 1) ) ){           // homozygous alternative allele
						binLocus[iBinGeno] |= oneBit_ << iInByteB;
						iInByteB++;
						aaCount += 2.0;
					} else {                                            // heterozygous
						uint8_t testBit = (oneBit_ << iInRand) & rand[iRandByte];
						if (testBit){
							binLocus[iBinGeno] |= oneBit_ << iInByteB;
						}
						aaCount += 1.0;
						iInRand++;
						if (iInRand == byteSize_){
							iInRand = 0;
							iRandByte++;
							if (iRandByte == llWordSize_){
								randW = rng_.ranInt();
								memcpy(rand.data(), &randW, llWordSize_);
								iRandByte = 0;
							}
						}
						iInByteB++;
					}
				} else {
					iInByteB++;
					// do not care what the next bit is: if the odd one is 0 (ie, homozygous or missing), the corresponding binary bit is zero
				}
				if (iInByteB == byteSize_){
					iInByteB = 0;
					iBinGeno++;
				}
			}
		}
		aaCount /= fNind;
		if (aaCount < -0.5){ // always want the alternative to be the minor allele
			for (auto &bg : binLocus){
				bg = ~bg;
				binGenotypes_.push_back(bg);
			}
			binGenotypes_.back() = binGenotypes_.back() << locusRemainder; // unset the remainder bits
			aaCount = 1.0 - aaCount;
		} else {
			for (const auto &bg : binLocus){
				binGenotypes_.push_back(bg);
			}
		}
		std::cout << aaCount << "\n";
		aaf_.push_back(aaCount);
		nLoci_++;
	}
	inStr.close();

	const size_t sketchSize = nIndividuals_ / kSketches + static_cast<bool>(nIndividuals_ % kSketches);
	if (sketchSize >= emptyBinToken_){
		throw string("ERROR: Number of sketches (") + to_string(kSketches) + string(") implies sketch size (") +
			to_string(sketchSize) + string(") that is larger than ") + to_string(emptyBinToken_) +
			string(", the largest allowed value in GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	sketchSize_ = static_cast<uint16_t>(sketchSize);
	kSketches_  = static_cast<uint16_t>( nIndividuals_ / static_cast<size_t>(sketchSize_) ) + static_cast<bool>(nIndividuals % sketchSize_);
	sketches_.resize(kSketches * nLoci_, emptyBinToken_);

	// generate the sequence of random integers; each column must be permuted the same
	vector<size_t> ranInts;
	size_t i = nIndividuals_;
	while (i >= 2UL){
		ranInts.push_back( rng_.ranInt() % i ); // need 0 <= j <= i, so i is actually i+1 (compared to the Wikipedia description)
		i--;
	}
	/*
	vector<uint32_t> seeds;
	seeds.push_back( static_cast<uint32_t>( rng_.ranInt() ) );
	for (size_t iLoc = 0; iLoc < nLoci_; iLoc++){
		makeSketches_(iLoc, ranInts, seeds);
	}
	*/
}

GenoTable::GenoTable(const vector<int8_t> &maCounts, const size_t &nIndividuals, const size_t &kSketches) : nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	} else if (nIndividuals > numeric_limits<size_t>::max() / nIndividuals ){ // a square will overflow
		throw string("ERROR: the number of individuals (") + to_string(nIndividuals) + string(") is too big to make a square relationship matrix in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	} else if (nLoci_ > numeric_limits<size_t>::max() / nLoci_ ){ // a square will overflow
		throw string("ERROR: the number of loci (") + to_string(nLoci_) + string(") is too big to make a square LD matrix in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	if (maCounts.size() % nIndividuals){
		throw string("ERROR: length of allele count vector (") + to_string( maCounts.size() ) + string(" is not divisible by the provided number of individuals (") +
			to_string(nIndividuals) + string(") in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	if ( maCounts.empty() ){
		throw string("ERROR: empty vector of minor allele counts in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	if (nIndividuals_ % 4){
		genotypes_.resize( (1 + nIndividuals_ / 4) * nLoci_, 0 );
	} else {
		genotypes_.resize( (nIndividuals_ / 4) * nLoci_, 0 );
	}

	uint8_t iInByte = 0; // index within the current byte
	size_t  iGeno   = 0; // genotype index in the genotype_ vector
	size_t  iIndv   = 0; // current individual index
	for (const auto &mac : maCounts){
		switch (mac){
			case 0:           // 00 for homozygous major allele
				{
					iInByte += 2;
					break;
				}
			case 1:           // 10 for heterozygous
				{
					genotypes_[iGeno] |= oneBit_ << iInByte;
					iInByte += 2;
					break;
				}
			case 2:           // 11 for homozygous minor allele
				{
					genotypes_[iGeno] |= oneBit_ << iInByte;
					iInByte++;
					genotypes_[iGeno] |= oneBit_ << iInByte;
					iInByte++;
					break;
				}
			case -9:          // 01 for missing genotype
				{
					iInByte++;
					genotypes_[iGeno] |= oneBit_ << iInByte;
					iInByte++;
					break;
				}
			default:
				throw string("ERROR: unknown value ") + to_string(mac) + string(" in GenoTable(const vector<int8_t> &, const size_t &) constructor");
		}
		iIndv++;
		if (iIndv == nIndividuals_){
			iInByte = 0;
			iIndv   = 0;
			iGeno++;
		} else {
			if (iInByte == byteSize_){
				iInByte = 0;
				iGeno++;
			}
		}
	}

	const size_t sketchSize = nIndividuals_ / kSketches + static_cast<bool>(nIndividuals_ % kSketches);
	if (sketchSize >= emptyBinToken_){
		throw string("ERROR: Number of sketches (") + to_string(kSketches) + string(") implies sketch size (") +
			to_string(sketchSize) + string(") that is larger than ") + to_string(emptyBinToken_) +
			string(", the largest allowed value in GenoTable(const vector<int8_t> &, const size_t &) constructor");
	}
	sketchSize_ = static_cast<uint16_t>(sketchSize);
	kSketches_  = static_cast<uint16_t>( nIndividuals_ / static_cast<size_t>(sketchSize_) ) + static_cast<bool>(nIndividuals % sketchSize_);
	locusSize_  = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	sketches_.resize(kSketches * nLoci_, emptyBinToken_);

	generateBinGeno_();
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
		makeSketches_(iLoc, ranInts, seeds);
	}
}

GenoTable::GenoTable(GenoTable &&in){
	if (this != &in){
		genotypes_    = move(in.genotypes_);
		binGenotypes_ = move(in.binGenotypes_);
		sketches_     = move(in.sketches_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;
		kSketches_    = in.kSketches_;
		sketchSize_   = in.sketchSize_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;
	}
}

GenoTable& GenoTable::operator=(GenoTable &&in){
	if (this != &in){
		genotypes_    = move(in.genotypes_);
		binGenotypes_ = move(in.binGenotypes_);
		sketches_     = move(in.sketches_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;
		kSketches_    = in.kSketches_;
		sketchSize_   = in.sketchSize_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;

	}
	return *this;
}

void GenoTable::saveGenoBed(const string &outFileName) const {
	fstream out;
	out.open(outFileName.c_str(), ios::out | ios::binary | ios::trunc);
	out.write( magicBytes_.data(), magicBytes_.size() );
	out.write( reinterpret_cast<const char*>( genotypes_.data() ), genotypes_.size() );
	out.close();
}

void GenoTable::saveGenoBinary(const string &outFileName) const {
	fstream out;
	out.open(outFileName.c_str(), ios::out | ios::binary | ios::trunc);
	out.write( reinterpret_cast<const char*>( binGenotypes_.data() ), binGenotypes_.size() );
	out.close();
}

void GenoTable::outputBits(const vector<uint8_t> &binVec, string &bitString) const {
	bitString.clear();
	uint8_t iInByte   = 0;
	size_t  iOfByte   = 0;
	size_t  iInSketch = 0;
	while ( iOfByte < binVec.size() ){
		if ( binVec[iOfByte] & (oneBit_ << iInByte) ){
			bitString += '1';
			iInByte++;
		} else {
			bitString += '0';
			iInByte++;
		}
		if (iInByte == byteSize_){
			iOfByte++;
			bitString += ' ';
			iInByte = 0;
		}
		iInSketch++;
		if (iInSketch == sketchSize_){
			bitString += "\n";
			iInSketch = 0;
		}
	}
}

void GenoTable::allSimilarity(vector<float> &LDmat) const {
	if ( (nLoci_ / 2) > ( numeric_limits<size_t>::max() / (nLoci_ - 1) ) ){ // too many loci to fit in the upper triangle
		throw string("ERROR: Number of loci (") + to_string(nLoci_) + string(") too large to calculate all by all LD");
	}
	LDmat.resize(nLoci_ * (nLoci_ - 1) / 2, 0.0);
	const float fNind = 1.0 / static_cast<float>(kSketches_);
	size_t resInd = 0;
	for (size_t iRow = 0; iRow < nLoci_; iRow++) {
		for (size_t jCol = iRow + 1; jCol < nLoci_; jCol++){
			float simVal = 0.0;
			size_t sketchRowInd = iRow * kSketches_;
			size_t sketchColInd = jCol * kSketches_;
			for (size_t iSk = 0; iSk < kSketches_; iSk++){
				if (sketches_[sketchRowInd + iSk] == sketches_[sketchColInd + iSk]){
					simVal += 1.0;
				}
			}
			simVal *= fNind;
			simVal -= aaf_[iRow] * aaf_[jCol]; // subtracting expected similarity
			LDmat[resInd] = simVal;
			resInd++;
		}
	}
}

void GenoTable::generateBinGeno_(){
	if (nIndividuals_ % 8){
		binGenotypes_.resize(nLoci_ * (1 + nIndividuals_ / 8), 0);
	} else {
		binGenotypes_.resize(nLoci_ * nIndividuals_ / 8, 0);
	}
	uint8_t iInByteB  = 0; // index within the current binary byte
	size_t  iBinGeno  = 0; // binGenotypes_ vector index
	array<uint8_t, llWordSize_> rand;
	uint64_t randW = rng_.ranInt();
	memcpy(rand.data(), &randW, llWordSize_);
	uint8_t iInRand    = 0;
	uint8_t iRandByte  = 0;
	const double dNind = 2.0 * static_cast<double>(nIndividuals_); // double for diploids
	double aaCount     = 0.0;
	for (const auto &g : genotypes_){
		for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
			if ( g & (oneBit_ << iInByteG) ){
				if ( g & ( oneBit_ << (iInByteG + 1) ) ){           // homozygous alternative allele
					binGenotypes_[iBinGeno] |= oneBit_ << iInByteB;
					iInByteB++;
					aaCount += 2.0;
				} else {                                            // heterozygous
					uint8_t testBit = (oneBit_ << iInRand) & rand[iRandByte];
					if (testBit){
						binGenotypes_[iBinGeno] |= oneBit_ << iInByteB;
					}
					aaCount += 1.0;
					iInRand++;
					if (iInRand == byteSize_){
						iInRand = 0;
						iRandByte++;
						if (iRandByte == llWordSize_){
							randW = rng_.ranInt();
							memcpy(rand.data(), &randW, llWordSize_);
							iRandByte = 0;
						}
					}
					iInByteB++;
				}
			} else {
				iInByteB++;
				// do not care what the next bit is: if the odd one is 0 (ie, homozygous or missing), the corresponding binary bit is zero
			}
			if (iInByteB == byteSize_){
				iInByteB = 0;
				iBinGeno++;
				if ( (iBinGeno % locusSize_) == 0){ // finished a locus
					aaf_.push_back(aaCount / dNind);
					aaCount = 0.0;
				}
			}
		}
	}
}

void GenoTable::makeSketches_(const size_t &locusIdx, const vector<size_t> &ranInts, vector<uint32_t> &seeds){
	// Start with a permutation
	size_t colInd = locusIdx * locusSize_;
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
	size_t iByte     = locusIdx * locusSize_;
	size_t colEnd    = iByte + locusSize_;
	size_t sketchBeg = locusIdx * kSketches_;
	size_t iSeed     = 0;                             // index into the seed vector
	uint8_t iInByte  = 0;
	// A possible optimization is to test a whole byte for 0
	// Will test later
	for (size_t iSketch = 0; iSketch < kSketches_; iSketch++){
		uint16_t firstSetBitPos = 0;
		while ( ( ( (oneBit_ << iInByte) & binGenotypes_[iByte] ) == 0 ) &&
				(firstSetBitPos < sketchSize_) ){
			iInByte++;
			if (iInByte == byteSize_){
				iInByte = 0;
				iByte++;
				if (iByte == colEnd){
					break;
				}
			}
			firstSetBitPos++;
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
		for (size_t iSk = 0; iSk < kSketches_; iSk++){ // this will overwrite the one assigned sketch, but the wasted operation should be swamped by the rest
			sketches_[sketchBeg + iSk] = sketches_[filledIndexes[0] + sketchBeg];
		}
	} else if (filledIndexes.size() != kSketches_){
		size_t emptyCount = kSketches_ - filledIndexes.size();
		while (emptyCount){
			for (const auto &f : filledIndexes){
				uint32_t newIdx = murMurHash_(f, seeds[iSeed]) % kSketches_ + sketchBeg;
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
