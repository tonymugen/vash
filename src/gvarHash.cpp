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
#include <bits/stdint-uintn.h>
#include <cstddef>
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>

#include <iostream>

#include "gvarHash.hpp"

using std::vector;
using std::array;
using std::string;
using std::to_string;
using std::move;
using std::fstream;
using std::ios;
using std::streampos;

using namespace BayesicSpace;

const array<char, 3> GenoTable::magicBytes_ = {0x6c, 0x1b, 0x01};
const uint8_t GenoTable::oneBit_            = 1;
const uint8_t GenoTable::byteSize_          = 8;
const uint8_t GenoTable::llWordSize_        = 8;
const uint8_t GenoTable::emptyBinToken_     = 0xff;
const uint32_t GenoTable::c1_               = 0xcc9e2d51;
const uint32_t GenoTable::c2_               = 0x1b873593;

// Constructors
GenoTable::GenoTable(const string &inputFileName, const size_t &nIndividuals, const size_t &kSketches) : nIndividuals_{nIndividuals} {
	if (nIndividuals <= 1){
		throw string("ERROR: number of individuals must be greater than 1 in the GenoTable(const string &, const size_t &) constructor");
	}
	size_t nBytes;
	if (nIndividuals_ % 4){
		nBytes = 1 + nIndividuals_ / 4;
	} else {
		nBytes = nIndividuals_ / 4;
	}
	fstream inStr;
	inStr.open(inputFileName.c_str(), ios::in | ios::binary);
	inStr.seekg(0, inStr.end);
	streampos inFileSize = inStr.tellg();
	inStr.seekg(0, inStr.beg);
	if ( inFileSize <= magicBytes_.size() ){
		throw string("ERROR: input .bed file has no genotypes in the GenoTable(const string &, const size_t &) constructor");
	}
	inFileSize -= magicBytes_.size();
	if (static_cast<size_t>(inFileSize) % nBytes){
		throw string("ERROR: The .bed file size not evenly divisible by the number of individuals provided (") + to_string(nIndividuals_) + string(")");
	}
	char magicBuf[magicBytes_.size()]{};
	inStr.read( magicBuf, magicBytes_.size() );
	if (magicBuf[0] != magicBytes_[0]){
		throw string("ERROR: first magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	} else if (magicBuf[1] != magicBytes_[1]){
		throw string("ERROR: second magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	} else if (magicBuf[2] != magicBytes_[2]){
		throw string("ERROR: third magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	}
	genotypes_.resize(static_cast<size_t>(inFileSize));
	inStr.read(reinterpret_cast<char *>( genotypes_.data() ), inFileSize);
	inStr.close();
	nLoci_ = static_cast<size_t>(inFileSize) / nBytes;

	sketchSize_ = static_cast<uint16_t>(nIndividuals_ / kSketches) + static_cast<bool>(nIndividuals_ % kSketches);
	kSketches_  = static_cast<uint16_t>( nIndividuals_ / static_cast<size_t>(sketchSize_) ) + static_cast<bool>(nIndividuals % sketchSize_);
	locusSize_  = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);

	generateBinGeno_();
	for (size_t iLoc = 0; iLoc < nLoci_; iLoc++) {
		permuteIndv_(iLoc);
		makeSketches_(iLoc);
	}
}

GenoTable::GenoTable(const vector<int8_t> &maCounts, const size_t &nIndividuals, const size_t &kSketches) : nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals} {
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

	sketchSize_ = static_cast<uint16_t>(nIndividuals_ / kSketches) + static_cast<bool>(nIndividuals_ % kSketches);
	kSketches_  = static_cast<uint16_t>( nIndividuals_ / static_cast<size_t>(sketchSize_) ) + static_cast<bool>(nIndividuals % sketchSize_);
	locusSize_  = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);

	generateBinGeno_();
	permuteIndv_(0);
	makeSketches_(0);
	string binStr;
	vector<uint8_t> subs(binGenotypes_.begin(), binGenotypes_.begin() + 25);
	outputBits(subs, binStr);
	std::cout << binStr << "\n";
	std::cout << kSketches_ << "\n";
	std::cout << sketchSize_ << "\n";
	std::cout << locusSize_ << "\n";
	//const uint32_t seed = 1;
	//uint32_t hash = murMurHash_(binGenotypes_, 50, 25, seed);
	//std::cout << hash << "\n";
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

/*
void GenoTable::outputBits(string &bitString){
	bitString.clear();
	uint8_t iInByte = 0;
	size_t  iOfByte = 0;
	for (size_t iInd = 0; iInd < nIndividuals_; iInd++){
		if (binGenotypes_[iOfByte] & (oneBit_ << iInByte) ){
			bitString += '1';
			iInByte++;
		} else {
			bitString += '0';
			iInByte++;
		}
		if (iInByte == byteSize_){
			iOfByte++;
			if (iOfByte % llWordSize_){
				bitString += ' ';
			} else {
				bitString += "\n";
			}
			iInByte = 0;
		}
	}
}
*/
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
	uint8_t iInRand   = 0;
	uint8_t iRandByte = 0;
	for (const auto &g : genotypes_){
		for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2){
			if ( g & (oneBit_ << iInByteG) ){
				if ( g & ( oneBit_ << (iInByteG + 1) ) ){           // homozygous minor allele
					binGenotypes_[iBinGeno] |= oneBit_ << iInByteB;
					iInByteB++;
				} else {                                            // heterozygous
					uint8_t testBit = (oneBit_ << iInRand) & rand[iRandByte];
					if (testBit){
						binGenotypes_[iBinGeno] |= oneBit_ << iInByteB;
					}
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
				// do not care what the next bit is: if the odd one is 0 (ie, homozygous or missing), the binary bit is zero
			}
			if (iInByteB == byteSize_){
				iInByteB = 0;
				iBinGeno++;
			}
		}
	}
}

void GenoTable::permuteIndv_(const size_t &locusIdx){
	// generate the sequence of random integers; each column must be permuted the same
	vector<size_t> ranInts;
	size_t i = nIndividuals_;
	while (i >= 2UL){
		ranInts.push_back( rng_.ranInt() % i ); // need 0 <= j <= i, so i is actually i+1 (compared to the Wikipedia description)
		i--;
	}
	size_t locusSize = nIndividuals_ / byteSize_ + static_cast<bool>(nIndividuals_ % byteSize_);
	size_t colInd = locusIdx * locusSize;
	size_t iIndiv = nIndividuals_ - 1UL; // safe b/c nIndividuals_ > 1 is checked at construction
	for (const auto &ri : ranInts){
		uint16_t firstIdx  = iIndiv % byteSize_;
		size_t firstByte   = (iIndiv / byteSize_) + colInd;
		uint16_t secondIdx = ri % byteSize_;
		size_t secondByte  = (ri / byteSize_) + colInd;
		uint16_t diff      = 8 * (firstByte != secondByte); // will be 0 if the same byte is being accessed; then need to swap bits within byte

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
}

void GenoTable::makeSketches_(const size_t &locusIdx){
	size_t colInd = locusIdx * locusSize_;
	size_t colEnd = colInd + locusSize_;
	vector<uint16_t> nonEmptySketches;                                                           // vector with non-empty sketches
	vector<size_t> emptyIndexes;                                                                 // indexes of the empty sketches in the sketches_ vector
	uint32_t seed           = 1; // TODO: this is just temporary
	size_t iByte            = colInd;
	uint8_t iInByte         = 0;
	for (size_t iSketch = 0; iSketch < kSketches_; iSketch++) {
		uint16_t firstSetBitPos = 0;
		while ( (firstSetBitPos < sketchSize_) && (iByte < colEnd) ){
			if ( (iInByte == 0) && (binGenotypes_[iByte] == 0) && ( (sketchSize_ - firstSetBitPos) < byteSize_) ){
				iByte++;
				firstSetBitPos += byteSize_;
			}
			if ( (oneBit_ << iInByte) & binGenotypes_[iByte] ){
				nonEmptySketches.push_back(firstSetBitPos);
				iInByte++;
				if (iInByte == byteSize_){
					iInByte = 0;
					iByte++;
				}
				break;
			}
			firstSetBitPos++;
			iInByte++;
			if (iInByte == byteSize_){
				iInByte = 0;
				iByte++;
			}
		}
		if (firstSetBitPos == sketchSize_){
			emptyIndexes.push_back(iSketch);
		}
	}
	/*
	while (iSketch < kSketches_){
		if ( (iInByte == 0) && (binGenotypes_[iByte] == 0) ){ // the whole next byte is zero, no reason to iterate through it
			firstSetBitPos += byteSize_;
			if (firstSetBitPos >= sketchSize_){
				firstSetBitPos -= sketchSize_;
				emptyIndexes.push_back(iSketch);
				iSketch++;
			}
			iByte++;
			continue;
		}
		while ( ( ( (oneBit_ << iInByte) & binGenotypes_[iByte] ) == 0 ) && (iSketch < kSketches_) ){
			iInByte++;
			firstSetBitPos++;
			if (firstSetBitPos == sketchSize_){
				firstSetBitPos = 0;
				emptyIndexes.push_back(iSketch);
				iSketch++;
			}
			if (iInByte == byteSize_){
				iInByte = 0;
				iByte++;
			}
		}
		if (firstSetBitPos){
			nonEmptySketches.push_back(firstSetBitPos);
		}
		firstSetBitPos = 0;
		iSketch++;
	}
	*/
	/*
	if (nonEmptySketches.size() == 1){
		for (const auto &e : emptyIndexes){
			sketches_[e] = nonEmptySketches[0]; // TODO: will need to be fixed for loci after the 1st
		}
	} else if (nonEmptySketches.size() > 1){
		size_t remainingEmpties = emptyIndexes.size();
		while (remainingEmpties){
			for (const auto &neb : nonEmptyBins){
				vector<uint8_t> key(sizeof(neb), 0);
				memcpy(key.data(), &neb, sizeof(neb));
				uint32_t hash = murMurHash_(key, 0, key.size(), seed);
				hash %= nBytes; // TODO: need k; check if this is right
			}
		}
	}
	*/
	std::cout << "empty indexes:\n";
	for (const auto &e : emptyIndexes){
		std::cout << e << " ";
	}
	std::cout << "\nnon-empty sketches:\n";
	for (const auto &s : nonEmptySketches){
		std::cout << static_cast<uint16_t>(s) << " ";
	}
	std::cout << "\nsketches:\n";
	/*
	for (const auto &ss : sketches_){
		std::cout << static_cast<uint16_t>(ss) << " ";
	}
	std::cout << "\n";
	*/
}

uint32_t GenoTable::murMurHash_(const vector<uint8_t> &key, const size_t &start, const size_t &size, const uint32_t &seed) const {
	const size_t nblocks = size / 4;

	uint32_t hash = seed;

	// body
	auto blocks = reinterpret_cast<const uint32_t *>(key.data() + start);

	for(size_t i = 0; i < nblocks; i++){
		uint32_t k1 = blocks[i];

		k1 *= c1_;
		k1 = (k1 << 15) | (k1 >> 17);
		k1 *= c2_;

		hash ^= k1;
		hash  = (hash << 13) | (hash >> 19);
		hash  = hash * 5 + 0xe6546b64;
	}

	// tail
	const uint8_t * tail = key.data() + start + nblocks * 4;
	uint32_t k1 = 0;

	switch(size & 3){
		case 3:
			k1 ^= tail[2] << 16;
		case 2:
			k1 ^= tail[1] << 8;
		case 1:
			k1   ^= tail[0];
			k1   *= c1_;
			k1    = (k1 << 15) | (k1 >> 17);
			k1   *= c2_;
			hash ^= k1;
	};

	// finalize
	hash ^= static_cast<uint32_t>(size);
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;

	return hash;
}
