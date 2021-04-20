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
#include <cstring>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <fstream>

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
const uint8_t GenoTable::oneBit_            = static_cast<uint8_t>(1);
const uint8_t GenoTable::byteSize_          = static_cast<uint8_t>(8);
const uint8_t GenoTable::llWordSize_        = static_cast<uint8_t>(8);

// Constructors
GenoTable::GenoTable(const string &inputFileName, const size_t &nIndividuals) : nIndividuals_{nIndividuals} {
	if (nIndividuals == 0){
		throw string("ERROR: number of individuals is 0 in the GenoTable(const string &, const size_t &) constructor");
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
	} else if (magicBuf[1] != magicBytes_[1]) {
		throw string("ERROR: second magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	} else if (magicBuf[2] != magicBytes_[2]) {
		throw string("ERROR: third magic byte in input .bed file is not the expected value in the GenoTable(const string &, const size_t &) constructor");
	}
	genotypes_.resize(static_cast<size_t>(inFileSize));
	inStr.read(reinterpret_cast<char *>( genotypes_.data() ), inFileSize);
	inStr.close();
	nLoci_ = static_cast<size_t>(inFileSize) / nBytes;
	generateBinGeno_();
}

GenoTable::GenoTable(const vector<int8_t> &maCounts, const size_t &nIndividuals) : nIndividuals_{nIndividuals}, nLoci_{maCounts.size() / nIndividuals} {
	if (nIndividuals == 0){
		throw string("ERROR: number of individuals is 0 in the GenoTable(const vector<int8_t> &, const size_t &) constructor");
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
		switch (mac) {
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
		iIndv %= nIndividuals_;
		if (iIndv == 0){
			iInByte = 0;
			iGeno++;
		} else {
			iInByte %= 8;      // TODO: replace with testing iInByte == 8
			if (iInByte == 0){
				iGeno++;
			}
		}
	}
	generateBinGeno_();
}

GenoTable::GenoTable(GenoTable &&in){
	if (this != &in){
		genotypes_    = move(in.genotypes_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;

		in.nIndividuals_ = 0;
		in.nLoci_        = 0;
	}
}

GenoTable& GenoTable::operator=(GenoTable &&in){
	if (this != &in){
		genotypes_    = move(in.genotypes_);
		nIndividuals_ = in.nIndividuals_;
		nLoci_        = in.nLoci_;

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

void GenoTable::generateBinGeno_(){
	if (nIndividuals_ % 8){
		binGenotypes_.resize(nLoci_ * (1 + nIndividuals_) / 8, 0);
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
		for (uint8_t iInByteG = 0; iInByteG < byteSize_; iInByteG += 2) {
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
