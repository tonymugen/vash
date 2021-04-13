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

using namespace BayesicSpace;

const array<uint8_t, 3> GenoTable::magicBytes_ = {0x6c, 0x1b, 0x01};

// Constructors
GenoTable::GenoTable(const string &inputFileName) {}

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
				iInByte += 2;
				return;
			case 1:           // 10 for heterozygous
				genotypes_[iGeno] |= static_cast<uint8_t>(1) << iInByte;
				iInByte += 2;
				return;
			case 2:           // 11 for homozygous minor allele
				genotypes_[iGeno] |= static_cast<uint8_t>(1) << iInByte;
				iInByte++;
				genotypes_[iGeno] |= static_cast<uint8_t>(1) << iInByte;
				iInByte++;
				return;
			case -9:          // 01 for missing genotype
				iInByte++;
				genotypes_[iGeno] |= static_cast<uint8_t>(1) << iInByte;
				iInByte++;
				return;
			default:
				throw string("ERROR: unknown value ") + to_string(mac) + string(" in GenoTable(const vector<int8_t> &, const size_t &) constructor");
				return;
		}
		iIndv++;
		iIndv %= nIndividuals_;
		if (iIndv == 0){
			iInByte = 0;
			iGeno++;
		} else {
			iInByte %= 8;
			if (iInByte == 0){
				iGeno++;
			}
		}
	}
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

GenoTable& GenoTable::operator=(GenoTable &&in) {
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
	out.write( reinterpret_cast<const char*>( magicBytes_.data() ), magicBytes_.size() );
	out.write( reinterpret_cast<const char*>( genotypes_.data() ), genotypes_.size() );
	out.close();
}
