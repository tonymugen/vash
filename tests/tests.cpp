#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>
#include <array>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <set>
#include <popcntintrin.h>

#include "gvarHash.hpp"
#include "vashFunctions.hpp"
#include "similarityMatrix.hpp"
#include "vashLogging.hpp"

#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

// Number of times tests of random events will be run
static constexpr uint16_t N_RAN_ITERATIONS{10};
// precision for float comparisons
static constexpr float FPREC{1e-4F};

namespace {
	// Write a minimal .bed file containing only homozygous (codes 00 and 11) and missing
	// (01) genotypes -- NO heterozygotes (code 10). Heterozygotes are the only source of
	// randomness in binarization, so a het-free file binarizes deterministically. With
	// identicalLoci == false the genotype of (locus, individual) varies with both indices;
	// with identicalLoci == true every locus is the same (so every OPH-Jaccard pair is
	// exactly 1.0, and any locus the chunked reader fails to write stands out).
	// NOLINTNEXTLINE(bugprone-easily-swappable-parameters) two dimension counts; this is test-only
	void writeHetFreeBed(const std::string &fileName, const size_t nIndividuals, const size_t nLoci, const bool identicalLoci) {
		std::fstream out(fileName, std::ios::out | std::ios::binary | std::ios::trunc);
		const std::array<char, 3> magicBytes{0x6c, 0x1b, 0x01};
		out.write( magicBytes.data(), static_cast<std::streamsize>( magicBytes.size() ) );
		const size_t bytesPerLocus{(nIndividuals + 3) / 4};
		const std::array<uint8_t, 3> codes{0b00, 0b01, 0b11};   // hom-major, missing, hom-minor (no het 0b10)
		std::vector<char> locusBytes(bytesPerLocus, 0);
		for (size_t iLocus = 0; iLocus < nLoci; ++iLocus) {
			std::fill( locusBytes.begin(), locusBytes.end(), static_cast<char>(0) );
			for (size_t iIndiv = 0; iIndiv < nIndividuals; ++iIndiv) {
				const size_t codeIdx{(identicalLoci ? iIndiv : iLocus + iIndiv) % codes.size()};
				const uint8_t code{codes.at(codeIdx)};
				const auto shifted{static_cast<uint8_t>( code << ((iIndiv % 4U) * 2U) )};
				locusBytes.at(iIndiv / 4) = static_cast<char>( static_cast<uint8_t>( locusBytes.at(iIndiv / 4) ) | shifted );
			}
			out.write( locusBytes.data(), static_cast<std::streamsize>(bytesPerLocus) );
		}
		out.close();
	}
} // anonymous namespace

TEST_CASE("Can count set bits correctly", "[countSetBits]") {
	constexpr uint16_t oneWord{0b11001110'01101001};
	constexpr uint16_t wCorrectCount{9};
	const uint16_t wordBitCount{BayesicSpace::countSetBits(oneWord)};
	REQUIRE(wordBitCount == wCorrectCount);
	constexpr size_t byteArraySize{14};
	constexpr std::array<uint8_t, byteArraySize> byteArray{
		0b11110111, 0b10100111, 0b01011000, 0b11111001, 0b11100110, 0b11111100, 0b00111111,
		0b01101011, 0b01011111, 0b01001001, 0b11001001, 0b11100010, 0b11010101, 0b01110001
	};
	constexpr uint64_t correctFullVectorCount{69};
	std::vector<uint8_t> byteVector( byteArray.begin(), byteArray.end() );
	const uint64_t fullVectorCount{BayesicSpace::countSetBits(byteVector)};
	REQUIRE(fullVectorCount == correctFullVectorCount);
	constexpr BayesicSpace::LocationWithLength byteVectorWindow{2, 4};
	constexpr uint64_t correctVectorWindowCount{20};
	const uint64_t vectorWindowCount{BayesicSpace::countSetBits(byteVector, byteVectorWindow)};
	REQUIRE(vectorWindowCount == correctVectorWindowCount);
}

TEST_CASE("Available RAM query returns a usable value", "[getAvailableRAM]") {
	// returns MemAvailable in bytes, or a non-zero fallback when /proc/meminfo is unavailable
	const size_t availableRAM{BayesicSpace::getAvailableRAM()};
	REQUIRE(availableRAM > 0);
	// the result is stable across repeated calls (allowing for normal small fluctuations)
	const size_t availableRAMagain{BayesicSpace::getAvailableRAM()};
	REQUIRE(availableRAMagain > 0);
}

TEST_CASE("MurMurHash works properly", "[MurMurHash]") {
	constexpr uint32_t mmHashSeed{2153025618};
	constexpr std::array<uint32_t, BayesicSpace::SIZE_OF_SIZET> arrayKey{335636695, 4242517348};
	constexpr uint32_t correctArrayHash{2730141477};
	const uint32_t arrayMMhash{BayesicSpace::murMurHash(arrayKey, mmHashSeed)};
	constexpr std::array<size_t, 11> idxArray{
		2437, 2444, 41116, 42353,
		45949, 58374, 75248, 80113,
		93649, 98640, 99638
	};
	std::vector<size_t> idxVector( idxArray.begin(), idxArray.end() );
	constexpr uint32_t correctIdxVectorHash{2643649892};
	const uint32_t idxVectorHash{BayesicSpace::murMurHash(idxVector, mmHashSeed)};
	constexpr std::array<uint32_t, 11> u32Array{
		2437, 2444, 41116, 42353,
		45949, 58374, 75248, 80113,
		93649, 98640, 99638
	};
	std::vector<uint32_t> u32Vector( u32Array.begin(), u32Array.end() );
	// the uint32_t and size_t overloads hash different byte widths, so their results differ by design
	constexpr uint32_t correctU32VectorHash{1671617805};
	const uint32_t u32VectorHash{BayesicSpace::murMurHash(u32Vector, mmHashSeed)};
	constexpr std::array<uint16_t, 13> array16bit{
		1256, 2117, 2866, 7434,
		11737, 16256, 22236, 39883,
		40023, 46299, 58123, 58167, 62187
	};
	std::vector<uint16_t> vector16bit( array16bit.begin(), array16bit.end() );
	constexpr BayesicSpace::LocationWithLength keyWindow{4, 5};
	constexpr uint32_t correct16bitHash{3760365877};
	const uint32_t v16bitHash{BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed)};
	constexpr BayesicSpace::LocationWithLength wholeKeySpan{0, array16bit.size()};
	constexpr uint32_t correctAll16bitHash{2280422248};
	const uint32_t all16bitHash{BayesicSpace::murMurHash(vector16bit, wholeKeySpan, mmHashSeed)};
	SECTION("MurMurHash correctness tests") {
		REQUIRE(arrayMMhash   == correctArrayHash);
		REQUIRE(idxVectorHash == correctIdxVectorHash);
		REQUIRE(u32VectorHash == correctU32VectorHash);
		REQUIRE(v16bitHash    == correct16bitHash);
		REQUIRE(all16bitHash  == correctAll16bitHash);
	}
	SECTION("MurMurHash sensitivity tests") {
		constexpr uint32_t mmHashSeed2{2153025619};
		constexpr std::array<uint32_t, BayesicSpace::SIZE_OF_SIZET> arrayKey2{335636696, 4242517348};
		REQUIRE(BayesicSpace::murMurHash(arrayKey, mmHashSeed2)  != correctArrayHash);
		REQUIRE(BayesicSpace::murMurHash(arrayKey2, mmHashSeed)  != correctArrayHash);
		REQUIRE(BayesicSpace::murMurHash(idxVector, mmHashSeed2) != correctIdxVectorHash);
		idxVector.at(1)++;
		REQUIRE(BayesicSpace::murMurHash(idxVector, mmHashSeed)  != correctIdxVectorHash);
		REQUIRE(BayesicSpace::murMurHash(u32Vector, mmHashSeed2) != correctU32VectorHash);
		u32Vector.at(1)++;
		REQUIRE(BayesicSpace::murMurHash(u32Vector, mmHashSeed)  != correctU32VectorHash);
		vector16bit.at(2)++;
		REQUIRE(BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed)  == correct16bitHash);
		REQUIRE(BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed2) != correct16bitHash);
		vector16bit.at(keyWindow.start + 2)++;
		REQUIRE(BayesicSpace::murMurHash(vector16bit, keyWindow, mmHashSeed)  != correct16bitHash);
	}
}

TEST_CASE(".bed related file and data parsing works", "[bedData]") {
	SECTION("Thread ranges") {
		constexpr BayesicSpace::CountAndSize threadSizes{4, 13};
		const std::vector< std::pair<size_t, size_t> > correctRanges{ {0, 13}, {13, 26}, {26, 39}, {39, 52} };
		const std::vector< std::pair<size_t, size_t> > threadRanges{BayesicSpace::makeThreadRanges(threadSizes)};
		REQUIRE(threadRanges.size() == threadSizes.count);
		REQUIRE(std::all_of(
				threadRanges.cbegin(),
				threadRanges.cend(),
				[](const std::pair<size_t, size_t> &eachRange) {return eachRange.first <= eachRange.second;}
			)
		);
		REQUIRE(std::equal(
				threadRanges.cbegin(),
				threadRanges.cend(),
				correctRanges.cbegin(),
				[](const std::pair<size_t, size_t> &pair1, const std::pair<size_t, size_t>&pair2) {
					return (pair1.first == pair2.first) && (pair1.second == pair2.second);
				}
			)
		);

		constexpr size_t nVecElements{35};
		constexpr size_t nChunks{4};
		constexpr std::array<size_t, nChunks> correctChunkSizes{9, 9, 9, 8};
		const std::vector<size_t> chunkSizes{BayesicSpace::makeChunkSizes(nVecElements, nChunks)};
		REQUIRE(std::equal(
				chunkSizes.cbegin(),
				chunkSizes.cend(),
				correctChunkSizes.cbegin()
			) 
		);
		constexpr size_t smallNelements{3};
		const std::vector<size_t> smallChunkSizes{BayesicSpace::makeChunkSizes(smallNelements, nChunks)};
		constexpr std::array<size_t, nChunks> correctSmallChunkSizes{1, 1, 1, 0};
		REQUIRE(std::equal(
				smallChunkSizes.cbegin(),
				smallChunkSizes.cend(),
				correctSmallChunkSizes.cbegin()
			)
		);

		constexpr std::array<uint32_t, nChunks> correctRowStarts{1, 4, 6, 7};
		constexpr std::array<uint32_t, nChunks> correctRowEnds{4, 6, 7, 8};
		constexpr std::array<uint32_t, nChunks> correctColStarts{0, 3, 3, 6};
		constexpr std::array<uint32_t, nChunks> correctColEnds{3, 3, 6, 7};
		std::vector< std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> > correctRowColPairs;
		size_t iChunk{0};
		while (iChunk < nChunks) {
			std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> tmpPair;
			tmpPair.first.iRow  = correctRowStarts.at(iChunk);
			tmpPair.first.jCol  = correctColStarts.at(iChunk);
			tmpPair.second.iRow = correctRowEnds.at(iChunk);
			tmpPair.second.jCol = correctColEnds.at(iChunk);
			correctRowColPairs.emplace_back(tmpPair);
			++iChunk;
		}
		BayesicSpace::LocationWithLength startAndNelements{};
		startAndNelements.start  = 0;
		startAndNelements.length = nVecElements;
		std::vector< std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> > rowColPairs{BayesicSpace::makeChunkRanges(startAndNelements, nChunks)};
		REQUIRE(std::equal(
				rowColPairs.cbegin(),
				rowColPairs.cend(),
				correctRowColPairs.cbegin(),
				[](const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair1,
							const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair2) {
					return  (pair1.first.iRow  == pair2.first.iRow) &&
							(pair1.first.jCol  == pair2.first.jCol) &&
							(pair1.second.iRow == pair2.second.iRow) &&
							(pair1.second.jCol == pair2.second.jCol);
				}
			)
		);
		// makeChunkRanges clamps the requested chunk count to the span length, so asking for more chunks
		// than elements yields one range per element (no empty trailing ranges) rather than nChunks ranges.
		BayesicSpace::LocationWithLength smallStartAndNelements{};
		smallStartAndNelements.start  = 0;
		smallStartAndNelements.length = smallNelements;
		std::vector< std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> > smallRowColPairs{BayesicSpace::makeChunkRanges(smallStartAndNelements, nChunks)};
		REQUIRE(smallRowColPairs.size() == smallNelements);
		constexpr std::array<uint32_t, smallNelements> smallCorrectRowStarts{1, 2, 2};
		constexpr std::array<uint32_t, smallNelements> smallCorrectRowEnds{2, 2, 3};
		constexpr std::array<uint32_t, smallNelements> smallCorrectColStarts{0, 0, 1};
		constexpr std::array<uint32_t, smallNelements> smallCorrectColEnds{0, 1, 0};
		std::vector< std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> > smallCorrectRowColPairs;
		iChunk = 0;
		while (iChunk < smallNelements) {
			std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> tmpPair;
			tmpPair.first.iRow  = smallCorrectRowStarts.at(iChunk);
			tmpPair.first.jCol  = smallCorrectColStarts.at(iChunk);
			tmpPair.second.iRow = smallCorrectRowEnds.at(iChunk);
			tmpPair.second.jCol = smallCorrectColEnds.at(iChunk);
			smallCorrectRowColPairs.emplace_back(tmpPair);
			++iChunk;
		}
		REQUIRE(std::equal(
				smallRowColPairs.cbegin(),
				smallRowColPairs.cend(),
				smallCorrectRowColPairs.cbegin(),
				[](const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair1,
							const std::pair<BayesicSpace::RowColIdx, BayesicSpace::RowColIdx> &pair2) {
					return  (pair1.first.iRow  == pair2.first.iRow) &&
							(pair1.first.jCol  == pair2.first.jCol) &&
							(pair1.second.iRow == pair2.second.iRow) &&
							(pair1.second.jCol == pair2.second.jCol);
				}
			)
		);
	}

	SECTION("Magic byte testing") {
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> correctBedBytes{0x6c, 0x1b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes1{ 0x6d, 0x1b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes2{ 0x6c, 0x0b, 0x01};
		constexpr std::array<char, BayesicSpace::N_BED_TEST_BYTES> wrongBedBytes3{ 0x6c, 0x1b, 0x11};
		REQUIRE_NOTHROW( BayesicSpace::testBedMagicBytes(correctBedBytes) );
		REQUIRE_THROWS_WITH(BayesicSpace::testBedMagicBytes(wrongBedBytes1), Catch::Matchers::StartsWith("ERROR: first magic byte in input .bed file") );
		REQUIRE_THROWS_WITH(BayesicSpace::testBedMagicBytes(wrongBedBytes2), Catch::Matchers::StartsWith("ERROR: second magic byte in input .bed file") );
		REQUIRE_THROWS_WITH(BayesicSpace::testBedMagicBytes(wrongBedBytes3), Catch::Matchers::StartsWith("ERROR: third magic byte in input .bed file") );
	}

	SECTION(".bim file reading") {
		const std::string bimFileName("../tests/ind197_397.bim");
		std::vector<std::string> locusNames{BayesicSpace::getLocusNames(bimFileName)};
		REQUIRE(locusNames.at(1)  == "14155618");
		REQUIRE(locusNames.back() == "14168708");
	}

	SECTION("Binarization and similarity groups") {
		// binarization
		constexpr size_t nBitsInByte{8};
		constexpr size_t nIndividuals{17};
		constexpr size_t nIndivPerBedByte{4};
		constexpr size_t nIndivPerBinByte{8};
		constexpr size_t nBedBytes{5};
		constexpr size_t nBinBytes{3};
		constexpr uint64_t correctMinUnion{6};
		constexpr std::array<uint8_t, nBedBytes> bedBytes{0b11001100, 0b00011011, 0b11001100, 0b00111001, 0b00000011};
		constexpr std::array<int, nIndividuals> macArray{0, 2, 0, 2, 2, 1, -9, 0, 0, 2, 0, 2, 1, -9, 2, 0, 2};
		BayesicSpace::LocationWithLength bedWindow{0, bedBytes.size()};
		BayesicSpace::LocationWithLength binWindow{0, nBinBytes};
		for (uint16_t iRanIt = 0; iRanIt < N_RAN_ITERATIONS; ++iRanIt) {
			std::vector<char> bedByteVec{bedBytes.cbegin(), bedBytes.cend()};
			std::vector<uint8_t> binBytesBed(nBinBytes, 0);
			BayesicSpace::binarizeBedLocus(bedWindow, bedByteVec, nIndividuals, binWindow, binBytesBed);
			REQUIRE(nIndividuals >= BayesicSpace::countSetBits(binBytesBed) * 2);
			std::vector<int> macVec{macArray.cbegin(), macArray.cend()};
			std::vector<uint8_t> binBytesMAC(nBinBytes, 0);
			BayesicSpace::binarizeMacLocus(macVec, binWindow, binBytesMAC);
			REQUIRE(nIndividuals >= BayesicSpace::countSetBits(binBytesMAC) * 2);
			uint32_t binBed{0};
			memcpy( &binBed, binBytesBed.data(), binBytesBed.size() );
			uint32_t binMAC{0};
			memcpy( &binMAC, binBytesMAC.data(), binBytesMAC.size() );
			const auto macBedUnion{binBed & binMAC};
			REQUIRE(_mm_popcnt_u32(macBedUnion) >= correctMinUnion);
		}

		// makeGroupRanges tests
		std::vector<BayesicSpace::HashGroup> groups;
		constexpr std::array<size_t, 3> groupSizes{7, 5, 11};
		constexpr size_t correctVGsize{86};
		groups.reserve( groupSizes.size() );
		size_t gStart{0};
		for (const auto &iGrpSize : groupSizes) {
			BayesicSpace::HashGroup currGrp;
			currGrp.locusIndexes.resize(iGrpSize);
			std::iota(currGrp.locusIndexes.begin(), currGrp.locusIndexes.end(), gStart);
			gStart                  += iGrpSize * (iGrpSize - 1) / 2;
			currGrp.cumulativeNpairs = gStart;
			groups.emplace_back(currGrp);
		}
		constexpr size_t groupStart{2UL};
		constexpr size_t testChunkSize{12UL};
		constexpr size_t correctSecondPC{4UL};
		const BayesicSpace::HashGroupItPairCount testStartCount{groupStart, std::next( groups.cbegin() )};
		const auto testGrpRange{BayesicSpace::makeGroupRanges(groups, testStartCount, testChunkSize)};
		REQUIRE(testGrpRange.first.pairCount == groupStart);
		REQUIRE( testGrpRange.first.hgIterator == std::next( groups.cbegin() ) );
		REQUIRE(testGrpRange.second.pairCount == correctSecondPC);
		REQUIRE( testGrpRange.second.hgIterator == std::next(testGrpRange.first.hgIterator) );
	}
}

TEST_CASE("Command line parsing works", "[commandLine]") {
	SECTION("parseCL flag extraction") {
		// a flag followed by a value, a value-less flag whose value is the next "--" token,
		// and a value-less flag at the very end of the argument list
		std::vector<std::string> args{
			"ldblocks", "--input-bed", "test.bed", "--n-individuals", "197",
			"--only-groups", "--hash-size", "100", "--add-locus-names"
		};
		std::vector<char *> argv;
		argv.reserve( args.size() );
		for (auto &eachArg : args) {
			argv.push_back( eachArg.data() );
		}
		int argc{static_cast<int>( argv.size() )};
		std::unordered_map<std::string, std::string> cli;
		BayesicSpace::parseCL(argc, argv.data(), cli);
		constexpr size_t correctNflags{5};
		REQUIRE( cli.size() == correctNflags );
		REQUIRE( cli.at("input-bed")     == "test.bed" );
		REQUIRE( cli.at("n-individuals") == "197" );
		REQUIRE( cli.at("hash-size")     == "100" );
		// a value-less flag followed by another flag is recorded as "set"
		REQUIRE( cli.at("only-groups")   == "set" );
		// a value-less flag at the end of the argument list is also recorded as "set"
		REQUIRE( cli.at("add-locus-names") == "set" );

		// nothing but the program name yields an empty map
		std::vector<std::string> noFlagArgs{"ldblocks"};
		std::vector<char *> noFlagArgv{ noFlagArgs.front().data() };
		int noFlagArgc{static_cast<int>( noFlagArgv.size() )};
		std::unordered_map<std::string, std::string> emptyCLI;
		BayesicSpace::parseCL(noFlagArgc, noFlagArgv.data(), emptyCLI);
		REQUIRE( emptyCLI.empty() );
	}

	SECTION("extractCLinfo defaults, overrides, and errors") {
		std::unordered_map<std::string, int>         intVariables;
		std::unordered_map<std::string, float>       floatVariables;
		std::unordered_map<std::string, std::string> stringVariables;

		// an empty parsed map is rejected
		const std::unordered_map<std::string, std::string> emptyCLI;
		REQUIRE_THROWS_WITH( BayesicSpace::extractCLinfo(emptyCLI, intVariables, floatVariables, stringVariables),
				Catch::Matchers::StartsWith("No command line flags specified") );

		// the required integer flag must be present
		const std::unordered_map<std::string, std::string> noNind{ {"input-bed", "test.bed"} };
		REQUIRE_THROWS_WITH( BayesicSpace::extractCLinfo(noNind, intVariables, floatVariables, stringVariables),
				Catch::Matchers::StartsWith("ERROR: n-individuals specification is required and must be an integer") );

		// the required integer flag must be parseable as an integer
		const std::unordered_map<std::string, std::string> badNind{ {"input-bed", "test.bed"}, {"n-individuals", "notAnInt"} };
		REQUIRE_THROWS_WITH( BayesicSpace::extractCLinfo(badNind, intVariables, floatVariables, stringVariables),
				Catch::Matchers::StartsWith("ERROR: n-individuals specification is required and must be an integer") );

		// the required string flag must be present
		const std::unordered_map<std::string, std::string> noBed{ {"n-individuals", "197"} };
		REQUIRE_THROWS_WITH( BayesicSpace::extractCLinfo(noBed, intVariables, floatVariables, stringVariables),
				Catch::Matchers::StartsWith("ERROR: input-bed specification is required") );

		// a minimal valid map fills in all defaults
		const std::unordered_map<std::string, std::string> minimalCLI{ {"input-bed", "test.bed"}, {"n-individuals", "197"} };
		REQUIRE_NOTHROW( BayesicSpace::extractCLinfo(minimalCLI, intVariables, floatVariables, stringVariables) );
		constexpr int correctNind{197};
		REQUIRE( intVariables.at("n-individuals")     == correctNind );
		REQUIRE( intVariables.at("hash-size")         == 0 );
		REQUIRE( intVariables.at("threads")           == -1 );
		REQUIRE( intVariables.at("n-rows-per-band")   == 0 );
		REQUIRE( floatVariables.at("min-similarity")  == 0.0F );
		REQUIRE( stringVariables.at("input-bed")      == "test.bed" );
		REQUIRE( stringVariables.at("log-file")       == "ldblocks.log" );
		REQUIRE( stringVariables.at("out-file")       == "ldblocksOut.tsv" );
		REQUIRE( stringVariables.at("only-groups")    == "unset" );
		REQUIRE( stringVariables.at("add-locus-names") == "unset" );

		// explicit values override the defaults
		const std::unordered_map<std::string, std::string> fullCLI{
			{"input-bed", "data.bed"}, {"n-individuals", "500"}, {"hash-size", "100"},
			{"threads", "4"}, {"n-rows-per-band", "5"}, {"min-similarity", "0.75"},
			{"log-file", "my.log"}, {"out-file", "my.tsv"}, {"only-groups", "set"}, {"add-locus-names", "set"}
		};
		REQUIRE_NOTHROW( BayesicSpace::extractCLinfo(fullCLI, intVariables, floatVariables, stringVariables) );
		constexpr int correctFullNind{500};
		constexpr int correctHashSize{100};
		constexpr int correctThreads{4};
		constexpr int correctNrows{5};
		constexpr float correctMinSim{0.75F}; // exactly representable in binary
		REQUIRE( intVariables.at("n-individuals")     == correctFullNind );
		REQUIRE( intVariables.at("hash-size")         == correctHashSize );
		REQUIRE( intVariables.at("threads")           == correctThreads );
		REQUIRE( intVariables.at("n-rows-per-band")   == correctNrows );
		REQUIRE( floatVariables.at("min-similarity")  == correctMinSim );
		REQUIRE( stringVariables.at("input-bed")      == "data.bed" );
		REQUIRE( stringVariables.at("log-file")       == "my.log" );
		REQUIRE( stringVariables.at("out-file")       == "my.tsv" );
		REQUIRE( stringVariables.at("only-groups")    == "set" );
		REQUIRE( stringVariables.at("add-locus-names") == "set" );

		// an unparseable optional integer silently falls back to its default
		const std::unordered_map<std::string, std::string> badOptional{
			{"input-bed", "test.bed"}, {"n-individuals", "197"}, {"hash-size", "notAnInt"}
		};
		REQUIRE_NOTHROW( BayesicSpace::extractCLinfo(badOptional, intVariables, floatVariables, stringVariables) );
		REQUIRE( intVariables.at("hash-size") == 0 );
	}
}

TEST_CASE("SimilarityMatrix methods work", "[SimilarityMatrix]") {
	constexpr std::array<uint32_t, 7> rowIndexes{4, 5, 5, 6, 6, 8, 8};
	constexpr std::array<uint32_t, 7> colIndexes{3, 1, 2, 2, 4, 1, 7};
	// for Jaccard calculation
	constexpr std::array<uint64_t, 7> nIsect{21, 81, 164, 551, 16, 8, 2};
	constexpr std::array<uint64_t, 7> nUnion{101, 254, 502, 1001, 35, 11, 5};
	// to test insertion
	constexpr std::array<uint32_t, 3> addRowIndexes{8, 3, 5};
	constexpr std::array<uint32_t, 3> addColIndexes{6, 0, 1};
	constexpr std::array<uint64_t, 3> addIsect{89, 44, 6};
	constexpr std::array<uint64_t, 3> addUnion{234, 51, 13};

	constexpr std::array<uint32_t, 9> correctRowIndexes{3, 4, 5, 5, 6, 6, 8, 8, 8};
	constexpr std::array<uint32_t, 9> correctColIndexes{0, 3, 1, 2, 2, 4, 1, 6, 7};
	constexpr std::array<float,    9> correctFloatValues{0.8627, 0.2078, 0.3176, 0.3255, 0.5490, 0.4549, 0.7255, 0.3765, 0.4000};
	constexpr uint64_t previousIdx{7};
	constexpr size_t nThreads{2};

	SECTION("Auxiliary functions") {
		constexpr uint64_t byteSize{8};
		constexpr std::array<uint64_t, 4> vecIdxArray{9, 11, 12, 17};
		constexpr uint64_t correctVecIdx{9};
		BayesicSpace::RowColIdx rowColValues{BayesicSpace::recoverRCindexes( vecIdxArray.at(0) )};
		REQUIRE( rowColValues.iRow == rowIndexes.at(0) );
		REQUIRE( rowColValues.jCol == colIndexes.at(0) );

		// the ascending-order cursor must agree with recoverRCindexes element for element
		BayesicSpace::RowColCursor cursor{ vecIdxArray.at(0) };
		bool cursorMatches{true};
		for (const auto &eachVecIdx : vecIdxArray) {
			const BayesicSpace::RowColIdx expected{BayesicSpace::recoverRCindexes(eachVecIdx)};
			const BayesicSpace::RowColIdx fromCursor{ cursor.advanceTo(eachVecIdx) };
			cursorMatches = cursorMatches && (fromCursor.iRow == expected.iRow) && (fromCursor.jCol == expected.jCol);
		}
		REQUIRE(cursorMatches);

		// a repeated index must not advance the cursor, and every index of a contiguous stretch
		// spanning several rows must resolve correctly
		BayesicSpace::RowColCursor repeatCursor{ vecIdxArray.at(0) };
		const BayesicSpace::RowColIdx firstVisit{ repeatCursor.advanceTo( vecIdxArray.at(0) ) };
		const BayesicSpace::RowColIdx secondVisit{ repeatCursor.advanceTo( vecIdxArray.at(0) ) };
		REQUIRE( firstVisit.iRow == secondVisit.iRow );
		REQUIRE( firstVisit.jCol == secondVisit.jCol );

		constexpr uint64_t nContiguous{200};
		BayesicSpace::RowColCursor sweepCursor{0};
		bool sweepMatches{true};
		for (uint64_t eachVecIdx = 0; eachVecIdx < nContiguous; ++eachVecIdx) {
			const BayesicSpace::RowColIdx expected{BayesicSpace::recoverRCindexes(eachVecIdx)};
			const BayesicSpace::RowColIdx fromCursor{ sweepCursor.advanceTo(eachVecIdx) };
			sweepMatches = sweepMatches && (fromCursor.iRow == expected.iRow) && (fromCursor.jCol == expected.jCol);
		}
		REQUIRE(sweepMatches);
	}
	SECTION("SimilarityMatrix methods") {
		std::array<BayesicSpace::RowColIdx, rowIndexes.size()> idxPairs{};
		std::array<BayesicSpace::JaccardPair, rowIndexes.size()> jaccPairs{};
		size_t vecIdx{0};
		while ( vecIdx < rowIndexes.size() ) {
			idxPairs.at(vecIdx).iRow = rowIndexes.at(vecIdx);
			idxPairs.at(vecIdx).jCol = colIndexes.at(vecIdx);

			jaccPairs.at(vecIdx).nIntersect = nIsect.at(vecIdx);
			jaccPairs.at(vecIdx).nUnion     = nUnion.at(vecIdx);
			++vecIdx;
		}

		std::array<BayesicSpace::RowColIdx, addRowIndexes.size()> addIdxPairs{};
		std::array<BayesicSpace::JaccardPair, addRowIndexes.size()> addJaccPairs{};
		vecIdx = 0;
		while ( vecIdx < addRowIndexes.size() ) {
			addIdxPairs.at(vecIdx).iRow = addRowIndexes.at(vecIdx);
			addIdxPairs.at(vecIdx).jCol = addColIndexes.at(vecIdx);

			addJaccPairs.at(vecIdx).nIntersect = addIsect.at(vecIdx);
			addJaccPairs.at(vecIdx).nUnion     = addUnion.at(vecIdx);
			++vecIdx;
		}

		BayesicSpace::SimilarityMatrix testMatrix;
		constexpr size_t initialSize{0};
		REQUIRE(testMatrix.objectSize() == initialSize);
		REQUIRE( testMatrix.elementSize() == sizeof(uint64_t) );
		vecIdx = 0;
		while ( vecIdx < rowIndexes.size() ) {
			testMatrix.insert( idxPairs.at(vecIdx), jaccPairs.at(vecIdx) );
			++vecIdx;
		}
		REQUIRE( testMatrix.nElements() == rowIndexes.size() );
		// test the last value insertion bypass
		testMatrix.insert( idxPairs.back(), jaccPairs.back() );
		REQUIRE( testMatrix.objectSize() == ( initialSize + ( sizeof(uint64_t) * idxPairs.size() ) ) );
		vecIdx = 0;
		while ( vecIdx < addRowIndexes.size() ) {
			testMatrix.insert( addIdxPairs.at(vecIdx), addJaccPairs.at(vecIdx) );
			++vecIdx;
		}
		REQUIRE( testMatrix.objectSize() == ( initialSize + ( sizeof(uint64_t) * correctFloatValues.size() ) ) );

		// test file save
		const std::string outputFileName("../tests/smallSimilarityMatrix.tsv");
		testMatrix.save(outputFileName, nThreads);
		std::fstream testSMoutfile(outputFileName, std::ios::in);
		std::string line;
		std::array<uint32_t, correctFloatValues.size()> rowsFromFile{};
		std::array<uint32_t, correctFloatValues.size()> colsFromFile{};
		std::array<float,    correctFloatValues.size()> floatsFromFile{};
		size_t arrayIdx{0};
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrayIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrayIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrayIdx) = stof(field);
			++arrayIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correctRowIndexes.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correctColIndexes.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correctFloatValues.cbegin() ) );

		const std::string bimFile("../tests/ind197_397.bim");
		testMatrix.save(outputFileName, nThreads, bimFile);
		arrayIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			floatsFromFile.at(arrayIdx) = stof(field);
			++arrayIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correctFloatValues.cbegin() ) );

		const std::string smallFile("../tests/small.bim");
		REQUIRE_THROWS_WITH(testMatrix.save(outputFileName, nThreads, smallFile),
			Catch::Matchers::StartsWith("ERROR: number of rows exceeds locus name count in") );

		// throwing tests
		BayesicSpace::RowColIdx wrongCombo{};
		BayesicSpace::SimilarityMatrix wrongMatrix;
		wrongCombo.iRow = 1;
		wrongCombo.jCol = 1;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(wrongCombo, jaccPairs.at(0)), 
			Catch::Matchers::StartsWith("ERROR: row and column indexes must be different in")
		);
		constexpr uint32_t tooLarge{379625063};
		wrongCombo.iRow = tooLarge;
		wrongCombo.jCol = tooLarge + 2;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(wrongCombo, jaccPairs.at(0)), 
			Catch::Matchers::StartsWith("ERROR: row or column index exceeds maximal allowable value in")
		);

		wrongCombo.iRow = 0;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(wrongCombo, jaccPairs.at(0)), 
			Catch::Matchers::StartsWith("ERROR: row index must be non-zero in")
		);
		BayesicSpace::JaccardPair wrongPair{};
		wrongPair.nIntersect = 0;
		wrongPair.nUnion     = 0;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(idxPairs.at(0), wrongPair), 
			Catch::Matchers::StartsWith("ERROR: union count cannot be 0 in")
		);
		wrongPair.nIntersect = 2;
		REQUIRE_THROWS_WITH(wrongMatrix.insert(idxPairs.at(0), wrongPair), 
			Catch::Matchers::StartsWith("ERROR: intersection count cannot be larger than the union count in")
		);
	}
	SECTION("SimilarityMatrix merge") {
		constexpr std::array<uint64_t, 12> vecIndexes1{3, 6, 8, 9, 12, 14, 15, 17, 19, 26, 28, 31};
		constexpr std::array<uint64_t, 12> nIsect1{87, 8, 77, 68, 96, 49, 97, 41, 122, 82, 23, 45};
		constexpr std::array<uint64_t, 12> nUnion1{254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254};
		constexpr std::array<uint64_t, 12> vecIndexes2{4, 6, 8, 9, 11, 14, 16, 20, 27, 29, 34, 35};
		// assigning different values to the same row/col pairs for debugging; in actual application they will be the same
		constexpr std::array<uint64_t, 12> nIsect2{217, 160, 228, 176, 167, 171, 228, 206, 174, 214, 201, 161};
		constexpr std::array<uint64_t, 12> nUnion2{254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254};
		constexpr std::array<uint64_t, 3> vecIndexes3{7, 10, 15};
		constexpr std::array<uint64_t, 3> nIsect3{140, 138, 136};
		constexpr std::array<uint64_t, 3> nUnion3{254, 254, 254};
		constexpr uint64_t valueSize{8};

		constexpr size_t nThreads{2};
		BayesicSpace::SimilarityMatrix matrix1;
		size_t arrIdx{0};
		while ( arrIdx < vecIndexes1.size() ) {
			BayesicSpace::RowColIdx tmp{BayesicSpace::recoverRCindexes(vecIndexes1.at(arrIdx))};
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion1.at(arrIdx);
			tmpJP.nIntersect = nIsect1.at(arrIdx);
			matrix1.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix matrix2;
		arrIdx = 0;
		while ( arrIdx < vecIndexes2.size() ) {
			BayesicSpace::RowColIdx tmp{BayesicSpace::recoverRCindexes(vecIndexes2.at(arrIdx))};
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion2.at(arrIdx);
			tmpJP.nIntersect = nIsect2.at(arrIdx);
			matrix2.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix matrix3;
		arrIdx = 0;
		while ( arrIdx < vecIndexes3.size() ) {
			BayesicSpace::RowColIdx tmp{BayesicSpace::recoverRCindexes(vecIndexes3.at(arrIdx))};
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion3.at(arrIdx);
			tmpJP.nIntersect = nIsect3.at(arrIdx);
			matrix3.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix tmp1 = matrix1;
		BayesicSpace::SimilarityMatrix tmp2 = matrix2;

		tmp1.merge(tmp2);
		const std::string outputFileName("../tests/mergeMatrix.tsv");
		tmp1.save(outputFileName, nThreads);
		constexpr std::array<uint32_t, 20> correct12mergeRow{3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 8};
		constexpr std::array<uint32_t, 20> correct12mergeCol{0, 1, 0, 2, 3, 1, 2, 4, 0, 1, 2, 4, 5, 5, 6, 0, 1, 3, 6, 7};
		constexpr std::array<float,    20> correct12mergeValues{
			0.3412, 0.8510, 0.0314, 0.3020, 0.2667, 0.6549, 0.3765, 0.1922, 0.3804, 0.8941,
			0.1608, 0.4784, 0.8078, 0.3216, 0.6824, 0.0902, 0.8392, 0.1765, 0.7882, 0.6314
		};

		std::fstream testSMoutfile(outputFileName, std::ios::in);
		std::string line;
		std::array<uint32_t, correct12mergeValues.size()> rowsFromFile{};
		std::array<uint32_t, correct12mergeValues.size()> colsFromFile{};
		std::array<float,    correct12mergeValues.size()> floatsFromFile{};
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// merge of a matrix with identical tail
		tmp2 = matrix2;
		tmp1.merge(tmp2);
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// merge of an identical matrix
		tmp2 = tmp1;
		tmp1.merge(tmp2);
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// merge of an empty matrix
		tmp1.merge(tmp2);
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		tmp2.merge(tmp1);
		tmp2.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues.cbegin() ) );

		// matrix with smaller indexes as second
		// correct values will change because of the contrived value assignment
		constexpr std::array<float,    20> correct12mergeValues2{
			0.3412, 0.8510, 0.6275, 0.8941, 0.6902, 0.6549, 0.3765, 0.6706, 0.3804, 0.8941,
			0.1608, 0.4784, 0.8078, 0.3216, 0.6824, 0.0902, 0.8392, 0.1765, 0.7882, 0.6314
		};
		tmp1 = matrix1;
		matrix2.merge(tmp1);
		matrix2.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile.cbegin(),   rowsFromFile.cend(),   correct12mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile.cbegin(),   colsFromFile.cend(),   correct12mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile.cbegin(), floatsFromFile.cend(), correct12mergeValues2.cbegin() ) );

		// a matrix completely within another
		constexpr std::array<uint32_t, 14> correct13mergeRow{3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 8};
		constexpr std::array<uint32_t, 14> correct13mergeCol{0, 0, 1, 2, 3, 0, 2, 4, 0, 2, 4, 5, 0, 3};
		constexpr std::array<float,    14> correct13mergeValues{
			0.3412, 0.0314, 0.5490, 0.3020, 0.2667, 0.5412, 0.3765,
			0.1922, 0.3804, 0.1608, 0.4784, 0.3216, 0.0902, 0.1765
		};
		std::array<uint32_t, correct13mergeValues.size()> rowsFromFile13{};
		std::array<uint32_t, correct13mergeValues.size()> colsFromFile13{};
		std::array<float,    correct13mergeValues.size()> floatsFromFile13{};
		tmp1 = matrix1;
		tmp2 = matrix3;
		tmp1.merge(tmp2);
		tmp1.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile13.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile13.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile13.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile13.cbegin(),   rowsFromFile13.cend(),   correct13mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile13.cbegin(),   colsFromFile13.cend(),   correct13mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile13.cbegin(), floatsFromFile13.cend(), correct13mergeValues.cbegin() ) );

		// similarity value preservation is asymmetric
		constexpr std::array<float,    14> correct13mergeValues2{
			0.3412, 0.0314, 0.5490, 0.3020, 0.2667, 0.5412, 0.3765,
			0.1922, 0.5333, 0.1608, 0.4784, 0.3216, 0.0902, 0.1765
		};
		matrix3.merge(matrix1);
		matrix3.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile13.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile13.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile13.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile13.cbegin(),   rowsFromFile13.cend(),   correct13mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile13.cbegin(),   colsFromFile13.cend(),   correct13mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile13.cbegin(), floatsFromFile13.cend(), correct13mergeValues2.cbegin() ) );

		// full append test
		constexpr std::array<uint32_t, 5> vecIndexes4{3, 6, 8, 9, 12};
		constexpr std::array<uint64_t, 5> nIsect4{4, 8, 13, 31, 41};
		constexpr std::array<uint64_t, 5> nUnion4{254, 254, 254, 254, 254};
		constexpr std::array<uint32_t, 7> vecIndexes5{14, 15, 17, 19, 26, 28, 31};
		constexpr std::array<uint64_t, 7> nIsect5{152, 160, 228, 176, 167, 161, 151};
		constexpr std::array<uint64_t, 7> nUnion5{254, 254, 254, 254, 254, 254, 254};

		constexpr std::array<uint32_t, 12> correct45mergeRow{3, 4, 4, 4, 5, 5, 6, 6, 6, 7, 8, 8};
		constexpr std::array<uint32_t, 12> correct45mergeCol{0, 0, 2, 3, 2, 4, 0, 2, 4, 5, 0, 3};
		constexpr std::array<float,    12> correct45mergeValues{
			0.0157, 0.0314, 0.0510, 0.1216, 0.1608, 0.5961,
			0.6275, 0.8941, 0.6902, 0.6549, 0.6314, 0.5922
		};
		std::array<uint32_t, correct45mergeValues.size()> rowsFromFile45{};
		std::array<uint32_t, correct45mergeValues.size()> colsFromFile45{};
		std::array<float,    correct45mergeValues.size()> floatsFromFile45{};

		BayesicSpace::SimilarityMatrix matrix4;
		arrIdx = 0;
		while ( arrIdx < vecIndexes4.size() ) {
			BayesicSpace::RowColIdx tmp{BayesicSpace::recoverRCindexes(vecIndexes4.at(arrIdx))};
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion4.at(arrIdx);
			tmpJP.nIntersect = nIsect4.at(arrIdx);
			matrix4.insert(tmp, tmpJP);
			++arrIdx;
		}
		BayesicSpace::SimilarityMatrix matrix5;
		arrIdx = 0;
		while ( arrIdx < vecIndexes5.size() ) {
			BayesicSpace::RowColIdx tmp{BayesicSpace::recoverRCindexes(vecIndexes5.at(arrIdx))};
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = nUnion5.at(arrIdx);
			tmpJP.nIntersect = nIsect5.at(arrIdx);
			matrix5.insert(tmp, tmpJP);
			++arrIdx;
		}
		matrix4.merge(matrix5);
		matrix4.save(outputFileName, nThreads);
		testSMoutfile.open(outputFileName, std::ios::in);
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFile45.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFile45.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFile45.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFile45.cbegin(),   rowsFromFile45.cend(),   correct45mergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFile45.cbegin(),   colsFromFile45.cend(),   correct45mergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFile45.cbegin(), floatsFromFile45.cend(), correct45mergeValues.cbegin() ) );

		// smaller matrix with an index far past the larger
		constexpr std::array<uint32_t, 3> largeIdxRows{4, 6, 40005};
		constexpr std::array<uint32_t, 3> largeIdxCols{1, 2, 30005};
		constexpr std::array<uint64_t, 3> largeNisect{246, 201, 251};
		constexpr std::array<uint64_t, 3> largeNunion{254, 254, 254};
		BayesicSpace::SimilarityMatrix matrixFar;
		arrIdx = 0;
		while ( arrIdx < largeIdxRows.size() ) {
			BayesicSpace::RowColIdx tmp{};
			tmp.iRow = largeIdxRows.at(arrIdx);
			tmp.jCol = largeIdxCols.at(arrIdx);
			BayesicSpace::JaccardPair tmpJP{};
			tmpJP.nUnion     = largeNunion.at(arrIdx);
			tmpJP.nIntersect = largeNisect.at(arrIdx);
			matrixFar.insert(tmp, tmpJP);
			++arrIdx;
		}
		matrix2.merge(matrixFar);
		matrix2.save(outputFileName, nThreads);
		constexpr std::array<uint32_t, 22> correctFarMergeRow{3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 8, 8, 8, 8, 8, 40005};
		constexpr std::array<uint32_t, 22> correctFarMergeCol{0, 1, 0, 1, 2, 3, 1, 2, 4, 0, 1, 2, 4, 5, 5, 6, 0, 1, 3, 6, 7, 30005};
		constexpr std::array<float,    22> correctFarMergeValues{
			0.3412, 0.8510, 0.6275, 0.9647, 0.8941, 0.6902, 0.6549, 0.3765,
			0.6706, 0.3804, 0.8941, 0.1608, 0.4784, 0.8078, 0.3216, 0.6824,
			0.0902, 0.8392, 0.1765, 0.7882, 0.6314, 0.9843
		};
		testSMoutfile.open(outputFileName, std::ios::in);
		std::array<uint32_t, correctFarMergeValues.size()> rowsFromFileFar{};
		std::array<uint32_t, correctFarMergeValues.size()> colsFromFileFar{};
		std::array<float,    correctFarMergeValues.size()> floatsFromFileFar{};
		arrIdx = 0;
		while ( std::getline(testSMoutfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			rowsFromFileFar.at(arrIdx) = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			colsFromFileFar.at(arrIdx) = stoi(field) - 1;
			lineStream >> field;
			floatsFromFileFar.at(arrIdx) = stof(field);
			++arrIdx;
		}
		testSMoutfile.close();
		std::remove( outputFileName.c_str() ); // NOLINT
		REQUIRE( std::equal( rowsFromFileFar.cbegin(),   rowsFromFileFar.cend(),   correctFarMergeRow.cbegin() ) );
		REQUIRE( std::equal( colsFromFileFar.cbegin(),   colsFromFileFar.cend(),   correctFarMergeCol.cbegin() ) );
		REQUIRE( std::equal( floatsFromFileFar.cbegin(), floatsFromFileFar.cend(), correctFarMergeValues.cbegin() ) );
	}

	SECTION("SimilarityMatrix unordered append and finalize") {
		constexpr size_t nPerMatrix{5};
		constexpr std::array<uint64_t, nPerMatrix> vecIndexesA{3, 8, 14, 19, 28};
		constexpr std::array<uint64_t, nPerMatrix> nIsectA{87, 77, 49, 122, 23};
		constexpr std::array<uint64_t, nPerMatrix> vecIndexesB{6, 9, 12, 17, 31};   // disjoint from A but interleaving its indexes
		constexpr std::array<uint64_t, nPerMatrix> nIsectB{160, 176, 167, 206, 161};
		constexpr uint64_t nUnionVal{254};
		constexpr size_t nThreads{2};
		const std::string appendFileName("../tests/appendMatrix.tsv");
		const std::string refFileName("../tests/appendRefMatrix.tsv");

		auto buildMatrix = [](const std::array<uint64_t, nPerMatrix> &indexes, const std::array<uint64_t, nPerMatrix> &isect, const uint64_t &nUnion) {
			BayesicSpace::SimilarityMatrix mat;
			for (size_t idx = 0; idx < indexes.size(); ++idx) {
				BayesicSpace::RowColIdx rowCol{BayesicSpace::recoverRCindexes(indexes.at(idx))};
				BayesicSpace::JaccardPair jPair{};
				jPair.nUnion     = nUnion;
				jPair.nIntersect = isect.at(idx);
				mat.insert(rowCol, jPair);
			}
			return mat;
		};
		auto readMatrixFile = [](const std::string &fileName) {
			std::vector<std::string> fileRows;
			std::fstream inFile(fileName, std::ios::in);
			std::string fileLine;
			while ( std::getline(inFile, fileLine) ) {
				fileRows.push_back(fileLine);
			}
			inFile.close();
			return fileRows;
		};

		const BayesicSpace::SimilarityMatrix matrixA{buildMatrix(vecIndexesA, nIsectA, nUnionVal)};
		const BayesicSpace::SimilarityMatrix matrixB{buildMatrix(vecIndexesB, nIsectB, nUnionVal)};

		// unordered append then finalize matches merge of the same disjoint matrices
		BayesicSpace::SimilarityMatrix appended{matrixA};
		BayesicSpace::SimilarityMatrix bToAppend{matrixB};
		appended.append(bToAppend);
		REQUIRE( bToAppend.nElements() == 0 );                                     // the source is cleared
		REQUIRE( appended.nElements() == vecIndexesA.size() + vecIndexesB.size() ); // no ordering/de-duplication yet
		appended.sortAndDeduplicate();
		appended.save(appendFileName, nThreads);

		BayesicSpace::SimilarityMatrix reference{matrixA};
		BayesicSpace::SimilarityMatrix bToMerge{matrixB};
		reference.merge(bToMerge);
		reference.save(refFileName, nThreads);

		const std::vector<std::string> appendedLines{readMatrixFile(appendFileName)};
		const std::vector<std::string> referenceLines{readMatrixFile(refFileName)};
		std::remove( appendFileName.c_str() ); // NOLINT
		std::remove( refFileName.c_str() );    // NOLINT
		REQUIRE( appendedLines.size() == vecIndexesA.size() + vecIndexesB.size() ); // sorted, no spurious drops
		REQUIRE( appendedLines == referenceLines );

		// duplicate indexes (identical values) collapse to a single entry on finalize
		BayesicSpace::SimilarityMatrix withDuplicates{matrixA};
		BayesicSpace::SimilarityMatrix duplicateA{matrixA};
		withDuplicates.append(duplicateA);
		REQUIRE( withDuplicates.nElements() == 2 * vecIndexesA.size() );
		withDuplicates.sortAndDeduplicate();
		REQUIRE( withDuplicates.nElements() == vecIndexesA.size() );
		withDuplicates.save(appendFileName, nThreads);
		const std::vector<std::string> dedupLines{readMatrixFile(appendFileName)};
		std::remove( appendFileName.c_str() ); // NOLINT
		matrixA.save(refFileName, nThreads);
		const std::vector<std::string> aloneLines{readMatrixFile(refFileName)};
		std::remove( refFileName.c_str() ); // NOLINT
		REQUIRE( dedupLines == aloneLines );
	}
	SECTION("SimilarityMatrix k-way merge of pre-sorted runs") {
		constexpr uint64_t nUnionVal{254};
		constexpr size_t nThreads{2};
		const std::string mergeFileName("../tests/kwayMergeMatrix.tsv");
		const std::string refFileName("../tests/kwayRefMatrix.tsv");

		auto buildRun = [](const std::vector<uint64_t> &indexes, const std::vector<uint64_t> &isect, const uint64_t &nUnion) {
			BayesicSpace::SimilarityMatrix mat;
			for (size_t idx = 0; idx < indexes.size(); ++idx) {
				BayesicSpace::RowColIdx rowCol{BayesicSpace::recoverRCindexes(indexes.at(idx))};
				BayesicSpace::JaccardPair jPair{};
				jPair.nUnion     = nUnion;
				jPair.nIntersect = isect.at(idx);
				mat.insert(rowCol, jPair);
			}
			return mat;
		};
		auto readMatrixFile = [](const std::string &fileName) {
			std::vector<std::string> fileRows;
			std::fstream inFile(fileName, std::ios::in);
			std::string fileLine;
			while ( std::getline(inFile, fileLine) ) {
				fileRows.push_back(fileLine);
			}
			inFile.close();
			return fileRows;
		};

		// three individually sorted runs; index 8 is shared by runA and runC, index 12 by runB and runC,
		// each duplicate carrying the same value it has in the other run (as the production paths do)
		const std::vector<uint64_t> idxA{3, 8, 14, 19, 28};
		const std::vector<uint64_t> isA {87, 77, 49, 122, 23};
		const std::vector<uint64_t> idxB{6, 9, 12, 17, 31};
		const std::vector<uint64_t> isB {160, 176, 167, 206, 161};
		const std::vector<uint64_t> idxC{8, 12, 20};
		const std::vector<uint64_t> isC {77, 167, 51};                             // 8 and 12 match runA/runB values

		auto makeRuns = [&]() {
			std::vector<BayesicSpace::SimilarityMatrix> runs;
			runs.emplace_back(buildRun(idxA, isA, nUnionVal));
			runs.emplace_back(buildRun(idxB, isB, nUnionVal));
			runs.emplace_back(buildRun(idxC, isC, nUnionVal));
			return runs;
		};

		std::vector<BayesicSpace::SimilarityMatrix> runs{makeRuns()};
		BayesicSpace::SimilarityMatrix merged{BayesicSpace::SimilarityMatrix::mergeSortedRuns(runs)};
		for (const auto &eachRun : runs) {
			REQUIRE( eachRun.nElements() == 0 );                                   // every run is cleared
		}
		REQUIRE( merged.nElements() == 11 );                                       // (5 + 5 + 3) - 2 shared
		merged.save(mergeFileName, nThreads);

		// reference: append all runs then finalize (the established path mergeSortedRuns replaces)
		std::vector<BayesicSpace::SimilarityMatrix> refRuns{makeRuns()};
		BayesicSpace::SimilarityMatrix reference;
		for (auto &eachRun : refRuns) {
			reference.append(eachRun);
		}
		reference.sortAndDeduplicate();
		reference.save(refFileName, nThreads);

		const std::vector<std::string> mergedLines{readMatrixFile(mergeFileName)};
		const std::vector<std::string> referenceLines{readMatrixFile(refFileName)};
		std::remove( mergeFileName.c_str() ); // NOLINT
		std::remove( refFileName.c_str() );   // NOLINT
		REQUIRE( mergedLines.size() == 11 );
		REQUIRE( mergedLines == referenceLines );

		// Parallel key-range path: a large multi-run input (past the serial/parallel threshold) with
		// cross-run duplicates, checked against the append + sortAndDeduplicate reference.
		constexpr size_t nBigRuns{5};
		constexpr size_t perRun{80000};
		constexpr uint64_t hotBase{1000000};                                       // shared "hot" indexes, above the strided range
		constexpr size_t overlapN{500};
		constexpr size_t forceParallelThreads{8};                                  // 5 * 80500 elements over 65536 -> parallel path
		constexpr uint64_t valueModulus{200};                                      // <= nUnionVal, keeps quantized values in range
		auto valueFor = [](const uint64_t vecIdx) { return vecIdx % valueModulus; };   // identical for a shared index
		auto buildBigRun = [&](const size_t runIndex) {
			std::vector<uint64_t> indexes;
			indexes.reserve(perRun + overlapN);
			for (size_t elem = 0; elem < perRun; ++elem) {
				indexes.push_back( static_cast<uint64_t>( (elem * nBigRuns) + runIndex ) );   // strided, disjoint across runs
			}
			for (size_t hot = 0; hot < overlapN; ++hot) {
				indexes.push_back( hotBase + hot );                                // shared by every run -> duplicates
			}
			std::vector<uint64_t> isect( indexes.size() );
			std::transform( indexes.cbegin(), indexes.cend(), isect.begin(), valueFor );
			return buildRun(indexes, isect, nUnionVal);                            // ascending indexes -> fast push_back path
		};
		auto makeBigRuns = [&]() {
			std::vector<BayesicSpace::SimilarityMatrix> bigRuns;
			bigRuns.reserve(nBigRuns);
			for (size_t run = 0; run < nBigRuns; ++run) {
				bigRuns.emplace_back( buildBigRun(run) );
			}
			return bigRuns;
		};

		std::vector<BayesicSpace::SimilarityMatrix> bigRuns{makeBigRuns()};
		BayesicSpace::SimilarityMatrix bigMerged{BayesicSpace::SimilarityMatrix::mergeSortedRuns(bigRuns, forceParallelThreads)};
		for (const auto &eachRun : bigRuns) {
			REQUIRE( eachRun.nElements() == 0 );                                   // every run is freed
		}
		REQUIRE( bigMerged.nElements() == (nBigRuns * perRun) + overlapN );        // strided all distinct, hot collapsed to one each

		const std::string bigMergeFile("../tests/kwayBigMerge.tsv");
		const std::string bigRefFile("../tests/kwayBigRef.tsv");
		bigMerged.save(bigMergeFile, nThreads);
		std::vector<BayesicSpace::SimilarityMatrix> bigRefRuns{makeBigRuns()};
		BayesicSpace::SimilarityMatrix bigReference;
		for (auto &eachRun : bigRefRuns) {
			bigReference.append(eachRun);
		}
		bigReference.sortAndDeduplicate();
		bigReference.save(bigRefFile, nThreads);
		const std::vector<std::string> bigMergedLines{readMatrixFile(bigMergeFile)};
		const std::vector<std::string> bigReferenceLines{readMatrixFile(bigRefFile)};
		std::remove( bigMergeFile.c_str() ); // NOLINT
		std::remove( bigRefFile.c_str() );   // NOLINT
		REQUIRE( bigMergedLines.size() == (nBigRuns * perRun) + overlapN );
		REQUIRE( bigMergedLines == bigReferenceLines );                            // parallel merge matches the serial reference

		std::vector<BayesicSpace::SimilarityMatrix> noRuns;
		REQUIRE( BayesicSpace::SimilarityMatrix::mergeSortedRuns(noRuns).nElements() == 0 ); // no runs
		std::vector<BayesicSpace::SimilarityMatrix> emptyRuns(3);
		REQUIRE( BayesicSpace::SimilarityMatrix::mergeSortedRuns(emptyRuns).nElements() == 0 ); // all empty
	}
	SECTION("SimilarityMatrixSink streams within a memory budget") {
		constexpr size_t nPerMatrix{5};
		constexpr std::array<uint64_t, nPerMatrix> vecIndexesA{3, 8, 14, 19, 28};
		constexpr std::array<uint64_t, nPerMatrix> nIsectA{87, 77, 49, 122, 23};
		constexpr std::array<uint64_t, nPerMatrix> vecIndexesB{6, 9, 12, 17, 31};   // disjoint from A but interleaving its indexes
		constexpr std::array<uint64_t, nPerMatrix> nIsectB{160, 176, 167, 206, 161};
		constexpr uint64_t nUnionVal{254};
		constexpr size_t nThreads{2};
		const std::string sinkFileName("../tests/sinkMatrix.tsv");
		const std::string refFileName("../tests/sinkRefMatrix.tsv");

		auto buildMatrix = [](const std::array<uint64_t, nPerMatrix> &indexes, const std::array<uint64_t, nPerMatrix> &isect, const uint64_t &nUnion) {
			BayesicSpace::SimilarityMatrix mat;
			for (size_t idx = 0; idx < indexes.size(); ++idx) {
				BayesicSpace::RowColIdx rowCol{BayesicSpace::recoverRCindexes(indexes.at(idx))};
				BayesicSpace::JaccardPair jPair{};
				jPair.nUnion     = nUnion;
				jPair.nIntersect = isect.at(idx);
				mat.insert(rowCol, jPair);
			}
			return mat;
		};
		auto readMatrixFile = [](const std::string &fileName) {
			std::vector<std::string> fileRows;
			std::fstream inFile(fileName, std::ios::in);
			std::string fileLine;
			while ( std::getline(inFile, fileLine) ) {
				fileRows.push_back(fileLine);
			}
			inFile.close();
			return fileRows;
		};

		const BayesicSpace::SimilarityMatrix matrixA{buildMatrix(vecIndexesA, nIsectA, nUnionVal)};
		const BayesicSpace::SimilarityMatrix matrixB{buildMatrix(vecIndexesB, nIsectB, nUnionVal)};

		// the merge of A and B is the reference full result (globally sorted)
		BayesicSpace::SimilarityMatrix reference{matrixA};
		BayesicSpace::SimilarityMatrix bToMerge{matrixB};
		reference.merge(bToMerge);
		reference.save(refFileName, nThreads);
		const std::vector<std::string> referenceLines{readMatrixFile(refFileName)};
		std::remove( refFileName.c_str() ); // NOLINT

		// a budget large enough to hold both blocks flushes only at finalize: file matches merge exactly
		std::remove( sinkFileName.c_str() ); // NOLINT (save appends, so start clean)
		{
			BayesicSpace::SimilarityMatrixSink sink(BayesicSpace::InOutFileNames{std::string(), sinkFileName}, BayesicSpace::WorkloadLimits{nThreads, 2 * nPerMatrix});
			BayesicSpace::SimilarityMatrix blockA{matrixA};
			BayesicSpace::SimilarityMatrix blockB{matrixB};
			sink.add(blockA);
			sink.add(blockB);
			REQUIRE( blockA.nElements() == 0 );                       // added blocks are consumed
			REQUIRE( blockB.nElements() == 0 );
			REQUIRE( sink.bufferedElements() == 2 * nPerMatrix );     // nothing flushed yet
			sink.finalize();
			REQUIRE( sink.bufferedElements() == 0 );                  // buffer emptied on finalize
		}
		const std::vector<std::string> bufferedLines{readMatrixFile(sinkFileName)};
		std::remove( sinkFileName.c_str() ); // NOLINT
		REQUIRE( bufferedLines == referenceLines );

		// a budget that holds one block but not two forces a flush between the adds; the streamed
		// file holds every pair (each flush is internally sorted, so the file is a set-equal permutation)
		std::remove( sinkFileName.c_str() ); // NOLINT
		{
			BayesicSpace::SimilarityMatrixSink sink(BayesicSpace::InOutFileNames{std::string(), sinkFileName}, BayesicSpace::WorkloadLimits{nThreads, nPerMatrix});
			BayesicSpace::SimilarityMatrix blockA{matrixA};
			BayesicSpace::SimilarityMatrix blockB{matrixB};
			sink.add(blockA);
			sink.add(blockB);                                         // triggers a flush of A before appending B
			REQUIRE( sink.bufferedElements() == nPerMatrix );        // only B remains buffered
			sink.finalize();
		}
		std::vector<std::string> streamedLines{readMatrixFile(sinkFileName)};
		std::remove( sinkFileName.c_str() ); // NOLINT
		std::vector<std::string> sortedReference{referenceLines};
		std::sort( streamedLines.begin(), streamedLines.end() );
		std::sort( sortedReference.begin(), sortedReference.end() );
		REQUIRE( streamedLines.size() == 2 * nPerMatrix );           // no pair dropped across the flush
		REQUIRE( streamedLines == sortedReference );

		// duplicate indexes within a single buffer collapse on flush
		std::remove( sinkFileName.c_str() ); // NOLINT
		{
			BayesicSpace::SimilarityMatrixSink sink(BayesicSpace::InOutFileNames{std::string(), sinkFileName}, BayesicSpace::WorkloadLimits{nThreads, 2 * nPerMatrix});
			BayesicSpace::SimilarityMatrix blockA{matrixA};
			BayesicSpace::SimilarityMatrix duplicateA{matrixA};
			sink.add(blockA);
			sink.add(duplicateA);
			sink.finalize();
		}
		const std::vector<std::string> dedupLines{readMatrixFile(sinkFileName)};
		std::remove( sinkFileName.c_str() ); // NOLINT
		matrixA.save(refFileName, nThreads);
		const std::vector<std::string> aloneLines{readMatrixFile(refFileName)};
		std::remove( refFileName.c_str() ); // NOLINT
		REQUIRE( dedupLines == aloneLines );
	}
	SECTION("Buffered save chunks within a byte budget") {
		constexpr size_t nPerMatrix{5};
		constexpr std::array<uint64_t, nPerMatrix> vecIndexesA{3, 8, 14, 19, 28};
		constexpr std::array<uint64_t, nPerMatrix> nIsectA{87, 77, 49, 122, 23};
		constexpr std::array<uint64_t, nPerMatrix> vecIndexesB{6, 9, 12, 17, 31};
		constexpr std::array<uint64_t, nPerMatrix> nIsectB{160, 176, 167, 206, 161};
		constexpr uint64_t nUnionVal{254};
		constexpr size_t nThreads{2};
		const std::string bufFileName("../tests/bufSaveMatrix.tsv");
		const std::string refFileName("../tests/bufSaveRefMatrix.tsv");

		auto buildMatrix = [](const std::array<uint64_t, nPerMatrix> &indexes, const std::array<uint64_t, nPerMatrix> &isect, const uint64_t &nUnion) {
			BayesicSpace::SimilarityMatrix mat;
			for (size_t idx = 0; idx < indexes.size(); ++idx) {
				BayesicSpace::RowColIdx rowCol{BayesicSpace::recoverRCindexes(indexes.at(idx))};
				BayesicSpace::JaccardPair jPair{};
				jPair.nUnion     = nUnion;
				jPair.nIntersect = isect.at(idx);
				mat.insert(rowCol, jPair);
			}
			return mat;
		};
		auto readMatrixFile = [](const std::string &fileName) {
			std::vector<std::string> fileRows;
			std::fstream inFile(fileName, std::ios::in);
			std::string fileLine;
			while ( std::getline(inFile, fileLine) ) {
				fileRows.push_back(fileLine);
			}
			inFile.close();
			return fileRows;
		};

		BayesicSpace::SimilarityMatrix matrix{buildMatrix(vecIndexesA, nIsectA, nUnionVal)};
		BayesicSpace::SimilarityMatrix second{buildMatrix(vecIndexesB, nIsectB, nUnionVal)};
		matrix.merge(second);   // ten elements, globally sorted

		// reference: the RAM-sized public overload writes the whole matrix in one pass
		std::remove( refFileName.c_str() ); // NOLINT
		matrix.save(refFileName, nThreads);
		const std::vector<std::string> referenceLines{readMatrixFile(refFileName)};
		std::remove( refFileName.c_str() ); // NOLINT

		// a tiny reserved string budget forces save() to write in many small chunks, reusing the
		// buffers; the output must be byte-for-byte identical to the single-pass reference
		std::remove( bufFileName.c_str() ); // NOLINT
		matrix.reserve(1);                  // string budget rounds down to a few bytes, ~one line per buffer
		matrix.save(bufFileName, nThreads);
		const std::vector<std::string> bufferedLines{readMatrixFile(bufFileName)};
		std::remove( bufFileName.c_str() ); // NOLINT

		REQUIRE( bufferedLines.size() == vecIndexesA.size() + vecIndexesB.size() ); // every element written
		REQUIRE( bufferedLines == referenceLines );                                // chunking preserves content and order
	}

	SECTION("Saved indexes are recovered across row boundaries") {
		// save() tracks the row incrementally rather than recovering it per element, re-seeding once
		// per slice; these indexes are the first and last column of each listed row, so slices begin
		// at row starts, at row ends and after jumps of many rows.
		constexpr std::array<uint64_t, 6> testRows{1, 2, 3, 100, 101, 4096};
		constexpr uint64_t nUnionVal{254};
		constexpr uint64_t isectStep{20};
		constexpr uint64_t isectBase{7};
		constexpr size_t nThreads{4};
		const std::string outFileName("../tests/rowBoundarySave.tsv");
		const std::string bimFile("../tests/ind197_397.bim");

		std::vector<uint64_t> vecIndexes;
		for (const auto &eachRow : testRows) {
			const uint64_t rowStart{eachRow * (eachRow - 1) / 2};
			vecIndexes.push_back(rowStart);                        // first column of the row
			if (eachRow > 1) {                                     // row 1 has a single element
				vecIndexes.push_back(rowStart + eachRow - 1);      // last column of the row
			}
		}
		BayesicSpace::SimilarityMatrix matrix;
		for (size_t idx = 0; idx < vecIndexes.size(); ++idx) {
			BayesicSpace::JaccardPair jPair{};
			jPair.nUnion     = nUnionVal;
			jPair.nIntersect = (static_cast<uint64_t>(idx) * isectStep) + isectBase;
			matrix.insert(BayesicSpace::recoverRCindexes( vecIndexes.at(idx) ), jPair);
		}
		REQUIRE( matrix.nElements() == vecIndexes.size() );

		auto readIndexPairs = [](const std::string &fileName) {
			std::vector<BayesicSpace::RowColIdx> filePairs;
			std::fstream inFile(fileName, std::ios::in);
			std::string fileLine;
			while ( std::getline(inFile, fileLine) ) {
				std::stringstream lineStream;
				lineStream.str(fileLine);
				std::string field;
				BayesicSpace::RowColIdx parsedPair{};
				lineStream >> field;
				parsedPair.iRow = static_cast<uint32_t>( std::stoul(field) ) - 1;  // the saved indexes are base-1
				lineStream >> field;
				parsedPair.jCol = static_cast<uint32_t>( std::stoul(field) ) - 1;
				filePairs.push_back(parsedPair);
			}
			inFile.close();
			return filePairs;
		};

		std::remove( outFileName.c_str() ); // NOLINT
		matrix.save(outFileName, nThreads);
		const std::vector<BayesicSpace::RowColIdx> savedPairs{ readIndexPairs(outFileName) };

		// the line layout save() assumes when sizing its buffers: two un-padded base-1 index
		// fields and a fixed-width value field, tab-separated
		constexpr size_t valueFieldWidth{6};
		std::fstream formatFile(outFileName, std::ios::in);
		std::string formatLine;
		bool allLinesWellFormed{true};
		size_t nFormatLines{0};
		while ( std::getline(formatFile, formatLine) ) {
			const size_t firstTab{ formatLine.find('\t') };
			const size_t secondTab{ formatLine.find('\t', firstTab + 1) };
			allLinesWellFormed = allLinesWellFormed
					&& (firstTab != std::string::npos) && (secondTab != std::string::npos)
					&& (formatLine.find('\t', secondTab + 1) == std::string::npos)
					&& (formatLine.size() - secondTab - 1 == valueFieldWidth)
					&& (formatLine.front() != '0') && (formatLine.at(firstTab + 1) != '0');
			++nFormatLines;
		}
		formatFile.close();
		std::remove( outFileName.c_str() ); // NOLINT
		REQUIRE( nFormatLines == vecIndexes.size() );
		REQUIRE(allLinesWellFormed);

		REQUIRE( savedPairs.size() == vecIndexes.size() );
		bool allPairsMatch{true};
		for (size_t idx = 0; idx < vecIndexes.size(); ++idx) {
			const BayesicSpace::RowColIdx expected{BayesicSpace::recoverRCindexes( vecIndexes.at(idx) )};
			allPairsMatch = allPairsMatch && (savedPairs.at(idx).iRow == expected.iRow)
											&& (savedPairs.at(idx).jCol == expected.jCol);
		}
		REQUIRE(allPairsMatch);                              // matches the per-element square-root recovery

		// a one-line-per-buffer budget re-seeds the row on every element
		std::remove( outFileName.c_str() ); // NOLINT
		matrix.reserve(1);
		matrix.save(outFileName, nThreads);
		const std::vector<BayesicSpace::RowColIdx> chunkedPairs{ readIndexPairs(outFileName) };
		std::remove( outFileName.c_str() ); // NOLINT
		REQUIRE( chunkedPairs.size() == savedPairs.size() );
		bool chunkedMatch{true};
		for (size_t idx = 0; idx < savedPairs.size(); ++idx) {
			chunkedMatch = chunkedMatch && (chunkedPairs.at(idx).iRow == savedPairs.at(idx).iRow)
										&& (chunkedPairs.at(idx).jCol == savedPairs.at(idx).jCol);
		}
		REQUIRE(chunkedMatch);

		// the same rows resolved to locus names, dropping the rows beyond the .bim file
		const std::vector<std::string> locusNames{ BayesicSpace::getLocusNames(bimFile) };
		BayesicSpace::SimilarityMatrix namedMatrix;
		std::vector<uint64_t> namedIndexes;
		for (const auto &eachIndex : vecIndexes) {
			if ( BayesicSpace::recoverRCindexes(eachIndex).iRow < locusNames.size() ) {
				namedIndexes.push_back(eachIndex);
				BayesicSpace::JaccardPair jPair{};
				jPair.nUnion     = nUnionVal;
				jPair.nIntersect = isectBase;
				namedMatrix.insert(BayesicSpace::recoverRCindexes(eachIndex), jPair);
			}
		}
		std::remove( outFileName.c_str() ); // NOLINT
		namedMatrix.save(outFileName, nThreads, bimFile);
		std::vector<std::string> nameFields;
		std::fstream namedFile(outFileName, std::ios::in);
		std::string namedLine;
		while ( std::getline(namedFile, namedLine) ) {
			std::stringstream lineStream;
			lineStream.str(namedLine);
			std::string field;
			lineStream >> field;
			nameFields.push_back(field);
			lineStream >> field;
			nameFields.push_back(field);
		}
		namedFile.close();
		std::remove( outFileName.c_str() ); // NOLINT
		REQUIRE( nameFields.size() == 2 * namedIndexes.size() );
		bool allNamesMatch{true};
		for (size_t idx = 0; idx < namedIndexes.size(); ++idx) {
			const BayesicSpace::RowColIdx expected{BayesicSpace::recoverRCindexes( namedIndexes.at(idx) )};
			allNamesMatch = allNamesMatch && ( nameFields.at(2 * idx) == locusNames.at(expected.iRow) )
										&& ( nameFields.at( (2 * idx) + 1 ) == locusNames.at(expected.jCol) );
		}
		REQUIRE(allNamesMatch);
	}
}

TEST_CASE("GenoTableBin methods work", "[gtBin]") {
	const std::string logFileName("../tests/binTest.log");
	const std::string inputBedName("../tests/ind197_397.bed");
	constexpr uint32_t nIndividuals{197};
	constexpr size_t nThreads{4};
	SECTION("Failed GenoTableBin constructors") {
		constexpr size_t smallNind{1};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(inputBedName, smallNind, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		const std::string absentFileName("../tests/noSuchFile.bed");
		const std::string noLociFile("../tests/threeByte.bed");
		const std::string wrongMagicBytes("../tests/wrongMB.bed");
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(absentFileName, nIndividuals, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: failed to open file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(noLociFile, nIndividuals, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: no genotype records in file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(wrongMagicBytes, nIndividuals, logFileName, nThreads),
				Catch::Matchers::StartsWith("ERROR: first magic byte in input .bed file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(inputBedName, nIndividuals, logFileName, nThreads, BayesicSpace::MemoryParameters{1}),
				Catch::Matchers::StartsWith("ERROR: the genotype table does not fit within the memory budget") );
		const std::vector<int> smallMACvec(13, 0);
		const std::vector<int> emptyMACvec{};
		constexpr size_t undivNind{5};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(emptyMACvec, nIndividuals, logFileName),
				Catch::Matchers::StartsWith("ERROR: empty vector of minor allele counts") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(smallMACvec, smallNind, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableBin(smallMACvec, undivNind, logFileName),
				Catch::Matchers::StartsWith("ERROR: length of allele count vector") );
	}
	SECTION("GenoTableBin constructors and methods with correct data") {
		constexpr size_t nChunks{3};
		BayesicSpace::GenoTableBin bedGTB(inputBedName, nIndividuals, logFileName, nThreads);
		const std::string ldFileName("../tests/tmpLDfile.tsv");
		std::fstream tmpLDfile;
		std::string line;
		BayesicSpace::InOutFileNames outAndBim{};
		outAndBim.outputFileName = ldFileName;
		bedGTB.allJaccardLD(outAndBim, nChunks);

		tmpLDfile.open(ldFileName, std::ios::in);
		std::vector<float> jaccValues;
		std::getline(tmpLDfile, line); // header
		while ( std::getline(tmpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.push_back( stof(field) );
		}
		tmpLDfile.close();
		std::remove( ldFileName.c_str() ); // NOLINT
		constexpr float upperCutOff{0.9F};
		constexpr uint32_t correctNlargeLD{2400};
		constexpr float lowerCutOff{0.1F};
		constexpr uint32_t correctNsmallLD{48000};
		uint32_t nLargeLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[upperCutOff](float eachValue) {return eachValue >= upperCutOff;}
		);
		REQUIRE(nLargeLD >= correctNlargeLD); // cannot test equality b/c of randomness
		uint32_t nSmallLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[lowerCutOff](float eachValue) {return eachValue <= lowerCutOff;}
		);
		REQUIRE(nSmallLD >= correctNsmallLD); // cannot test equality b/c of randomness
		REQUIRE( nSmallLD + nLargeLD < jaccValues.size() );

		const std::string alleleCountsFile("../tests/alleleCounts.txt");
		std::fstream inAlleleCounts;
		std::string eachLine;
		inAlleleCounts.open(alleleCountsFile, std::ios::in);
		std::vector<int> macVector;
		while ( std::getline(inAlleleCounts, eachLine) ) {
			macVector.push_back( std::stoi(eachLine) );
		}
		inAlleleCounts.close();
		BayesicSpace::GenoTableBin macGTB(macVector, nIndividuals, logFileName, nThreads);
		outAndBim.outputFileName = ldFileName;
		macGTB.allJaccardLD(outAndBim, nChunks);

		tmpLDfile.open(ldFileName, std::ios::in);
		jaccValues.clear();
		std::getline(tmpLDfile, line); // header
		while ( std::getline(tmpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.push_back( stof(field) );
		}
		tmpLDfile.close();
		std::remove( ldFileName.c_str() ); // NOLINT
		nLargeLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[upperCutOff](float eachValue) {return eachValue >= upperCutOff;}
		);
		REQUIRE(nLargeLD >= correctNlargeLD); // cannot test equality b/c of randomness
		nSmallLD = std::count_if(
			jaccValues.cbegin(),
			jaccValues.cend(), 
			[lowerCutOff](float eachValue) {return eachValue <= lowerCutOff;}
		);
		REQUIRE(nSmallLD >= correctNsmallLD); // cannot test equality b/c of randomness
		REQUIRE( nSmallLD + nLargeLD < jaccValues.size() );
	}
	SECTION("saveGenoBinary round-trip") {
		constexpr uint32_t nIndividualsSmall{17};
		constexpr size_t nLociSmall{4};
		constexpr size_t binLocusSize{(nIndividualsSmall / 8) + static_cast<size_t>( (nIndividualsSmall % 8) > 0 )};
		// deterministic minor-allele counts: no heterozygotes (value 1), so binarization is reproducible
		const std::vector<int> macVector{
			0, 2, 0, 2, 2, 0, -9,  0, 0, 2, 0, 2, 0, -9, 2, 0, 2,  // locus 0
			2, 2, 2, 0, 0, 2,  2, -9, 0, 2, 0, 0, 2,  2, 0, 2, 0,  // locus 1
			0, 0, 0, 2, 2, 0,  0,  2, -9, 2, 2, 0, 2,  0, 0, 2, 2, // locus 2
			2, 0, 2, 2, -9, 2, 0,  0, 2, 2, 2, 0, 2,  2, 0, 2, 2   // locus 3
		};
		REQUIRE( macVector.size() == nLociSmall * nIndividualsSmall );

		BayesicSpace::GenoTableBin gtb(macVector, nIndividualsSmall, logFileName, nThreads);
		const std::string binFileName("../tests/tmpGenoBinary.bin");
		gtb.saveGenoBinary(binFileName);

		// read the saved bytes back
		std::fstream binIn(binFileName, std::ios::in | std::ios::binary);
		const std::vector<char> savedBytes( (std::istreambuf_iterator<char>(binIn)), std::istreambuf_iterator<char>() );
		binIn.close();
		std::remove( binFileName.c_str() ); // NOLINT

		// independently reconstruct the expected binary buffer locus by locus
		std::vector<uint8_t> expected(nLociSmall * binLocusSize, 0);
		for (size_t iLocus = 0; iLocus < nLociSmall; ++iLocus) {
			const std::vector<int> macLocus(
				std::next( macVector.cbegin(), static_cast<std::ptrdiff_t>(iLocus * nIndividualsSmall) ),
				std::next( macVector.cbegin(), static_cast<std::ptrdiff_t>( (iLocus + 1) * nIndividualsSmall ) )
			);
			const BayesicSpace::LocationWithLength binWindow{iLocus, binLocusSize};
			BayesicSpace::binarizeMacLocus(macLocus, binWindow, expected);
		}

		REQUIRE( savedBytes.size() == expected.size() );
		REQUIRE( std::equal(
				savedBytes.cbegin(),
				savedBytes.cend(),
				expected.cbegin(),
				[](char savedByte, uint8_t expectedByte) {
					return static_cast<uint8_t>(savedByte) == expectedByte;
				}
			)
		);
	}
	SECTION("The log is flushed to file on destruction") {
		const std::string slfLogName("../tests/saveLogTestBin.log");
		std::remove( slfLogName.c_str() ); // start from a clean slate // NOLINT
		{
			BayesicSpace::GenoTableBin logBin(inputBedName, nIndividuals, slfLogName, nThreads);
		} // destructor flushes the accumulated log here
		std::fstream logIn(slfLogName, std::ios::in);
		REQUIRE( logIn.good() );
		const std::string logContents( (std::istreambuf_iterator<char>(logIn)), std::istreambuf_iterator<char>() );
		logIn.close();
		std::remove( slfLogName.c_str() ); // NOLINT
		// the constructor accumulates log messages, so the saved file must be non-empty
		REQUIRE( !logContents.empty() );
	}
	SECTION("An empty log file name disables log writing") {
		const std::string slfLogName("../tests/saveLogTestBinEmpty.log");
		std::remove( slfLogName.c_str() ); // start from a clean slate // NOLINT
		{
			BayesicSpace::GenoTableBin logBin(inputBedName, nIndividuals, std::string(), nThreads);
		} // destructor must not create a file
		std::fstream logIn(slfLogName, std::ios::in);
		REQUIRE( !logIn.good() );
	}
	SECTION("Chunked .bed reading matches single-chunk reading") {
		// Regression test for the chunk-boundary index accounting in bed2bin_. Force the
		// .bed to be read in several memory chunks with a remainder that is not divisible
		// by the thread count; het-free input makes binarization deterministic, so the
		// chunked result must be byte-identical to the single-chunk result (a stale
		// `locusInd += excessLoci` overshoots, leaving gaps and overflowing binGenotypes_).
		const std::string chunkBedName("../tests/tmpChunked.bed");
		constexpr uint32_t chunkNind{37};
		constexpr size_t   chunkNloci{100};
		constexpr size_t   maxLociPerChunk{30};   // -> 3 full chunks + 10 remainder; 30 % nThreads(4) != 0
		writeHetFreeBed(chunkBedName, chunkNind, chunkNloci, false);

		BayesicSpace::GenoTableBin singleChunk(chunkBedName, chunkNind, logFileName, nThreads);
		BayesicSpace::GenoTableBin multiChunk(chunkBedName, chunkNind, logFileName, nThreads, BayesicSpace::MemoryParameters{0, maxLociPerChunk});
		const std::string singleFile("../tests/tmpSingleChunk.bin");
		const std::string multiFile("../tests/tmpMultiChunk.bin");
		singleChunk.saveGenoBinary(singleFile);
		multiChunk.saveGenoBinary(multiFile);

		std::fstream singleIn(singleFile, std::ios::in | std::ios::binary);
		const std::vector<char> singleBytes( (std::istreambuf_iterator<char>(singleIn)), std::istreambuf_iterator<char>() );
		singleIn.close();
		std::fstream multiIn(multiFile, std::ios::in | std::ios::binary);
		const std::vector<char> multiBytes( (std::istreambuf_iterator<char>(multiIn)), std::istreambuf_iterator<char>() );
		multiIn.close();
		std::remove( chunkBedName.c_str() ); // NOLINT
		std::remove( singleFile.c_str() );   // NOLINT
		std::remove( multiFile.c_str() );    // NOLINT

		const size_t binLocusSize{(static_cast<size_t>(chunkNind) / 8) + static_cast<size_t>( (chunkNind % 8) > 0 )};
		REQUIRE( singleBytes.size() == chunkNloci * binLocusSize );
		REQUIRE( singleBytes == multiBytes );
	}
	SECTION("allJaccardLD emits every pair exactly once under a tight budget") {
		// A small RAM budget forces many small sink flushes; combined with several threads this drives
		// jaccardBlock_ ranges that lie within a single row. Each unordered pair must still appear
		// exactly once in the output (no duplicates across flush boundaries, none dropped).
		constexpr uint32_t nLoci{397};                                   // ../tests/ind197_397.bed
		constexpr size_t expectedPairs{ static_cast<size_t>(nLoci) * (nLoci - 1) / 2 };
		constexpr size_t tightBudget{12000};                            // just above the ~9925-byte table -> tiny residual
		BayesicSpace::GenoTableBin tightGTB(inputBedName, nIndividuals, logFileName, nThreads, BayesicSpace::MemoryParameters{tightBudget});
		const std::string ldFileName("../tests/tmpTightLD.tsv");
		BayesicSpace::InOutFileNames outAndBim{};
		outAndBim.outputFileName = ldFileName;
		tightGTB.allJaccardLD(outAndBim);

		std::fstream ldIn(ldFileName, std::ios::in);
		std::string ldLine;
		std::getline(ldIn, ldLine);                                     // header
		size_t totalLines{0};
		std::set< std::pair<uint32_t, uint32_t> > uniquePairs;
		while ( std::getline(ldIn, ldLine) ) {
			std::stringstream lineStream(ldLine);
			uint32_t row{0};
			uint32_t col{0};
			lineStream >> row >> col;
			uniquePairs.emplace(row, col);
			++totalLines;
		}
		ldIn.close();
		std::remove( ldFileName.c_str() ); // NOLINT
		REQUIRE( totalLines == expectedPairs );          // no duplicated lines
		REQUIRE( uniquePairs.size() == expectedPairs );  // and every pair present
	}
}

TEST_CASE("GenoTableHash methods work", "[gtHash]") {
	const std::string logFileName("../tests/binTest.log");
	const std::string inputBedName("../tests/ind197_397.bed");
	constexpr uint32_t nIndividuals{197};
	constexpr uint16_t kSketches{29};
	constexpr size_t nRowsPerBand{5};
	constexpr float invKlowBound{0.05};
	constexpr float invKhighBound{0.95};
	constexpr uint32_t lowCountMin{5000};
	constexpr uint32_t highCountMin{1000};
	constexpr size_t nThreads{4};
	constexpr uint32_t nLoci{397};
	constexpr size_t totNpairs{static_cast<size_t>(nLoci) * (static_cast<size_t>(nLoci) - 1) / 2};
	const std::string alleleCountsFile("../tests/alleleCounts.txt");
	std::fstream inAlleleCounts;
	std::string eachLine;
	inAlleleCounts.open(alleleCountsFile, std::ios::in);
	std::vector<int> macVector;
	while ( std::getline(inAlleleCounts, eachLine) ) {
		macVector.push_back( std::stoi(eachLine) );
	}
	inAlleleCounts.close();
	constexpr BayesicSpace::IndividualAndSketchCounts sketchParameters{nIndividuals, kSketches};
	SECTION("Failed GenoTableHash constructors") {
		constexpr size_t smallNind{1};
		constexpr size_t smallSketch{1};
		constexpr size_t kGtN{nIndividuals + 1};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{smallNind, smallSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{nIndividuals, smallSketch}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be at least three") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, BayesicSpace::IndividualAndSketchCounts{nIndividuals, kGtN}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be smaller than the number of individuals") );
		const std::string absentFileName("../tests/noSuchFile.bed");
		const std::string noLociFile("../tests/threeByte.bed");
		const std::string wrongMagicBytes("../tests/wrongMB.bed");
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(absentFileName, sketchParameters, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: failed to open file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(noLociFile, sketchParameters, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: no genotype records in file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(wrongMagicBytes, sketchParameters, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: first magic byte in input .bed file") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(inputBedName, sketchParameters, nThreads, logFileName, BayesicSpace::MemoryParameters{1}),
				Catch::Matchers::StartsWith("ERROR: the genotype table does not fit within the memory budget") );
		const std::vector<int> smallMACvec(13, 0);
		const std::vector<int> emptyMACvec{};
		constexpr size_t undivNind{5};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(emptyMACvec, sketchParameters, logFileName),
				Catch::Matchers::StartsWith("ERROR: empty vector of minor allele counts") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(smallMACvec, BayesicSpace::IndividualAndSketchCounts{smallNind, kSketches}, logFileName),
				Catch::Matchers::StartsWith("ERROR: number of individuals must be greater than 1") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(smallMACvec, BayesicSpace::IndividualAndSketchCounts{undivNind, kSketches}, logFileName),
				Catch::Matchers::StartsWith("ERROR: length of allele count vector") );
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(macVector, BayesicSpace::IndividualAndSketchCounts{nIndividuals, kGtN}, nThreads, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch number must be smaller than the number of individuals") );
		// the minor-allele-count constructor has its own (distinct) sketch-count checks
		// fewer than three sketches is rejected (a size-13 vector implies 13 individuals, so the count is divisible)
		const std::vector<int> macVec13(13, 0);
		constexpr size_t tinySketch{2};
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(macVec13, BayesicSpace::IndividualAndSketchCounts{13, tinySketch}, logFileName),
				Catch::Matchers::StartsWith("ERROR: sketch size must be at least three") );
		// a sketch count that reaches the empty-bin sentinel (uint16_t max) implies too large a sketch size
		constexpr uint32_t sentinelNind{65535};
		constexpr uint16_t sentinelSketch{65535}; // == emptyBinToken_
		const std::vector<int> sentinelMACvec(sentinelNind, 0);
		REQUIRE_THROWS_WITH( BayesicSpace::GenoTableHash(sentinelMACvec, BayesicSpace::IndividualAndSketchCounts{sentinelNind, sentinelSketch}, logFileName),
				Catch::Matchers::StartsWith("ERROR: Number of sketches") );
	}
	SECTION("The log is flushed to file on destruction") {
		const std::string slfLogName("../tests/saveLogTest.log");
		std::remove( slfLogName.c_str() ); // start from a clean slate // NOLINT
		{
			BayesicSpace::GenoTableHash logHSH(inputBedName, sketchParameters, nThreads, slfLogName);
		} // destructor flushes the accumulated log here
		std::fstream logIn(slfLogName, std::ios::in);
		REQUIRE( logIn.good() );
		const std::string logContents( (std::istreambuf_iterator<char>(logIn)), std::istreambuf_iterator<char>() );
		logIn.close();
		std::remove( slfLogName.c_str() ); // NOLINT
		// the constructor accumulates log messages, so the saved file must be non-empty
		REQUIRE( !logContents.empty() );
	}
	SECTION("An empty log file name disables log writing") {
		const std::string slfLogName("../tests/saveLogTestEmpty.log");
		std::remove( slfLogName.c_str() ); // start from a clean slate // NOLINT
		{
			BayesicSpace::GenoTableHash logHSH(inputBedName, sketchParameters, nThreads, std::string());
		} // destructor must not create a file
		std::fstream logIn(slfLogName, std::ios::in);
		REQUIRE( !logIn.good() );
	}
	SECTION("GenoTableHash .bed file constructor and methods with correct data") {
		BayesicSpace::GenoTableHash bedHSH(inputBedName, sketchParameters, nThreads, logFileName);
		const std::string tmpJacFile("../tests/tmpJac.tsv");
		BayesicSpace::InOutFileNames tmpFileGrp{};
		tmpFileGrp.outputFileName = tmpJacFile;
		tmpFileGrp.inputFileName  = "";
		constexpr size_t forcedChunks{3};
		constexpr float cutOff{0.0};
		bedHSH.allHashLD(cutOff, tmpFileGrp, forcedChunks);
		std::fstream hashLDfile(tmpJacFile, std::ios::in);
		std::vector<float> jaccValues; 
		std::string line;
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.emplace_back( stof(field) );
		}
		hashLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE(jaccValues.size() == totNpairs);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKlowBound](float value) {return value <= invKlowBound;}
			 ) >= lowCountMin
		);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKhighBound](float value) {return value >= invKhighBound;}
			) >= highCountMin
		);

		// a non-zero cutoff drops pairs below it
		constexpr float nonZeroCutOff{0.5F};
		bedHSH.allHashLD(nonZeroCutOff, tmpFileGrp, forcedChunks);
		hashLDfile.open(tmpJacFile, std::ios::in);
		jaccValues.clear();
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.emplace_back( stof(field) );
		}
		hashLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE(jaccValues.size() < totNpairs); // some pairs are below the cutoff
		REQUIRE(std::all_of(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&nonZeroCutOff](float value) {return value >= nonZeroCutOff;}
			)
		);
		// pairs at or above the high bound (>= invKhighBound > nonZeroCutOff) must be retained
		REQUIRE(jaccValues.size() >= highCountMin);

		std::vector<BayesicSpace::HashGroup> groups{bedHSH.makeLDgroups(nRowsPerBand)};
		REQUIRE(std::is_sorted(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
					return grp1.cumulativeNpairs < grp2.cumulativeNpairs;
				}
			)
		);
		REQUIRE(std::is_sorted(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
					return grp1.locusIndexes.at(0) < grp2.locusIndexes.at(0);
				}
			)
		);
		REQUIRE(std::all_of(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &eachGroup) {
					return std::is_sorted( eachGroup.locusIndexes.cbegin(), eachGroup.locusIndexes.cend() );
				}
			)
		);

		// file-output makeLDgroups overload, base-1 indexes (no .bim file).
		// each call re-seeds the hash, so structural invariants are checked rather than equality with the vector above
		const std::string grpFileName("../tests/tmpGroups.tsv");
		BayesicSpace::InOutFileNames grpFiles{};
		grpFiles.outputFileName = grpFileName;
		grpFiles.inputFileName  = "";
		bedHSH.makeLDgroups(nRowsPerBand, grpFiles);
		std::fstream grpFile(grpFileName, std::ios::in);
		std::string grpLine;
		std::getline(grpFile, grpLine); // header
		REQUIRE(grpLine == "groupID\tlocusIdx");
		std::vector<uint32_t> fileGroupIDs;
		std::vector< std::vector<uint32_t> > fileGroups; // locus indexes per group, in file order
		std::vector<uint32_t> allFileIdx;
		bool allRowsPrefixedG{true};
		while ( std::getline(grpFile, grpLine) ) {
			std::stringstream lineStream(grpLine);
			std::string idField;
			std::string idxField;
			lineStream >> idField;
			lineStream >> idxField;
			allRowsPrefixedG = allRowsPrefixedG && (idField.front() == 'G');
			const auto groupID{static_cast<uint32_t>( std::stoul( idField.substr(1) ) )};
			const auto locusIdx{static_cast<uint32_t>( std::stoul(idxField) )};
			if ( fileGroupIDs.empty() || (fileGroupIDs.back() != groupID) ) {
				fileGroupIDs.push_back(groupID);
				fileGroups.emplace_back();
			}
			fileGroups.back().push_back(locusIdx);
			allFileIdx.push_back(locusIdx);
		}
		grpFile.close();
		std::remove( grpFileName.c_str() ); // NOLINT

		REQUIRE( !fileGroups.empty() );
		REQUIRE( allRowsPrefixedG );
		// group IDs are a contiguous, 1-based sequence
		std::vector<uint32_t> expectedIDs( fileGroupIDs.size() );
		std::iota( expectedIDs.begin(), expectedIDs.end(), 1U );
		REQUIRE( std::equal( fileGroupIDs.cbegin(), fileGroupIDs.cend(), expectedIDs.cbegin() ) );
		// every group has at least two loci, with strictly increasing indexes
		REQUIRE(std::all_of(
				fileGroups.cbegin(),
				fileGroups.cend(),
				[](const std::vector<uint32_t> &grp) {
					return ( grp.size() >= 2 )
						&& ( std::adjacent_find( grp.cbegin(), grp.cend(),
								[](uint32_t lhs, uint32_t rhs){ return lhs >= rhs; } ) == grp.cend() );
				}
			)
		);
		// groups are ordered by their first locus index
		REQUIRE(std::is_sorted(
				fileGroups.cbegin(),
				fileGroups.cend(),
				[](const std::vector<uint32_t> &grpOne, const std::vector<uint32_t> &grpTwo){ return grpOne.front() < grpTwo.front(); }
			)
		);
		// indexes are base-1 and within range
		REQUIRE( *std::min_element( allFileIdx.cbegin(), allFileIdx.cend() ) >= 1U );
		REQUIRE( *std::max_element( allFileIdx.cbegin(), allFileIdx.cend() ) <= nLoci );

		// the same overload with a .bim file emits locus names instead of indexes
		const std::string grpBimFile("../tests/ind197_397.bim");
		grpFiles.inputFileName = grpBimFile;
		bedHSH.makeLDgroups(nRowsPerBand, grpFiles);
		const std::vector<std::string> grpLocusNames{BayesicSpace::getLocusNames(grpBimFile)};
		grpFile.open(grpFileName, std::ios::in);
		std::getline(grpFile, grpLine); // header
		REQUIRE(grpLine == "groupID\tlocusIdx");
		size_t namedRows{0};
		bool allNamesKnown{true};
		while ( std::getline(grpFile, grpLine) ) {
			std::stringstream lineStream(grpLine);
			std::string idField;
			std::string nameField;
			lineStream >> idField;
			lineStream >> nameField;
			allNamesKnown = allNamesKnown
				&& ( std::find( grpLocusNames.cbegin(), grpLocusNames.cend(), nameField ) != grpLocusNames.cend() );
			++namedRows;
		}
		grpFile.close();
		std::remove( grpFileName.c_str() ); // NOLINT
		REQUIRE( namedRows >= 2 );
		REQUIRE( allNamesKnown );

		constexpr float grpCutOff{0.75};
		BayesicSpace::SparsityParameters sparsity{};
		sparsity.similarityCutOff = grpCutOff;
		sparsity.nRowsPerBand     = nRowsPerBand;
		std::string smFileName("../tests/smTest.tsv");
		tmpFileGrp.outputFileName = smFileName;
		tmpFileGrp.inputFileName  = "";
		bedHSH.ldInGroups(sparsity, tmpFileGrp, forcedChunks);

		std::vector< std::pair<uint32_t, uint32_t> > locusPairs;
		std::vector<float> ldValues;
		std::fstream grpLDfile(smFileName, std::ios::in);
		std::getline(grpLDfile, line);             // get rid of the header
		while ( std::getline(grpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			std::pair<uint32_t, uint32_t> curLocusPair{};
			lineStream >> field;
			curLocusPair.first = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			curLocusPair.second = stoi(field) - 1; // the saved indexes are base-1
			locusPairs.emplace_back(curLocusPair);
			lineStream >> field;
			ldValues.push_back( stof(field) );
		}
		grpLDfile.close();
		std::remove( smFileName.c_str() ); // NOLINT
		REQUIRE(std::all_of(
				locusPairs.cbegin(),
				locusPairs.cend(),
				[](const std::pair<uint32_t, uint32_t> &eachPair) {
					return eachPair.first > eachPair.second;
				}
			)
		);
		REQUIRE(std::is_sorted(
				locusPairs.cbegin(),
				locusPairs.cend(),
				[](const std::pair<uint32_t, uint32_t> &pairOne, const std::pair<uint32_t, uint32_t> &pairTwo) {
					return pairOne.first < pairTwo.first;
				}
			)
		);
		REQUIRE(std::all_of(
				ldValues.cbegin(),
				ldValues.cend(),
				[&grpCutOff](const float &eachLDval) {
					return eachLDval >= grpCutOff;
				}
			)
		);
	}
	SECTION("GenoTableHash mac vector constructor and methods with correct data") {
		BayesicSpace::GenoTableHash vecHSH(macVector, sketchParameters, nThreads, logFileName);
		const std::string tmpJacFile("../tests/tmpJac.tsv");
		BayesicSpace::InOutFileNames tmpFileGrp{};
		tmpFileGrp.outputFileName = tmpJacFile;
		tmpFileGrp.inputFileName  = "";
		constexpr size_t forcedChunks{3};
		constexpr float cutOff{0.0};
		vecHSH.allHashLD(cutOff, tmpFileGrp, forcedChunks);
		std::fstream hashLDfile(tmpJacFile, std::ios::in);
		std::vector<float> jaccValues; 
		std::string line;
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.emplace_back( stof(field) );
		}
		hashLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE(jaccValues.size() == totNpairs);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKlowBound](float value) {return value <= invKlowBound;}
			 ) >= lowCountMin
		);
		REQUIRE(std::count_if(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&invKhighBound](float value) {return value >= invKhighBound;}
			) >= highCountMin
		);

		// a non-zero cutoff drops pairs below it
		constexpr float nonZeroCutOff{0.5F};
		vecHSH.allHashLD(nonZeroCutOff, tmpFileGrp, forcedChunks);
		hashLDfile.open(tmpJacFile, std::ios::in);
		jaccValues.clear();
		std::getline(hashLDfile, line);             // get rid of the header
		while ( std::getline(hashLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			lineStream >> field;
			lineStream >> field;
			lineStream >> field;
			jaccValues.emplace_back( stof(field) );
		}
		hashLDfile.close();
		std::remove( tmpJacFile.c_str() ); // NOLINT
		REQUIRE(jaccValues.size() < totNpairs); // some pairs are below the cutoff
		REQUIRE(std::all_of(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[&nonZeroCutOff](float value) {return value >= nonZeroCutOff;}
			)
		);
		// pairs at or above the high bound (>= invKhighBound > nonZeroCutOff) must be retained
		REQUIRE(jaccValues.size() >= highCountMin);

		std::vector<BayesicSpace::HashGroup> groups{vecHSH.makeLDgroups(nRowsPerBand)};
		REQUIRE(std::is_sorted(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
					return grp1.cumulativeNpairs < grp2.cumulativeNpairs;
				}
			)
		);
		REQUIRE(std::is_sorted(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &grp1, const BayesicSpace::HashGroup &grp2) {
					return grp1.locusIndexes.at(0) < grp2.locusIndexes.at(0);
				}
			)
		);
		REQUIRE(std::all_of(
				groups.cbegin(),
				groups.cend(),
				[](const BayesicSpace::HashGroup &eachGroup) {
					return std::is_sorted( eachGroup.locusIndexes.cbegin(), eachGroup.locusIndexes.cend() );
				}
			)
		);
		constexpr float grpCutOff{0.75};
		BayesicSpace::SparsityParameters sparsity{};
		sparsity.similarityCutOff = grpCutOff;
		sparsity.nRowsPerBand     = nRowsPerBand;
		std::string smFileName("../tests/smTest.tsv");
		tmpFileGrp.outputFileName = smFileName;
		tmpFileGrp.inputFileName  = "";
		vecHSH.ldInGroups(sparsity, tmpFileGrp, forcedChunks);

		std::vector< std::pair<uint32_t, uint32_t> > locusPairs;
		std::vector<float> ldValues;
		std::fstream grpLDfile(smFileName, std::ios::in);
		std::getline(grpLDfile, line);             // get rid of the header
		while ( std::getline(grpLDfile, line) ) {
			std::stringstream lineStream;
			lineStream.str(line);
			std::string field;
			std::pair<uint32_t, uint32_t> curLocusPair{};
			lineStream >> field;
			curLocusPair.first = stoi(field) - 1; // the saved indexes are base-1
			lineStream >> field;
			curLocusPair.second = stoi(field) - 1; // the saved indexes are base-1
			locusPairs.emplace_back(curLocusPair);
			lineStream >> field;
			ldValues.push_back( stof(field) );
		}
		grpLDfile.close();
		std::remove( smFileName.c_str() ); // NOLINT
		REQUIRE(std::all_of(
				locusPairs.cbegin(),
				locusPairs.cend(),
				[](const std::pair<uint32_t, uint32_t> &eachPair) {
					return eachPair.first > eachPair.second;
				}
			)
		);
		REQUIRE(std::is_sorted(
				locusPairs.cbegin(),
				locusPairs.cend(),
				[](const std::pair<uint32_t, uint32_t> &pairOne, const std::pair<uint32_t, uint32_t> &pairTwo) {
					return pairOne.first < pairTwo.first;
				}
			)
		);
		REQUIRE(std::all_of(
				ldValues.cbegin(),
				ldValues.cend(),
				[&grpCutOff](const float &eachLDval) {
					return eachLDval >= grpCutOff;
				}
			)
		);
	}
	SECTION("Chunked .bed reading is correct in the hash path") {
		// Regression test for the chunk-boundary index accounting in bed2oph_ (the OPH
		// analogue of the bed2bin_ fix). The OPH permutation is re-randomized per
		// construction, so we cannot compare against a single-chunk run; instead every
		// locus is made identical, so a correct read yields OPH-Jaccard == 1.0 for every
		// pair. A stale `locusInd += excessLoci` skips loci across chunk boundaries,
		// leaving them at the emptyBinToken_ default, which drops their pairs well below
		// 1.0 (and overflows sketches_).
		const std::string chunkBedName("../tests/tmpChunkedHash.bed");
		constexpr uint32_t chunkNind{37};
		constexpr size_t   chunkNloci{100};
		constexpr uint16_t chunkSketches{5};
		constexpr size_t   maxLociPerChunk{30};   // -> 3 full chunks + 10 remainder
		writeHetFreeBed(chunkBedName, chunkNind, chunkNloci, true);   // every locus identical
		constexpr BayesicSpace::IndividualAndSketchCounts chunkParams{chunkNind, chunkSketches};

		BayesicSpace::GenoTableHash multiChunkHash(chunkBedName, chunkParams, nThreads, logFileName, BayesicSpace::MemoryParameters{0, maxLociPerChunk});
		const std::string tmpJacFile("../tests/tmpChunkedHashJac.tsv");
		BayesicSpace::InOutFileNames outNames{};
		outNames.outputFileName = tmpJacFile;
		outNames.inputFileName  = "";
		constexpr float  zeroCutOff{0.0};
		constexpr size_t forcedChunks{2};
		multiChunkHash.allHashLD(zeroCutOff, outNames, forcedChunks);

		std::fstream jacIn(tmpJacFile, std::ios::in);
		REQUIRE( jacIn.good() );
		std::vector<float> jaccValues;
		std::string line;
		std::getline(jacIn, line);   // discard the header
		while ( std::getline(jacIn, line) ) {
			std::stringstream lineStream(line);
			std::string field;
			lineStream >> field >> field >> field;   // locus1, locus2, jaccard
			jaccValues.emplace_back( std::stof(field) );
		}
		jacIn.close();
		std::remove( chunkBedName.c_str() ); // NOLINT
		std::remove( tmpJacFile.c_str() );   // NOLINT
		// at a zero cutoff every unique locus pair is reported
		REQUIRE( jaccValues.size() == chunkNloci * (chunkNloci - 1) / 2 );
		// every locus is identical, so every pair must be a perfect match; a chunk-boundary
		// gap would leave some loci unwritten and pull their pairs below 1.0
		REQUIRE( std::all_of(
				jaccValues.cbegin(),
				jaccValues.cend(),
				[](const float value) { return value >= 1.0F - FPREC; }
			)
		);
	}
}

TEST_CASE("VashLog works", "[VashLog]") {
	const std::string logName("../tests/testVashLog.log");
	const std::string headerMsg("Test logging");

	// slurp a whole text file into a string
	auto slurp = [](const std::string &fileName) {
		std::fstream inStream(fileName, std::ios::in);
		std::stringstream buffer;
		buffer << inStream.rdbuf();
		return buffer.str();
	};

	SECTION("Header and entries are flushed on destruction") {
		{
			BayesicSpace::VashLog log( BayesicSpace::LogFileNameWithMessage{logName, headerMsg} );
			log.add("first event");
			log.add("second event");
		} // destructor flushes here
		const std::string contents{slurp(logName)};
		std::remove( logName.c_str() ); // NOLINT
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring(headerMsg + " started on") );
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring("first event")  );
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring("second event") );
	}

	SECTION("Constructor truncates a pre-existing file") {
		{
			std::fstream stale(logName, std::ios::out | std::ios::trunc);
			stale << "STALE GARBAGE\n";
		}
		{
			BayesicSpace::VashLog log( BayesicSpace::LogFileNameWithMessage{logName, headerMsg} );
			log.add("fresh");
		}
		const std::string contents{slurp(logName)};
		std::remove( logName.c_str() ); // NOLINT
		REQUIRE_THAT( contents, !Catch::Matchers::ContainsSubstring("STALE GARBAGE") );
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring("fresh") );
	}

	SECTION("Each entry carries an elapsed [m:ss] prefix and newline") {
		{
			BayesicSpace::VashLog log( BayesicSpace::LogFileNameWithMessage{logName, headerMsg} );
			log.add("hello");
		}
		std::fstream inStream(logName, std::ios::in);
		std::string line;
		std::string lastLine;
		while ( std::getline(inStream, line) ) {
			lastLine = line; // the entry is the final line after the header
		}
		inStream.close();
		std::remove( logName.c_str() ); // NOLINT
		REQUIRE( lastLine.front() == '[' );                      // bracketed [m:ss] stamp
		const auto colonPos       = lastLine.find(':');
		const auto closeBracketPos = lastLine.find("] ");
		REQUIRE( colonPos != std::string::npos );
		REQUIRE( closeBracketPos != std::string::npos );
		REQUIRE( colonPos > 1 );                                 // at least one minutes digit after '['
		REQUIRE( colonPos < closeBracketPos );
		// minutes digits sit between '[' and ':'
		const bool minutesDigits = std::all_of(
			std::next( lastLine.cbegin() ),
			std::next( lastLine.cbegin(), static_cast<std::string::difference_type>(colonPos) ),
			[](const char chr){ return std::isdigit( static_cast<unsigned char>(chr) ) != 0; }
		);
		REQUIRE( minutesDigits );
		// seconds remainder is always padded to exactly two digits
		REQUIRE( closeBracketPos - colonPos == 3 );              // ':' followed by two digits
		const bool secondsDigits = std::all_of(
			std::next( lastLine.cbegin(), static_cast<std::string::difference_type>(colonPos + 1) ),
			std::next( lastLine.cbegin(), static_cast<std::string::difference_type>(closeBracketPos) ),
			[](const char chr){ return std::isdigit( static_cast<unsigned char>(chr) ) != 0; }
		);
		REQUIRE( secondsDigits );
		REQUIRE( lastLine.substr(closeBracketPos + 2) == "hello" ); // "] " then the message, no trailing junk
	}

	SECTION("Move construction transfers ownership; the log is written exactly once") {
		{
			BayesicSpace::VashLog src( BayesicSpace::LogFileNameWithMessage{logName, headerMsg} );
			src.add("before move");
			BayesicSpace::VashLog dst( std::move(src) );
			dst.add("after move");
		} // both destructors run; only dst (toSave_ == true) must write
		const std::string contents{slurp(logName)};
		std::remove( logName.c_str() ); // NOLINT
		const auto firstHeader = contents.find("started on");
		REQUIRE( firstHeader != std::string::npos );
		// the moved-from object has toSave_ == false, so the header must not be duplicated
		REQUIRE( contents.find("started on", firstHeader + 1) == std::string::npos );
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring("before move") );
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring("after move")  );
	}

	SECTION("Move assignment transfers ownership without double-writing") {
		{
			BayesicSpace::VashLog src( BayesicSpace::LogFileNameWithMessage{logName, headerMsg} );
			src.add("alpha");
			BayesicSpace::VashLog dst;        // default-constructed: toSave_ == false, no file
			dst = std::move(src);             // dst takes over the file
			dst.add("beta");
		}
		const std::string contents{slurp(logName)};
		std::remove( logName.c_str() ); // NOLINT
		const auto firstHeader = contents.find("started on");
		REQUIRE( firstHeader != std::string::npos );
		REQUIRE( contents.find("started on", firstHeader + 1) == std::string::npos );
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring("alpha") );
		REQUIRE_THAT( contents, Catch::Matchers::ContainsSubstring("beta")  );
	}

	SECTION("Default-constructed log neither writes nor throws") {
		REQUIRE_NOTHROW( [](){
			BayesicSpace::VashLog log;        // toSave_ == false, no file opened
			log.add("orphan entry");
		}() );
	}

	SECTION("An unwritable path is handled without throwing") {
		REQUIRE_NOTHROW( [](){
			BayesicSpace::VashLog log(
				BayesicSpace::LogFileNameWithMessage{"../tests/no_such_dir/cannot_open.log", "Test logging"} );
			log.add("entry into the void");
		}() );
	}
}

