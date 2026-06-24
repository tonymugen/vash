# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`vash` is a C++ library (plus the `ldblocks` command-line tool) that summarizes genetic polymorphism statistics from PLINK [`.bed` binary variant files](https://www.cog-genomics.org/plink/1.9/input#bed). Linkage disequilibrium (LD) between locus pairs is estimated as Jaccard similarity, approximated via Locality-Sensitive Hashing (a one-permutation hash, OPH). Loci are hashed, the hashes are banded, and loci sharing a band land in the same hash-table bucket; LD is then estimated only within buckets, avoiding the full N(N-1)/2 pairwise computation when high-LD pairs are rare.

Requires `x86_64` (uses BMI/BMI2/POPCNT intrinsics — see `-mbmi -mbmi2 -mpopcnt` in `CMakeLists.txt`) and has only been tested on Linux. C++17 (`CMAKE_CXX_STANDARD` is 17; recently migrated from C++14). ARM64 support is planned but not yet available.

## Build, test, docs

CMake ≥ 3.21. The build fetches the [`bayesicUtilities`](https://github.com/tonymugen/bayesicUtilities) dependency (`brutilities`) automatically via `FetchContent`; tests additionally fetch Catch2, and these are only built when the project is top-level.

Build types are `Release` (default), `Debug`, `Profile` (adds `-p` profiling), and `Test` (`-g -O3`). All commands run from the `build/` directory.

```sh
# Release build + install (install may need root)
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .

# Tests (compiled with address+leak sanitizers enabled via BUILD_TESTS)
cmake -DCMAKE_BUILD_TYPE=Test -DBUILD_TESTS=ON ..
cmake --build .
./tests                          # run all
./tests "[gtHash]"               # run one TEST_CASE by Catch2 tag
ctest                            # run via CTest (tests are auto-discovered)

# Docs (needs Doxygen + pdflatex on PATH)
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCS=ON ..
cmake --build .
```

When `BUILD_TESTS=ON`, AddressSanitizer and LeakSanitizer are compiled into `vash`, `ldblocks`, and `tests`. Compiler warnings are aggressive (`-Wall -Wextra -Wconversion -Wpedantic -Wshadow -Wold-style-cast -Wsign-conversion`, plus more for GCC). `clang-tidy` config is in `.clang-tidy`.

## Architecture

All code lives in `namespace BayesicSpace`. Three source/header pairs plus the app:

- **`gvarHash` (`include/gvarHash.hpp`, `src/gvarHash.cpp`)** — the core. Two main classes:
  - `GenoTableBin` — one-bit binarized genotype table. Converts `.bed` data (or a minor-allele-count `std::vector<int>`) to one-bit encoding (major allele/missing → 0, minor → 1, heterozygote → 1 with 50% probability) and computes **full** pairwise Jaccard LD via `allJaccardLD`. This is the `--hash-size 0` / no-hashing path.
  - `GenoTableHash` — the OPH/LSH path. Builds OPH sketches, groups loci into LD buckets (`makeLDgroups`), and estimates LD within groups (`ldInGroups`) or across all hash collisions (`allHashLD`). Sketch count and individual count are passed together via `IndividualAndSketchCounts`.
  - Both classes are move-only (copy deleted), take a log-file name, parallelize across threads (default `std::thread::hardware_concurrency()`), and process data in memory-bounded chunks (`BedDataStats`, `suggestNchunks`) — the code tracks available RAM and aims to use ~half of it so whole-genome inputs don't exhaust memory. Many small POD structs (`LocationWithLength`, `SparsityParameters`, `InOutFileNames`, `HashGroup`, etc.) are used to pass grouped parameters.

- **`similarityMatrix` (`include/similarityMatrix.hpp`, `src/similarityMatrix.cpp`)** — `SimilarityMatrix`, a compact representation of a (possibly sparse) similarity matrix. Jaccard values (0.0–1.0) are **quantized into 256 bins** (`uint8_t`), and only included pairs are stored. Index helpers convert between `RowColIdx` (row/column) and the vectorized full triangular-matrix index (`FullIdxValue`). `JaccardPair` holds intersection/union counts. Hashing blocks return `SimilarityMatrix` objects that are merged.

- **`vashFunctions` (`include/vashFunctions.hpp`, `src/vashFunctions.cpp`)** — free helper functions used by the classes: bit counting (`countSetBits`), MurMurHash3 mixer/finalizer, `.bed`/`.bim` file parsing and magic-byte checks, thread range partitioning, and `getAvailableRAM` (reads `procfs`, falls back to 2 GiB).

- **`apps/ldblocks.cpp`** — CLI wrapper. Parses string/int flags into maps and dispatches to the appropriate `GenoTable*` constructor and method. Key flags: `--hash-size` (0 = full Jaccard, no hashing), `--n-rows-per-band` (banding/sparsity; 0 = all pairs), `--only-groups` (emit groups, skip LD), `--add-locus-names` (use `.bim` names instead of base-1 indexes). Running with no flags prints usage.

The `random.cpp` source from `brutilities` is compiled directly into the `vash` library target.

## Tests

`tests/tests.cpp` is a single Catch2 file organized by `TEST_CASE` with tags (`[countSetBits]`, `[MurMurHash]`, `[bedData]`, `[SimilarityMatrix]`, `[gtBin]`, `[gtHash]`) and nested `SECTION`s. `tests/alleleCounts.txt` is test input data.
