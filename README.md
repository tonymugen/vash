# Overview

A C++14 library and software to efficiently summarize genetic polymorphism statistics from [binary variant files](https://www.cog-genomics.org/plink/1.9/input#bed). Linkage disequilibrium (LD) between pairs of loci is estimated using Locality-Sensitive Hashing ([LSH](https://en.wikipedia.org/wiki/Locality-sensitive_hashing)), in particular a variant of the one-permutation hash [Mai _et al._, 2020](https://auai.org/uai2019/proceedings/papers/302.pdf). These hashes can be used to efficiently group loci by similarity in hash tables. [Jaccard similarity](https://en.wikipedia.org/wiki/Jaccard_index) is used as the LD statistic because this is the measure approximated by LSH collision probability.

The gist of the method is that each locus (or individual, although the latter functionality is not yet available) data are converted to one-bit encoding by setting the major allele (as well as missing genotype values) to 0, minor allele to 1, and heterozygotes to 1 with 50% probability. The resulting binary data are hashed, the hashes banded (as described, e.g., in the [Mining massive data sets book](http://www.mmds.org/), Chapter 3.4), and loci assigned to buckets in an unordered hash table. Loci end up in the same bucket if at least one hash band is identical between them. LD is then estimated within buckets only. Thus, if high-LD pairs are rare in a data set, many fewer than \f$N(N - 1)/2 \f$ pairwise LD values are estimated.

# Dependencies

Building the library and binaries requires a C++14 compiler. The build process requires `cmake` version 3.11 or later, but everything can be compiled by hand if desired. Current implementation requires `x86_64` processors and has only been tested on Linux. I plan to add 64-bit ARM (e.g., Apple M processor) support.

# Download and install

The repository comes with a submodule, so to clone use

```sh
git clone --recurse-submodules https://github.com/tonymugen/vash
```
Next, create a build directory

```sh
cd vash
mkdir build
```
Finally, run `cmake` to build and install the software

```sh
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .
```
Installation may require root privileges.

# Use ldblocks

`lblocks` is a command line tool that estimates LD among loci in a `plink` `.bed` file. Running it without command line flags prints the flags and their possible values. Most flags are self-explanatory, but setting values to some of them requires special consideration.

    --hash-size       - bigger values result in better estimates of LD, but at the expense of speed.
        Based on some experimentation, minimal recommended value is 40, and 100 seems sufficient.
        Setting this flag to 0 leads to full (not hash-based) Jaccard estimates among all pairs.

    --n-rows-per-band - is the banding parameter that controls sparsity.
        Setting it high (1/5 the hash size or more) leads to lower probability of
        inclusion of low to moderately similar pairs.
        Setting this flag 0 leads to all pairwise estates to be calculated.
        Beware, since this can result in huge output files.
        The software does its best to not use up RAM, but this has only been tested on Linux.

    --only-groups     - no LD calculations are performed and only groups
        and the loci they contain is saved to a file.
    --add-locus-names - locus names from the corresponding `.bim` file are outputinstead
        instead the default (base-1) locus indexes. This may result in larger output files.

Output files are tab-delimited and include group IDs, locus pair indexes, and Jaccard similarity estimates. If full Jaccard estimates are produced, group IDs are not included, but \f$r^2\f$ as well as Jaccard similarity estimates are output.

Running the software on whole genomes with millions of loci should not tax RAM (the software keeps track of free memory and only uses about half available RAM), but can still tax disk space. In addition, some pairs can be assigned to more than one group. Removal of these duplicates requires in-memory operations, so if partial results are written to disk some of these duplicates are retained. I recommend using the `--only-groups` flags in preliminary runs to get a sense of the number of locus pairs that will result given a set of parameters. Analyzing a single chromosome at a time also speeds up the analyses.

Library interface documentation can be found [here](https://www.bayesicresearch.org/softwareDocs/vash/html/index.html).
