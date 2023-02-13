# Overview

A C++14 library and software to efficiently summarize genetic polymorphism statistics from [binary variant files](https://www.cog-genomics.org/plink/1.9/input#bed). Linkage disequilibrium (LD) between pairs of loci is estimated using Locality-Sensitive Hashing ([LSH](https://en.wikipedia.org/wiki/Locality-sensitive_hashing)), in particular a variant of the one-permutation hash (Mai _et al._, 2020). These hashes can be used to efficiently group loci by similarity in hash tables.

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
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
cmake --install .
```
Installation may require root privileges.

# Use `ldblocks`

`lblocks` is a command line tool that estimates LD among loci in a `plink` `.bed` file. Running it without command line flags prints the flags and their possible values.

This software is still in development. Implemented features do appear to work, but testing is ongoing.
