# Overview

A C++14 library and software to efficiently summarize genetic polymorphism statistics from [binary variant files](https://www.cog-genomics.org/plink/1.9/input#bed). Linkage disequilibrium between pairs of loci and relationships among individuals are estimated using Locality-Sensitive Hashing ([LSH](https://en.wikipedia.org/wiki/Locality-sensitive_hashing)).

The repository comes with a submodule, so to clone use

```sh
git clone --recurse-submodules https://github.com/tonymugen/vash
```

This software is still in development. Implemented features do appear to work, but testing is ongoing. The makefile is also still a work in progress. It may be necessary to directly compile the test programs. This will be fixed shortly, however.
