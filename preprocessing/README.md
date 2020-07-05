# Overview

Here we build experssion matrix, allele-specific count matrix, and PEER factors for each tissue.
The steps are:

1. `wasp2matrix`: convert WASP output to AS matrix.
2. `expression`: construct the expression matrix (total count matrix).
3. `peer`: compute the PEER factor of expression matrix.
