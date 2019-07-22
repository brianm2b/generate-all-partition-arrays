# GenerateAllPartitionArray

The GenerateAllPartitionArray backtracking algorithm generates a musical structure created by Milton Babbitt (1916–2011) known as an all-partition array. 
An all-parititon array is an N-by-M regular matrix, A, containing (N * M) / 12 integers each of [0, 11] that can be "covered" by the K distinct 
partitions of 12 into N parts or fewer. Each of these K partitions contains 12 distinct matrix elements [0, 11] whose N 
or fewer parts each consist of consecutive row locations in A. Parts from each of these K partitions are "packed" or ordered 
sequentially in A from left to right. The number of elements in these partitions, K * 12, required to cover A exceeds the 
number of elements in A. For this reason, an all-partition array will have (K * 12) - (N * M) overlapping locations in A 
between contiguous partitions. Such partitions may overlap in at most one row location in A for each of their contiguous 
parts. This row location in A must be the right-most element of some part in a partition k and the left-most element of 
a part in the same row in A of some contiguous partition greater than k (i.e., to its right).

## Installation

Clone or download the generate-all-partition-arrays repository to your machine and place it in a directory accessible by [Julia](https://julialang.org).

Be sure you have Julia version 1.1.0 or above.


## Usage

The GenerateAllPartitionArray algorithm can be run from the command line in the Julia REPL environment with two commands:

```julia

cd("/user/name/path/to/generate-all-partition-arrays/src")

include("GenerateAllPartitionArray.jl")

```

Note that functionality depends on the Combinatorics, DelimitedFiles, and Statistics packages, which may need to be installed in Julia prior to running.


# Authors and acknowledgment

Work for this project is largely based on publications which can be found [here](https://vbn.aau.dk/en/persons/131453/publications/). 