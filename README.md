# PSYCHIC
Code for finding putative enhancers using Hi-C data

### Usage
python htad-chain \<config file\>

For running the example use
`python htad-chain examples/himr90.chr20.conf`
from the repo directory.

#### Config file format
- `res` resolution of the Hi-C file, in bases (40000)
- `win` interaction distance cutoff in bases, usually 2000000
- `chrname` chromosome name, used for genes, size and output. should use 'chr1'..'chrX'
- `chrsize` path to bed file of chromsome lengths (examples/hg19.size.bed)
- `output\_prefix` output prefix for this conf file (hIMR90)
- `output\_dir` path in which to store the output files (examples/output)
- `input\_matrix` path to input Hi-C matrix, should be in either CSV file or Bingren & Dixon format (examples/hIMR90.chr20). _This matrix must by symmetric_.
- `genes\_file` path to bed file describing genes (examples/hg19.genes.bed)

For a functioning example config file consult `examples/himr90.chr20.conf`

### Requirements
The code was used on a Linux machine.
It has scripts in matlab, python and perl, so the minimal requirements would be - 
- matlab
- python2.7
- perl
- \*nix tools - cut, sed
- a shell that supports pushd and popd
- \*nix file paths (/example/of/a/path)

### Filesystem
- **`htad-chain.py`**
Main command line interface for PSYCHIC
- **`matlab/`**
Main matlab files
- **`domaincall_software/`**
Slightly adapted files from [Hi-C DomainCaller](http://chromosome.sdsc.edu/mouse/hi-c/download.html)
- **`examples/`**
Example files, contains config, Hi-C matrix, chromsome sizes and gene list
