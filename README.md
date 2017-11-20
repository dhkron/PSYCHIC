# PSYCHIC
Code for finding putative enhancers using Hi-C data

### Usage
python htad-chain.py `\<config file\>`

For running the example use
`python htad-chain.py examples/himr90.chr20.conf`
from the repo directory.

#### Config file format
- `res` resolution of the Hi-C file, in bases (40000)
- `win` interaction distance cutoff in bases, usually 2000000
- `chrname` chromosome name, used for genes, size and output. should use 'chr1'..'chrX'
- `chrsize` path to bed file of chromsome lengths (examples/hg19.size.bed)
- `output_prefix` output prefix for this conf file (hIMR90)
- `output_dir` path in which to store the output files (examples/output)
- `input_matrix` path to input Hi-C matrix for the chromosome (see format specifications below)
- `genes_file` path to bed file describing genes (examples/hg19.genes.bed)

For a functioning example config file consult `examples/himr90.chr20.conf`

### Requirements
The code was used on a Linux machine.
It has scripts in matlab, python and perl, so the minimal requirements would be - 
- Matlab
- python2.7
- Perl
- Unix tools - cut, sed, pushd, popd (typically installed by default)

### Filesystem
- **`htad-chain.py`**
Main command line interface for PSYCHIC
- **`matlab/`**
Main matlab files
- **`domaincall_software/`**
Slightly adapted files from [Hi-C DomainCaller](http://chromosome.sdsc.edu/mouse/hi-c/download.html)
- **`insulation/`**
Domain caller from Crane et al. and additional scripts [Insulation Score DomainCaller](https://github.com/dekkerlab/crane-nature-2015)
- **`examples/`**
Example files, contains config, Hi-C matrix, chromsome sizes and gene list

### Input matrix
`input_matrix` should be in a _symmetric_ csv or tab-delimitered file, specifying the Hi-C data (for a given chomromsome). The first column could be either empty or contain the names of each genomic segment (in a fixed stride, as specified in the `res` variable).
Then, each cell (\< i,j \>) should contain the number of contacts between the matching segments in the chromosome. These data could be previously normalized to account for various Hi-C biases, and is assumed to be _symmetric_. See example under `examples/hIMR90.chr20.matrix.txt`.

### Output files
