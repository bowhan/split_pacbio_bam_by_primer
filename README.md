# Split PacBio subreads BAM file using a custom smrtbell primer sequence

## Pre-requisite

    1. Linux environment
    2. C++ compiler that supports C++11
    3. Boost `apt-get install -y libboost-all-dev`
    4. [Anaconda](https://www.anaconda.com/)
    5. pbbam `conda install -c bioconda pbbam`

### Install

```bash
git clone --recursive git@github.com:bowhan/split_pacbio_bam_by_primer.git
mkdir -p split_pacbio_bam_by_primer/build && cd split_pacbio_bam_by_primer/build
PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:~/anaconda3/lib/pkgconfig/ cmake -DCMAKE_BUILD_TYPE=RELEASE .. && make -j4 
make install || sudo make install
```

### Usage

Shell pipeline usage:

```bash
# You can provide arbitrary number of inputs from either RS II or Seqeuel
# in the following example,
# A01_1 is the primary analysis output of RS II
# sample.subreads.bam is the primary output of Sequel

split_pacbio_bam.sh split \
    -i A01_1 \
    -i sample.subreads.bam \
    -p ATCTCTCTCAATTTTTTTTTTTTTTTTTTTTTTTAAGAGAGAGAT \
    -o output_dir
```

C++ program usage:

```bash
# -p AT...AT: custom smrtbell sequence
# -t 8 : use 8 threads
# -o out.subreads.bam: output bam file
# test.subreads.bam: input subreads bam file

split_primer_from_pbbam \
    -p ATCTCTCTCAATTTTTTTTTTTTTTTTTTTTTTTAAGAGAGAGAT \
    -t 8 \
    -o out.subreads.bam \
    test.subreads.bam

# advanced usage:
# -m: Reads with SW-score lower than 70 will be dumped
# -f: If the SW-score difference between the two best alignments are
#     lower than 10, this read is dumped
# -M: score for match in SW alignments
# -S: penalty (positive integer) for a mismatch in SW alignments
# -O: gap open penalty (positive integer)
# -E: gap extension penalty (positive integer)
split_primer_from_pbbam \
    -p ATCTCTCTCAATTTTTTTTTTTTTTTTTTTTTTTAAGAGAGAGAT \
    -t 8 \
    -o out.subreads.bam \
    -m 70 \
    -f 10 \
    -M 2 \
    -S 1 \
    -O 3 \
    -E 1 \
    test.subreads.bam
```
