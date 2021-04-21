# Big-XBWT

The code is a modification of the repository [Big-BWT](https://gitlab.com/manzai/Big-BWT)
by Giovanni Manzini. Given a genome and an aligned set of read to this genome, it creates a tree structure and build the XBWT of that tree which indexes the reads and the genome.
 The modification were done by [Garance Gourdel](https://github.com/fnareoh) with the help and direction of [Giovanni Manzini](https://people.unipmn.it/manzini/) and [Travis Gagie](https://www.dal.ca/faculty/computerscience/faculty-staff/travis-gagie.html) and the participation of [Jan Studený](https://github.com/jendas1).

## Dependencies

Those two library are needed:

* **Bamtools** a C++ library to handle [BAM](https://genome.sph.umich.edu/wiki/BAM) file [[Github](https://github.com/pezmaster31/bamtools)][[Wiki](https://github.com/pezmaster31/bamtools/wiki)].
* **SDSL** a succint data structure library for C++, used for sparse bitvector [[Github](https://github.com/simongog/sdsl-lite)].

After installing those library you might need to modify the path to those library in the makefile.

## Usage

To compute the structure for a genome stored in [FASTA](https://zhanglab.ccmb.med.umich.edu/FASTA/) file format and reads in [BAM](https://genome.sph.umich.edu/wiki/BAM) file format, the instructions are the following:

```
make all
./scan_BAM_READER.x data/ref.fasta data/reads.bam
./bwtparse.x data/ref.fasta
./pfbwt64.x data/ref.fasta
```

The result is then in `data/ref.fasta.bwt`. This will also add numerous temporary file to the path where the reference was. For convenience, a simple script `run_bam.sh` regroups those 3 steps.

## Example

```
./run.sh data/toy_example/ref.in data/toy_example/reads.in
```

Computes the xbwt of a synthetic dataset, then we check that the occurrences of each char is the same in the xbwt and in the original input.

## Input format

The input format for **the reference file** can be either:
 - The [FASTA](https://zhanglab.ccmb.med.umich.edu/FASTA/) file format with
   only one sequence: the reference genome. The filename needs to end by
   ".fasta".
 - A text file with only the genome sequence and no other character (in
   particular no `\n`), as in the example file `data/ref.in`

For **the read file** we support two format:
 - A text file with a line per read : first the reference position in the
  genome and then the sequence of the read  `data/reads.in`.
  ```
  make
  ./scan.x data/ref.in data/reads.in
  ./bwtparse.x data/ref.in
  ./pfbwt.x data/ref.in
  ```
 - The [BAM](https://genome.sph.umich.edu/wiki/BAM) file format, but to parse
   it we depend on the [bamtools](https://github.com/pezmaster31/bamtools)
   library that you need to have installed in you home. The file name extension has to be ".bam".

## Reproducibility

To reproduce our experiments you will need to install:

* [Samtools](https://github.com/samtools/samtools) a software to handle SAM and BAM files.
* [BWA](https://github.com/lh3/bwa) a software package for mapping DNA sequences against a large reference genome.
* [ropebwt2](https://github.com/lh3/ropebwt2) a tool for constructing the FM-index for a collection of DNA sequences.
* [SPRING](https://github.com/shubhamchandak94/Spring) a compression tool for Fastq files.

To use our experimental pipeline you will then need to modify the path in `experiments/makefile` to fit your installation path. Then usage of that pipeline is given a path to a genome in fasta format `data/ecoli/ref.fasta` and a path to the reads in fastq format `data/ecoli/reads.fastq`:
```
cd experiments
make DIR=../data/ecoli NAME=ecoli
cat ecoli.runs_summary
```
