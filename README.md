# Big-XBWT

The code is a modification of the repository [Big-BWT](https://gitlab.com/manzai/Big-BWT)
by Giovanni Manzini, that given a genome and an alligned set of read to this genome, creates a tree structure and build the XBWT of that tree that indexes the reads and the genome.
 The modification were done by [Garance Gourdel](https://github.com/fnareoh) with the help and direction of [Giovanni Manzini](https://people.unipmn.it/manzini/) and [Travis Gagie](https://www.dal.ca/faculty/computerscience/faculty-staff/travis-gagie.html) and the participation of [Jan Studený](https://github.com/jendas1).

## Dependencies


## Usage

```
make all
./scan_BAM_READER.x data/ref.fasta data/reads.bam
./bwtparse.x data/ref.fasta
./pfbwt64.x data/ref.fasta
```

The result is then in `data/ref.fasta.bwt`

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
