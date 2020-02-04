# Big-XBWT

The code is a modification of the repository [Big-BWT](https://gitlab.com/manzai/Big-BWT)
by Giovanni Manzini, to obtain something close to a XBWT of a set of reads and it's assembled genome. The modification were done by [Jan Studený](https://github.com/jendas1) and [Garance Gourdel](https://github.com/fnareoh) with the help and direction of [Giovanni Manzini](https://people.unipmn.it/manzini/) and [Travis Gagie](https://www.dal.ca/faculty/computerscience/faculty-staff/travis-gagie.html).

## Usage

```
make all
./newscan.x data/ref.in data/reads.in
./bwtparse.x data/ref.in
./pfbwt.x data/ref.in
```

The result is then in `data/ref.in.bwt`

The input format for the reference can be either:
 - The [FASTA](https://zhanglab.ccmb.med.umich.edu/FASTA/) file format with
   only one sequence: the reference genome. The filename needs to end by
   ".fasta".
 - A text file with only the genome sequence and no other character (in
   particular no `\n`), as in the example file `data/ref.in`

For the read file we support two format:
 - A text file with a line per read : first the reference position in the
  genome and then the sequence of the read, as in the example `data/reads.in`.
 - The [BAM](https://genome.sph.umich.edu/wiki/BAM) file format, but to parse
   it we depend on the [bamtools](https://github.com/pezmaster31/bamtools)
   library that you need to have installed in you home. The file name needs to
   end by ".bam" and then the usage of the `newscan` step becomes
   ```
   ./newscan_BAM_READER.x <input file for the reference> <input file for the
   reads>
   ```

## Structure

- `newscan.cpp`
- `bwtparse.cpp`
- `pfbwt.cpp`

## Description of the intermediate files

- `file.parse_old`
- `file.parse`
    - first a parse of the genome, separator (PRIME+1) and a sequence of triples, starting position of the read as 64 bit unsigned int in a parse, parse of a particular read as 64 bit ints, and a separator (PRIME+1) as 64 bit int
        - position is most likely 1-based
- `file.dict`
- `file.occ`
- `file.extended_input`
- `file.ilist`
- `file.full_children`
- `file.bwt`
