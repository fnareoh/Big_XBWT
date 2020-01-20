# Big-XBWT

The code is a modification of the repository [Big-BWT](https://gitlab.com/manzai/Big-BWT)
by Giovanni Manzini.

## Usage

```
make newscanNT.x
./newscanNT.x data/ref.in data/reads.in
# For bwtparse
mkdir build; cd build; cmake ..; make 
```

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
   end by ".bam" and then the usage becomes
   ```
   make newscanNT_BAM_READER.x
   ./newscanNT_BAM_READER.x <input file for the reference> <input file for the
   reads>
   ```

## Description of the intermediate files

- file.parse
    - first a parse of the genome, separator (PRIME+1) and a sequence of triples, starting position of the read as 64 bit unsigned int in a parse, parse of a particular read as 64 bit ints, and a separator (PRIME+1) as 64 bit int
        - position is most likely 1-based
