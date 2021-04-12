import pysam
import sys

""" Starting from the SAM format, filter out the reads mapped to the reversed strand and output to BAM format """

nb_filtered = 0
tot_nb_reads = 0

samfile = pysam.AlignmentFile(sys.argv[1], "rb")
forward_reads = pysam.AlignmentFile(sys.argv[1] + ".no_reverse", "wb", template=samfile)

for read in samfile.fetch():
    flag = read.flag
    tot_nb_reads += 1
    # print(flag)
    if (
        read.flag != 16
        and read.pos != -1
        and read.seq != ""
        and read.seq != "*"
        and read.pos + len(read.seq) < len(ref)
        and read.mapq == 60
        and not ("N" in read.seq)
    ):
        forward_reads.write(read)
    else:
        nb_filtered += 1

forward_reads.close()
samfile.close()

percentage_filtered = 100 * nb_filtered / tot_nb_reads
print("Number of filtered reads:", nb_filtered)
print("Percentage of filtered reads:", percentage_filtered)
