import pysam
import math
import sys
import matplotlib.pyplot as plt

path = ""
name_fasta, name_sam = "", ""
if len(sys.argv) == 3:
    name_fasta = sys.argv[1]
    name_sam = sys.argv[2]
elif len(sys.argv) == 2:
    path = sys.argv[1]
    name_fasta = path + "ref.fasta"
    name_sam = path + "aln-se.sam.filtered"

fasta_file = open(name_fasta, "r")
id_fasta = fasta_file.readline()
line = fasta_file.readline()
ref = ""
while line:
    ref += line.strip()
    assert line[0] != ">"  # A single sequence
    line = fasta_file.readline()
fasta_file.close()

# we assume this bamfile does not contain read match to the reversed strand
samfile = pysam.AlignmentFile(name_sam, "rb")


def round_half_up(n, decimals=0):
    multiplier = 10 ** decimals
    return math.floor(n * multiplier + 0.5) / multiplier


def test_print(read, ref):
    print(read.pos)
    print(read.seq)
    print(ref[read.pos : read.pos + len(read.seq)])
    input()


nb_read = 0
size_read = 0
sum_first_error = 0
position_first_error = [0] * 252
mapq_stats = [0] * 61
nb_err = 0
nb_reads_no_err = 0
for read in samfile.fetch():
    # assert read.flag != 16  # The read must be map forward
    if (
        read.flag != 16
        and read.pos != -1
        and read.seq != ""
        and read.seq != "*"
        and read.pos + len(read.seq) < len(ref)
        # and read.mapq == 60
        # and not ("N" in read.seq)
    ):
        nb_read += 1
        size_read = max(size_read, len(read.seq))
        if position_first_error == []:
            position_first_error = [0] * (size_read + 1)
        error = False
        for i in range(len(read.seq)):
            if read.seq[i] != ref[read.pos + i]:
                if not error:
                    error = True
                    # print(i)
                    # test_print(read, ref)
                    if read.mapq < 61:
                        mapq_stats[read.mapq] += 1
                    else:
                        print(read.mapq)
                    sum_first_error += len(read.seq) - i
                    position_first_error[i] += 1
                else:
                    nb_err += 1
                # break
        if not error:
            nb_reads_no_err += 1
            position_first_error[len(read.seq)] += 1
samfile.close()
print(f"Number of reads (forward only): {nb_read}")
print(f"Length or read: {size_read}")
print(f"coverage: {round_half_up(nb_read*size_read/len(ref),2)}")
print(
    f"Average distance from the first sequencing error position to the end: {sum_first_error}/{nb_read} = {round_half_up(sum_first_error/nb_read,2)}"
)
print(
    f"Reads without any sequencing err: {nb_reads_no_err}/{nb_read} = {round_half_up(nb_reads_no_err*100/nb_read,2)}%"
)
print(
    f"Percentage of position that are error: {round_half_up(nb_err/(nb_read*size_read),2)}%"
)
# print(f"Number of read with first position i:\n {position_first_error}")

plt.plot(position_first_error[:size_read])
plt.xlabel("Position of the first sequencing error")
plt.ylabel("Number of reads")
plt.title("Plot of the position of the first sequencing error if there is one")
if path != "":
    plt.savefig(path + "fig" + "_position_first_error.png")

# plt.clf()
# plt.plot(mapq_stats)
# plt.xlabel("Mapq of the read")
# plt.ylabel("Number of reads")
# plt.savefig("mapq.png")
