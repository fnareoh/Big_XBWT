import os
import sys

name_fasta = sys.argv[1]
head_tail = os.path.split(name_fasta)
fasta_file = open(name_fasta, "r")
output_file = open(head_tail[0] + "/longest_" + head_tail[1], "w")
id = fasta_file.readline()
seq = fasta_file.readline()
assert id and seq and id[0] == ">"
max_contig = seq
max_contig_id = id

while id and seq:
    assert id[0] == ">"
    if len(seq) > len(max_contig):
        max_contig = seq
        max_contig_id = id
    l = fasta_file.readline()
    if l and l[0] == ">":
        id = l
        seq = fasta_file.readline()
    elif l:
        seq += l
    else:
        id = l

fasta_file.close()
output_file.write(max_contig_id)
output_file.write(max_contig)
output_file.close()
