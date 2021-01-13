import pybam
import sys

f = open(sys.argv[1]+".txt", "w")
for alignment in pybam.read(sys.argv[1]):
    s = alignment.sam_seq
    if (alignment.sam_pos0 != -1 and s !="" and s!="*"):
        f.write(s+"\n")
f.close()

