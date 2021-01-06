import random
import sys

ref_file = sys.argv[1]
size_ref = int(sys.argv[2])
size_read = int(sys.argv[3])
nb_read = int(sys.argv[4])

file = open(ref_file, "r")
entire_ref = file.read()
file.close()

start_ref = random.randint(0,len(entire_ref)-size_ref-1)
ref = entire_ref[start_ref : start_ref+size_ref]

file = open(ref_file+".ref", "w")
file.write(ref)
file.close()

file = open(ref_file+".read", "w")
for i in range(nb_read):
    start_read = random.randint(0,len(ref)-size_read-1)
    read = ref[start_read: start_read+size_read]
    print(start_read, read)
    file.write(str(start_read)+ " " +read)
    if i != nb_read-1 :
        file.write("\n")
file.close()