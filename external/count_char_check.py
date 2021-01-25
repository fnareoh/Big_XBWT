import sys
import struct

assert(len(sys.argv)==3)
name_file = [sys.argv[1],sys.argv[2]]
name_file_bwt = [sys.argv[1]+".bwt"]
name_end = sys.argv[1]+".is_end"

def count(name_files):
    dico = {}
    for name in name_files:
        file = open(name, "r")
        if name.split('.')[-1] == "fasta":
            fasta = file.readline()
            print(fasta)
        print(name, end= " ")
        while True:
            c = file.read(1)
            if not c:
                break
            if c < 'A' or c > 'Z':
                pass
            elif c not in dico:
                dico[c]=1
            else:
                dico[c]+=1
    print()
    print(dico)
    tot=0
    for v in dico.values():
        tot+=v
    return tot

def count_rle(name_files):
    dico = {}
    for name in name_files:
        file = open(name, "r")
        bwt = file.read().split("\n")[:-1]
        print(name, end= " ")
        for s in bwt:
            c = s[0]
            nb = int(s[1:])
            if c < 'A' or c > 'Z':
                pass
            elif c not in dico:
                dico[c]=nb
            else:
                dico[c]+=nb
    print()
    print(dico)
    tot=0
    for v in dico.values():
        tot+=v
    return tot

def bits(f):
    while True:
        b = f.read(1)
        if not b: break
        b = ord(b)
        for i in reversed(range(8)):
            yield b >> i & 1

def count_is_end(name):
    tot_is_end = 0
    nb_char = 0
    with open(name, "rb") as file:
         for bit in bits(file):
             #print(bit, end="")
             tot_is_end += bit
             nb_char +=1
    print()
    print("Number of 1 in is_end vector: ",tot_is_end)
    #print(nb_char)


tot = count(name_file)
tot_bwt = count_rle(name_file_bwt)
count_is_end(name_end)

print(tot,"=?",tot_bwt)
if (tot== tot_bwt):
    sys.exit(0)
else :
    sys.exit(20)
