import sys
import struct

assert(len(sys.argv)==3)
name_file = [sys.argv[1],sys.argv[2]]
name_file_bwt = [sys.argv[1]+".bwt"]

def count(name_file):
    dico = {}
    for name in name_file:
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

def count_rle(name_file):
    dico = {}
    for name in name_file:
        file = open(name, "rb")
        print(name, end= " ")
        while True:
            bin_s = file.read(1)
            if not bin_s:
                break
            c = struct.unpack('c', bin_s)[0].decode()
            nb = struct.unpack('i', file.read(4))[0]
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

tot = count(name_file)
tot_bwt = count_rle(name_file_bwt)

print(tot,"=?",tot_bwt)
if (tot== tot_bwt):
    sys.exit(0)
else :
    sys.exit(20)