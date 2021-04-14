import sys

name = sys.argv[1]
in_file = open(name, "r")
out_file = open(name + ".upper_case","w")

if name.split('.')[-1] == "fasta":
    fasta = in_file.readline()
    out_file.write(fasta+'\n')

while True:
    c = in_file.read(1)
    if not c:
        break
    if ('A' <= c <= 'Z') or c =='\n':
        out_file.write(c)
    elif ('a' <= c <= 'z'):
        out_file.write(c.upper())
    else:
        print('non alphabet char :',c)

out_file.close()
in_file.close()