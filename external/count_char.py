import sys

dico = {}
for i in range(1,len(sys.argv)):
    file = open(sys.argv[i], "r")
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
print(dico)
tot=0
for v in dico.values():
    tot+=v
print(tot)
