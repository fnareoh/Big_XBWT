import sys
import bisect
from collections import defaultdict

assert(len(sys.argv)==3)
pattern = sys.argv[1]
name_file_bwt = sys.argv[2]

file_bwt = open(name_file_bwt, "r")
bwt = file_bwt.read().split(" ")[:-1]
nb_char = [0]*(len(bwt)+1)
seen_char_bwt = [defaultdict(lambda: 0)]
occ = defaultdict(lambda: 0)

for i in range(len(bwt)):
    c = bwt[i][0]
    nb = int(bwt[i][1:])
    bwt[i] = (c,nb)
    nb_char[i+1] = nb_char[i]+ nb
    seen_char_bwt.append(seen_char_bwt[i].copy())
    if c in seen_char_bwt[i+1]:
        seen_char_bwt[i+1][c]+=nb
    else:
        seen_char_bwt[i+1][c]=nb
    if c in occ:
        occ[c]+= nb
    else:
        occ[c]=nb

tot = 0
F = {}
for key in sorted(occ):
    F[key] = (tot,tot+occ[key])
    tot += occ[key]

def compute_char_before(pos):
    if pos==nb_char[-1]:
        return seen_char_bwt[-1]
    i_prev_run = bisect.bisect(nb_char,pos)-1
    seen = seen_char_bwt[i_prev_run].copy()
    c,nb = bwt[i_prev_run]
    dist_to_start_run = pos - nb_char[i_prev_run]
    if c in seen :
        seen[c] += dist_to_start_run
    else :
        seen[c] = dist_to_start_run
    return seen

interval_F = F[pattern[-1]]
for i in range(len(pattern)-2,-1,-1):
    a,b = interval_F
    seen_a = compute_char_before(a)
    seen_b = compute_char_before(b)
    start_letter_in_F = F[pattern[i]][0]
    interval_F = (start_letter_in_F + seen_a[pattern[i]],
    start_letter_in_F + seen_b[pattern[i]])

a,b = interval_F
print("counting for "+pattern+":", b-a+1)