import sys

naive = ""
prefix_free = ""

with open(sys.argv[1]+'.bwt', 'r') as file:
    prefix_free = file.read().replace('\n', '')

with open(sys.argv[1]+'.extended_input.bwt', 'r') as file:
    naive = file.read().replace('\n', '')

assert(len(naive) == len(prefix_free))

diff = [i for i in range(len(naive)) if prefix_free[i] != naive[i]]

print("=================== Comparison for", sys.argv[1])
print("mismatches indexes and chars: i naive[i] prefix_free[i]")
for i in diff:
    print(i,naive[i],prefix_free[i])
