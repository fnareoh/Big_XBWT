import sys,random
from Bio import SeqIO

def extract_random(ref_seq):
    for _ in range(100):
        pos = random.randint(0,len(ref_seq)-100)
        read = ref_seq[pos:pos+100]
        print(pos, read)

def main():
    if len(sys.argv) !=2:
        print("Usage: %s ref.fasta " % sys.argv[0],)
        return 1
    input_file = sys.argv[1]
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    for fasta in fasta_sequences:
        extract_random(str(fasta.seq))
        break
    return 0

if __name__ == '__main__':
    main()
