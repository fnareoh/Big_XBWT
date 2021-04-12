#!/usr/bin/env python3
import sys, argparse, io, struct, collections

Description = """
Tool to estimate the size of Run Lengh Encoding
using Gamma, Delta and 7x8 byte encoding.
No real compression is done, just reports the number of runs
and the cost of the encoding.
Slow: for files of hundred of MBs use the C tool gammarle"""


# generic encoder class to be initialized with a name
# and an ecoding function 
class encoder:

  def __init__(self,lenfun,name):
    self.lenfun = lenfun
    self.bits = 0
    self.ints = 0
    self.name = name

  def reset(self):
    self.bits = 0
    self.ints = 0
    
  def encode(self,n):
    self.bits += self.lenfun(n)
    self.ints += 1


# compute size of rle encoding when input file is 
# in input format with symbols+runs separated by space
# as in "A4 B23 A21 C77"  
def ascii_rle(encoders,symbols,c,args):
  tot_ibytes  = 0
  tot_runs = 0

  with open(args.infile) as f:
    for s in f:  # read one group at a time
      symbols.add(s[0])
      n = int(s[1:])
      c[s[0]] += n
      for e in encoders: # encode n 
        e.encode(n)
      tot_ibytes += n    # update total length and number runs 
      tot_runs +=1
      if tot_runs==args.r:
        break
  return tot_ibytes,tot_runs


def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('infile', help='input file name', type=str)
  parser.add_argument('-z', help="ignore character \0", action='store_true')
  parser.add_argument('-a', help="ascii input", action='store_true')
  parser.add_argument('-r', help="# runs to extract (def all)", default=-1, type=int)
  parser.add_argument('-v', help="verbose", action='store_true')
  args = parser.parse_args()
  if args.v:
    show_command_line(sys.stderr)
  # ---
  encoders = [encoder(gammal,"gamma"),encoder(deltal,"delta"),encoder(vb7x8l,"vb")]
  # --- init counters
  tot_ibytes  = 0
  tot_runs = 0
  symbols = set()  # distinct symbols 
  cur = 0          # chars in current run
  last = -1        # char of current run 
  c = collections.Counter()
  
  if args.a:
    (tot_ibytes,tot_runs) = ascii_rle(encoders,symbols,c,args)
  else:
    with open(args.infile,"rb") as f:
      assert isinstance(f, io.BufferedIOBase)
      while True:
        bb = f.read(1) # read single byte
        if(len(bb)==0): break  # end of file
        tot_ibytes +=1   # one more input byte
        if args.v and (tot_ibytes % (100*2**20))==0:
          print(tot_ibytes//2**20, "MBs read", file=sys.stderr)
        z = bb[0]        # integer in range(256)
        c[z] += 1        # update counter
        if args.z and z==0:
          continue 
        if(z!=last):     # end of current run
          if(cur>0):     # not first run there is something to encode
            if tot_runs +1 == args.r:
              break
            tot_runs += 1
            for e in encoders:
              e.encode(cur)
          symbols.add(z)
          last= z
          cur = 1
        else:
          cur += 1
    # encode last run   
    tot_runs += 1
    for e in encoders:
      e.encode(cur)

  # output 
  print("Input bytes: %19s" % tot_ibytes)
  print("Distinct bytes: ", len(symbols))
  for sym in c:
    count = c[sym]
    # sometimes sym is a char sometimes is an ascii code
    if type(sym)==str: sym = ord(sym)
    print("%3s %3s  %s  %20s" % (sym,hex(sym)[2:], chr(sym) if (sym>=32) else " ", count))
  print("Number of runs: ",tot_runs)
  print("Output size:")
  bitxsymbol = 0 if len(symbols)==1 else bitsize(len(symbols)-1)
  for e in encoders:
    assert e.ints == tot_runs, "Runs count error"
    totbits = e.bits + bitxsymbol*tot_runs
    outbytes = (totbits+7)//8
    print("%9s: %d  %.2f%%" % (e.name,outbytes,100*outbytes/tot_ibytes))
  return  
      




def show_command_line(f):
  f.write("Command line: ") 
  for x in sys.argv:
     f.write(x+" ")
  f.write("\n")   

# compute number of bits of a positive integer
def bitsize(n):
  assert n>0, "Invalid input"
  tot = 0
  while n>0:
    tot +=1
    n = n//2
  return tot

# ----------- gamma ------------------
  
# len of gamma code   
def gammal(n):  
  assert n>0, "Invalid input: " + str(n)
  x = bitsize(n)
  return 2*x -1

# gamma code  
def gamma(n):
  assert n>0, "Invalid input: " + str(n)
  x = bitsize(n)    
  s = "0"*(x-1) + bin(n)[2:]
  return s
  
# gamma code inversion  
def gammai(s):
  assert len(s)%2==1, "Invalid input: " + s
  assert s.count("0") + s.count("1") == len(s), "Invalid input: " + s
  x = (len(s)+1)//2
  n =  int(s[len(s)-x:],2)
  return n

# --- delta -------------

def deltal(n):
  assert n>0, "Invalid input: " + str(n)
  x = bitsize(n)
  return gammal(n) + x-1
  
def delta(n):
  assert n>0, "Invalid input: " + str(n)
  x = bitsize(n)
  return gamma(n) + bin(n)[3:]
  
def deltai(s):  
  assert len(s)>0, "Invalid input: " + s
  assert s.count("0") + s.count("1") == len(s), "Invalid input: " + s
  c0 = s.find("1") # number of leading 0s
  x = len(s) - (2*c0+1)
  n = int("1" + s[len(s)-x:],2)
  return n

# ---- variable length 7x8 -----------------

def vb7x8l(n):
  assert n>0, "Invalid input: " + str(n)
  x = bitsize(n)
  b = (x+6)//7   # number of bytes
  return 8*b
  
def vb7x8(n):
  assert n>=0, "Invalid input: " + str(n)
  bits = bin(n)[2:]
  s = ""
  while len(bits)>7:
    s += "0" + bits[-7:]
    bits = bits[:-7]
  s += "1" + "0"*(7-len(bits))+bits
  return s  
  
def vb7x8i(s):  
  assert len(s)%8==0, "Invalid input: " + s
  assert s.count("0") + s.count("1") == len(s), "Invalid input: " + s
  n = 0
  shift = 0
  while True:
    bits = s[:8]
    stop = (s[0]=='1')
    n += int(bits[1:],2) << shift
    shift += 7
    s = s[8:]
    if stop: break 
  assert len(s)==0
  return n


if __name__ == '__main__':
  main()


