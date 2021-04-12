#!/usr/bin/env python3

Description = """
Convert a collection of reads (one read per line) into
fasta format adding a dummy header line before each read
"""

import sys, argparse

def show_command_line(f):
  f.write("Command line: ") 
  for x in sys.argv:
     f.write(x+" ")
  f.write("\n")   

def main():
  show_command_line(sys.stderr)
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('infile', help='input file name', type=str)
  parser.add_argument('outfile', nargs='?', help='output file name', default=None)
  parser.add_argument('-v', help="verbose", action='store_true')
  args = parser.parse_args()
  # --- init counters
  written = 0
  symbols = 0
  if args.outfile:
    g = open(args.outfile,"w")
  else:
    g = sys.stdout
  
  with open(args.infile) as f:
    for s in f:
      print(">", written, file=g)
      print(s,end="", file=g)
      written += 1
      symbols += len(s)
  print("Input lines:",written,file=sys.stderr)
  print("Total symbols (including newlines):",symbols,file=sys.stderr)  
  return        

if __name__ == '__main__':
  main()
