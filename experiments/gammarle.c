/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   gammarle.c
   ver 1.0 26-jan-21
   Estimate cost of rle encoding
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
  void count(FILE *, int);
  extern char *optarg;
  extern int optind, opterr, optopt;
  char *fnam;
  FILE *f;
  int c,num_opt,skip=-1;

  /* ------------- read options from command line ----------- */
  num_opt = opterr = 0;
  while ((c=getopt(argc, argv, "s:")) != -1) {
    switch (c) 
      {
      case 's':
        skip=atoi(optarg); break;
      case '?':
        fprintf(stderr,"Unknown option: %c\n", optopt);
        exit(1);
      }
    num_opt++;
  }
 
 if(optind<argc)
    fnam=argv[optind];
  else {
    fprintf(stderr, "Usage:  %s [-s a] filename\n",argv[0]);
    fprintf(stderr,"\t-s a         ascii code to skip, def none\n\n");    
    return 0;
  }

  /* -------- read input file ------------- */
   if (! (f=fopen(fnam, "rb"))) {
      perror(fnam);
      return 1;
   }

   count(f,skip);

   fclose(f);
   return 0;
}

int bitsize(unsigned long n)
{
  assert(n>0);
  int tot=0;
  while(n>0) {
    tot += 1;
    n = n/2;
  }
  return tot;
}

int gammaenc(unsigned long n)
{
  assert(n>0);
  int size = bitsize(n);
  return size + size-1;
}
  

void count(FILE *f,int skip)
{
  unsigned long count[256],tot=0,runs=0,cur=0,bits=0;
  int c,last=257,distinct=0;

  // ---- clear counters
  for(c=0;c<256;c++) 
    count[c]=0;
  

  // --- count
  while(1) {
    c = getc(f);
    if(c==EOF) break;
    // update counts for all chars
    if(count[c]++ == 0) distinct+=1;
    tot++;
    if(c==skip) continue; // do not update runs for skip
    if(c!=last) {         // end of current run
      if(cur>0) {         // not first run: there is something to encode
        runs += 1;
        bits += gammaenc(cur);
      }
      last= c;
      cur = 1;
    }
    else cur += 1;
  }
  // pending run 
  runs += 1;
  bits += gammaenc(cur);
  
  // --- print stats
  printf("# of input symbols:  %11lu\n",tot);
  for(c=0; c<256; c++)
    if(count[c]) {
      printf("%3d  %2x  ", c, c);
      if(c>31 && c<127)
        putchar(c);
      else
        putchar(' ');
      printf("  %20ld\n",count[c]);
    }
  // remove skip from total
  tot -= (skip<0) ? 0 : count[skip];
  distinct -= (skip<0 || count[skip]==0) ? 0 : 1;
  printf("# of output symbols: %11lu\n",tot);
  printf("Number of runs: %lu\n",runs);
  printf("Distinct symbols: %d\n",distinct);
  int bitxsymbol = (distinct==1) ? 0 : bitsize(distinct-1);
  unsigned long outbytes = (7+bits + bitxsymbol*runs)/8;
  printf("Gamma code: %lu (%.2f%%)\n",outbytes,100.0*outbytes/tot);
}



// ----- this function prints any char in a readable form
void pretty_putchar(int c)
{
  
  if(c>=32 && c<127)      // printable char
    printf("  %c", c);
  else if(c=='\n')
    printf(" \\n");        // \n
  else if(c=='\t')
    printf(" \\t");        // \t
  else     
    printf(" %02x", c);      // print hex code
}



