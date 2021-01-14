#script to extract a test instance at random from ref.in and check that it runs
#and restitute the correct number of chars
while true ; do
  python external/generate_tests.py data/ecoli/ref.in 50000 100 1000
  REF=data/ecoli/ref.in.ref
  READS=data/ecoli/ref.in.read
  make
  ./scan.x $REF $READS
  ./bwtparse.x $REF
  ./pfbwt.x $REF
  if ! python external/count_char_check.py $REF $READS;
  then exit;
  fi
done