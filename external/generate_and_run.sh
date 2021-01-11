#script to extract a test instance at random from ref.in and check that it runs
#and restitute the correct number of chars
while true ; do
  python external/generate_tests.py data/ref.in 50000 1000 500
  REF=data/ref.in.ref
  READS=data/ref.in.read
  make
  ./scan.x $REF $READS
  ./bwtparse.x $REF
  ./pfbwt.x $REF
  if ! python external/count_char_check.py $REF $READS;
  then exit;
  fi
done