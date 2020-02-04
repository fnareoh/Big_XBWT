REF=(data/small_ref.in data/medium_ref.in data/ref.in)
READS=(data/small_reads.in data/medium_reads.in data/reads.in)


make all
g++ test/bwt.cpp -o test/xbwt_of_reference.x
for i in ${!REF[@]}; do
  ./newscan_extended.x ${REF[$i]} ${READS[$i]}
  ./bwtparse.x ${REF[$i]}
  ./pfbwt.x ${REF[$i]}
  ./test/xbwt_of_reference.x ${REF[$i]}.extended_input
  python test/compare.py ${REF[$i]}
done
