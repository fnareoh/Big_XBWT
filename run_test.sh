# usage: ./run_test.sh at the root of the project
# it computes the bwt of the ref and reads in REF and READS
# according to a naive implementation and the prefix free parsing
# There can be differences (they are ouputed by test/compare.py)
# when the context is precisely the same, because the order then
# depends on the parsing which the naive implementation doesn't have.
# But there should be relativly few (0 for small, 2 for medium, 8 for ref.in).
REF=(data/small_ref.in data/medium_ref.in data/ref.in)
READS=(data/small_reads.in data/medium_reads.in data/reads.in)


make all
make newscan_extended.x
g++ test/bwt.cpp -o test/xbwt_of_reference.x
for i in ${!REF[@]}; do
  ./newscan_extended.x ${REF[$i]} ${READS[$i]}
  ./bwtparse.x ${REF[$i]}
  ./pfbwt.x ${REF[$i]}
  ./test/xbwt_of_reference.x ${REF[$i]}.extended_input
  python test/compare.py ${REF[$i]}
done
