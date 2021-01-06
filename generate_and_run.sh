python external/generate_tests.py data/ref.in 30 10 5
REF=data/ref.in.ref
READS=data/ref.in.read
make
./scan.x $REF $READS
./bwtparse.x $REF
./pfbwt.x $REF
python external/count_char.py $REF $READS
python external/count_char.py $REF.bwt
