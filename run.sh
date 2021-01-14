REF=$1
READS=$2
make
./scan.x $REF $READS
./bwtparse.x $REF
./pfbwt.x $REF
python external/count_char_check.py $REF $READS
