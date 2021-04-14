REF=$1
READS=$2
make all
./scan_BAM_READER.x $REF $READS
./bwtparse.x $REF
./pfbwt64.x $REF
#python external/count_char_check.py $REF $READS.txt
