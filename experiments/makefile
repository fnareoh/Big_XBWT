# pipeline for the comparison of the number of runs
# using several eBWT variants and the XBWT

# directory containing the input files
DIR=~/Documents/Big_XBWT/data/hydrophilia
NAME=hydrophilia_raw_ctg
# use ln -s to create the links ref.fasta -> reference  reads.fastq -> reads
REF=ref.fasta
READS=reads.fastq

# executables used in the pipeline
SAMTOOLS=samtools
BWA=~/Documents/bwa/bwa
ROPE=~/Documents/ropebwt2/ropebwt2
SPRING=~/Documents/SPRING/build/spring-reorder
XBWTDIR=~/Documents/Big_XBWT
TXT2FA=$(XBWTDIR)/experiments/txt2fa.py
XBWTEXE=run_bam.sh
RLE1=$(XBWTDIR)/experiments/multirle.py
REV=$(XBWTDIR)/experiments/sam_filter_reversed.py
GAMMARLE=$(XBWTDIR)/experiments/gammarle

# used internally no need to change
XBWT=xref.fasta.bwt
BAM=aln-se.bam
SAM=aln-se.sam



all: $(DIR)/$(NAME).runs_summary

$(GAMMARLE):
	gcc $(GAMMARLE).c -o $(GAMMARLE)

# compute xwbt statistics
$(DIR)/xbwt.runs: $(DIR)/$(XBWT)
	$(RLE1) -a $< > $@

# compute ebwt statistics: reads ordering
$(DIR)/ebwt.runs: $(DIR)/$(BAM).fasta
	$(ROPE) -R -o tmp.bwt $^
	$(GAMMARLE) tmp.bwt       > $@
	$(GAMMARLE) tmp.bwt -s 36 >> $@

# compute ebwt statistics: RLO ordering
$(DIR)/ebwtRLO.runs: $(DIR)/$(BAM).fasta
	$(ROPE) -R -s -o tmpRLO.bwt $^
	$(GAMMARLE) tmpRLO.bwt       > $@
	$(GAMMARLE) tmpRLO.bwt -s 36 >> $@

# compute ebwt statistics: SPRING ordering
$(DIR)/ebwtSP.runs: $(DIR)/$(BAM).fasta
	$(SPRING) --fasta-input -i $(DIR)/$(BAM).fasta -o tmp.spring.fasta
	$(ROPE) -R -o tmpSP.bwt tmp.spring.fasta
	$(GAMMARLE) tmpSP.bwt       > $@
	$(GAMMARLE) tmpSP.bwt -s 36 >> $@

# create summary
$(DIR)/$(NAME).runs_summary: $(DIR)/ebwtSP.runs $(DIR)/ebwtRLO.runs $(DIR)/ebwt.runs $(DIR)/xbwt.runs
	echo "-------- ebwt" >> $(DIR)/$(NAME).runs_summary
	cat $(DIR)/ebwt.runs >> $(DIR)/$(NAME).runs_summary
	echo "-------- ebwt SPRING" >> $(DIR)/$(NAME).runs_summary
	cat $(DIR)/ebwtSP.runs >> $(DIR)/$(NAME).runs_summary
	echo "-------- ebwt RLO" >> $(DIR)/$(NAME).runs_summary
	cat $(DIR)/ebwtRLO.runs >> $(DIR)/$(NAME).runs_summary
	echo "-------- xbwt" >> $(DIR)/$(NAME).runs_summary
	cat $(DIR)/xbwt.runs >> $(DIR)/$(NAME).runs_summary


# create fasta file containing aligned reads from bam.txt file
$(DIR)/$(BAM).fasta: $(DIR)/$(BAM).txt
	$(TXT2FA) $< $@

# compute bam.txt and xbwt file from bam and xref
$(DIR)/$(BAM).txt: $(DIR)/$(BAM) $(DIR)/x$(REF)
	cd $(XBWTDIR) && ./$(XBWTEXE) $(DIR)/x$(REF) $(DIR)/$(BAM)

# make xbwt dependant from bam.txt
$(DIR)/$(XBWT): $(DIR)/$(BAM).txt

# create xref from ref
$(DIR)/x$(REF): $(DIR)/$(REF)
		ln -s $(DIR)/$(REF) $(DIR)/x$(REF)

# compute bam file from sam
$(DIR)/$(BAM): $(DIR)/$(SAM)
	python $(REV) $<
	$(SAMTOOLS)  view -S -b $<.filtered > $@
	#$(SAMTOOLS)  view -S -b $< > $@

# compute sam file from ref and reads
$(DIR)/$(SAM): $(DIR)/$(REF) $(DIR)/$(READS)
	$(BWA) index $(DIR)/$(REF)
	$(BWA) mem $(DIR)/$(REF) $(DIR)/$(READS) > $(DIR)/$(SAM)


clean:
	rm -fr gammarle tmp.bwt tmpRLO.bwt tmpSP.bwt tmp.spring.fasta $(DIR)/*runs $(DIR)/$(BAM).fasta

