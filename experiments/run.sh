#Run all experiment necessary for the real world graph

make DIR=~/Documents/Big_XBWT/data/ecoli NAME=ecoli
make DIR=~/Documents/Big_XBWT/data/staphylococus NAME=staphylococus
make DIR=~/Documents/Big_XBWT/data/R_sphaeroides/HiSeq/raw NAME=R_sphaeroides_Hiseq_ctg_raw
make DIR=~/Documents/Big_XBWT/data/R_sphaeroides/MiSeq/raw NAME=R_sphaeroides_Miseq_ctg_raw
make DIR=~/Documents/Big_XBWT/data/R_sphaeroides/HiSeq/trimmed NAME=R_sphaeroides_Hiseq_ctg_trimmed
make DIR=~/Documents/Big_XBWT/data/R_sphaeroides/MiSeq/trimmed NAME=R_sphaeroides_Miseq_ctg_trimmed

python create_csv_plot.py r_sphaeroides_dataset_path.csv
python create_csv_plot.py real_world_dataset_path.csv