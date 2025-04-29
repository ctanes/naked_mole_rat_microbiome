#!/bin/bash

# set file paths
viral_contigs=
IMG_VR_database=
viral_contigs_blastdb=

## Blast contigs to the IMG_VR database and determine average nucleotide identity

blastn -query $viral_contigs -db $IMG_VR_database -out IMG_VR_blast.tsv -outfmt '6 std qlen slen'  -perc_identity 90 -num_threads 40

python ~/programs/MGV/ani_cluster/blastani.py -i IMG_VR_blast.tsv  -o contig_imgvr_ani.tsv

## Blast all contigs against each other for clustering into viral OTUs

blastn -query $viral_contigs -db $viral_contigs_blastdb -out  ll_v_all_blast.m6 -outfmt '6 std qlen slen'  -perc_identity 90 -num_threads 40

python ~/programs/MGV/ani_cluster/blastani.py -i all_v_all_blastn.txt -o blast_ani.tsv

~/programs/MGV/ani_cluster/cluster.py --fna $viral_contigs --ani blast_ani.tsv --out ANI_clusters --min_ani 95 --min_tcov 85


