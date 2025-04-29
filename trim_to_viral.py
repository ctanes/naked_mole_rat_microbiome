#!/home/carter/mambaforge/bin/python

from Bio import SeqIO
import pandas as pd


## all file paths are changed to indicate which file, rather than actual location

# Read in full contig file
contig_path = "/path/to/contig/file.fna"
contigs = {}
for record in SeqIO.parse(contig_path, "fasta"):
    contigs[record.id] = record

out = open("checkv_trimmed_contigs.fasta", 'w')


complete_genomes = pd.read_csv("/path/to/checkv/complete_genomes.tsv", sep = '\t')
circular = complete_genomes[ complete_genomes["prediction_type"].isin(["ITR", "DTR"])]
circular_contigs = set(circular["contig_id"].tolist())

used_contigs = set()
## Go through and trim each contig based on contamination.tsv file from CheckV
with open("/path/to/checkv/contamination.tsv", 'r') as filein:
    header = True
    for line in filein:
        if header:
            header = False
        else:
            row = line.split("\t")
            contig = row[0]
            provirus = row[5] == "Yes"
            if not provirus:
##    If it's not a provirus, we just have to filter based on he viral and bacerial genes
##    >3 viral genes and <3 bacterial
                viral_genes = int(row[3])
                host_genes = int(row[4])
                if viral_genes > 3 and host_genes < 3 or contig in circular_contigs:
                    out.write(">" + contig + "\n")
                    out.write( str(contigs[contig].seq) + "\n") 
                    if contig in used_contigs:
                        print(contig)
                        assert False
                    used_contigs.add(contig)
            else:
##    if it is a provirus, we want to just look at the viral region
##
                regions = row[8].split(',')
                boundaries = row[10].split(',')
                viral_genes_region = row[12].split(',')
                host_genes_region = row[13].split(',')
                for i,r in enumerate(regions):
##
##      Check each region with the same criteria as the non-provirus regions
                    boundary = boundaries[i]
                    start,end = boundary.split('-')
                    n_viral = int(viral_genes_region[i])
                    n_host = int(host_genes_region[i])
#
                    if n_viral > 3 and n_host < 3:
            #            print(contig, n_viral, n_host)
                        ## write the contig region
                        sequence = contigs[contig].seq
                        trimmed = sequence[int(start):int(end)]
##
##     in rare cases a contig includes multiple viruses, like virus, host, virus, etc.
##     I just split these up and give the later ones letter modifiers
##
                        if contig in used_contigs:
                            print(contig)
                            new_name = contig + "_" + ["b", "c", "d", "e", "f"][i] 
                        else:
                            new_name = contig
                        used_contigs.add(contig)
                        out.write(">" + new_name + "\n")
                        out.write(str(trimmed) + "\n")


out.close()
