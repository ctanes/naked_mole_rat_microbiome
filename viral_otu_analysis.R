library(tidyverse)

## turn ANI-level vOTUs into a tidy table
clusters = read.table("ANI_clusters.txt", sep = '\t')

clusters_list = lapply(rownames(clusters), function(i){

    seed = clusters[i,1]
    cl = clusters[i,2]
    contigs = strsplit(cl, ",")[[1]]

    return(data.frame( species = seed, contig = contigs))
})

clusters_df = do.call("rbind", clusters_list)
write.csv(clusters_df, "tidy_species.csv", row.names = F, quote = F)

## read metadata, simplify
meta = read.table("20221213_calico_metadata.txt", sep = "\t", header = T)
meta$SampleID = gsub("\\.", "_", meta$SampleID)

meta = meta[,c("SampleID", "SampleType", "HostSpecies", "sex", "breeding_status",
                "age", "cage", "nonhost")]


## ---- get the sample ID from the contig. Maybe could have made this easier on myself ----
clusters_df$sample = sapply(clusters_df$contig, function(contig){

    if (grepl("_k", contig)){ ## genomad contigs still have this name
        return( strsplit(contig, "_k")[[1]][1]) 
    }
    else {      ## CenoteTaker3 contigs re-write the name and are harder to parse
        split_at = strsplit(contig, "@")[[1]][1]
        return( strsplit(split_at, "_(?!.*_)", perl = T)[[1]][1])
    }

})
write.table(unique(clusters_df$species), "all_votu_seeds.txt",
                row.names = F, col.names = F, quote = F)

## --- merge with metadata
votu = clusters_df[,c("species", "sample")] %>% unique()

length(unique(votu$species))

merged = merge(votu, meta, by.x = "sample", by.y = "SampleID")

blanks = merged[ merged$SampleType %in% c("DNA freewater", "Blank tip",
                                        "Elution buffer", "Empty well"),]
feces = merged[merged$SampleType == "Feces",]
write.csv(feces, "../data/feces_votu_presence.csv", row.names = F)

## ---- shared/unique molerat vs mouse ----
range = group_by(feces, species, HostSpecies) %>% summarize(n = n()) %>%
        ungroup() %>%
        spread(HostSpecies, n, fill = 0) %>%
        mutate(range = ifelse(Mouse > 0 & `Naked mole rat` > 0, "Both", 
                          ifelse(Mouse > 0, "Mouse", "NMR"))) %>%
        as.data.frame()
table(range$range)
write.table(range[range$range == "Both","species"], "shared_mouse_molerat.txt", 
            row.names = F, col.names = F, quote = F)
shared = range[range$range == "Both", "species"]

## ---- CheckV stats -----
checkv = read.table("/path/to/checkv/quality_summary.tsv", sep = "\t", header = T)
checkv = checkv[ checkv$contig_id %in% feces$species,]

table(checkv$checkv_quality)

checkv[ checkv$contig_id %in% shared,]

## ---- commonness ---- 
molerat = merged[ merged$HostSpecies == "Naked mole rat" & merged$SampleType == "Feces",]
frequency = molerat %>% group_by(species) %>%
                summarize(n_samples = n()) %>%
                ungroup() %>% as.data.frame()

pdf("../figures/n_samples.pdf")
ggplot(frequency, aes(x = n_samples)) +
    theme_classic() + 
    geom_histogram() + 
    theme(text = element_text(size = 16))
dev.off()

common = frequency %>% filter(n_samples > 30)
write.table(common$species, "../common_votus.txt", row.names = F, col.names = F, quote = F)



## ----------------------------------
breeder_range = filter(feces, breeding_status %in% c("Breeder Not Pregnant", "Breeder Pregnant")) %>%
        group_by(species, breeding_status) %>% summarize(n = n()) %>%
        ungroup() %>%
        filter(!(is.na(breeding_status))) %>%
        spread(breeding_status, n, fill = 0) %>%
        mutate(group = ifelse( `Breeder Not Pregnant` > 0, 
                        ifelse(`Breeder Pregnant` > 0, "Both", "Non-Pregnant"), "Pregnant")) %>%
        as.data.frame()
table(breeder_range$group)


## --- imgvr ---
imgvr = read.table("../ANI/contig_imgvr_ani.tsv", sep = "\t", header = T)
imgvr = imgvr[ imgvr$pid > 95 & imgvr$qcov > 85,]

nmr_feces = feces[ feces$HostSpecies == "Naked mole rat",]
species = unique(nmr_feces$species)

length(species)
sum(species %in% imgvr$qname)

mouse_feces = feces[ feces$HostSpecies == "Mouse",]
mouse_species = unique(mouse_feces$species)

length(mouse_species)
sum(mouse_species %in% imgvr$qname)

sum(mouse_species %in% species)
sum( intersect(mouse_species, species) %in% imgvr$qname)

## ------ checkv length and quality ----
checkv$log_length = log10(checkv$contig_length)

pdf("../figures/length.pdf", height = 3, width = 6)
ggplot(checkv, aes(x = log_length, fill = checkv_quality)) + 
    theme_classic() + 
    geom_histogram(color = "grey10", linewidth = .5) + 
    theme(text = element_text(size = 10))
dev.off()

good_q = checkv[ checkv$checkv_quality %in% c("Complete", "High-quality", "Medium-quality"),"contig_id"]

nmr_good = intersect(species, good_q)
mouse_good = intersect(mouse_species, good_q)
imgvr_good = intersect(imgvr$qname, good_q)

sum( intersect(nmr_good, mouse_good) %in% imgvr_good)
sum( !(nmr_good %in% union(mouse_good, imgvr_good)))
sum( !(nmr_good %in% union(mouse_good, imgvr_good))) / length(nmr_good)
sum( !(mouse_good %in% union(nmr_good, imgvr_good))) / length(mouse_good)








