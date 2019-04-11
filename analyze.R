# round 1
library(tidyverse)

source("scripts/gff_parser.R")

read_gff("data/output/razers3/simulated_miRNA_BF.tsv.gz")
read_gff("data/output/seqbuster/simulated_miRNA_BF.tsv.gz")
read_gff("data/output/bwa/bwa.simulated.tsv.gz")
read_gff("data/output/output/bowtiev1.0.1/mirtop.tsv.gz")
read_gff("data/output/output/gsnap.m10/mirtop.tsv.gz")
read_gff("data/output/output/gsnap.m3/mirtop.tsv.gz")
read_gff("data/output/output/star/mirtop.tsv.gz")
read_gff("data/output/output/shortstack/mirtop.tsv.gz")
read_gff("data/output/output/bwa/bwa.simulated.tsv.gz")
read_gff("data/output/manatee/simulated_miRNA_BF_Manatee.tsv.gz")



parse = . %>% 
    select(isomir, correct_mirna, n, nts) %>% 
    unite("uid", isomir, nts, sep = "#") %>% 
    spread(correct_mirna, n) %>% 
    separate(uid, c("isomirs", "nts"), sep = "#") %>% 
    rowwise() %>% 
    mutate(true_positive_rate = correct/(correct+incorrect), 
           mapping_rate = 1-missed/(correct+incorrect+missed)) %>% 
    select(-correct, -incorrect, -missed)

df = bind_rows(
    read_csv("data/output/razers3/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "razers3"),
    read_csv("data/output/seqbuster/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "seqbuster"),
    read_csv("data/bwa/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "bwa"),
    read_csv("data/output/output/bowtiev1.0.1/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "bowtiev1.0.1"),
    read_csv("data/output/output/gsnap.m10/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "gsnap.m10"),
    read_csv("data/output/output/gsnap.m10/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "gsnap.m3"),
    read_csv("data/output/output/star/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "star"),
    read_csv("data/output/output/shortstack/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "shortstack"),
    read_csv("data/output/manatee/parsed.tsv") %>% 
        parse() %>%  
        mutate(tool = "manatee")
) 

df %>% 
    gather(metric, accuracy, -isomirs, -nts, -tool) %>% 
    ggplot(aes(x = nts, y = accuracy, fill = metric)) +
        geom_bar(stat = "identity", position = "dodge") +
        facet_grid(tool~isomirs, scales = "free_x") +
    theme_minimal() +
    theme(legend.position="bottom", strip.text = element_text(size=7)) +
    scale_fill_manual(values = c("black", "orange")) +
    ggsave("figures/summary.pdf", width = 12, height = 9)


# round 2

# precision
parser = . %>% 
    set_names(c("class", "hit", "isomir", "counts")) %>%
    mutate(isomir=ifelse(isomir=="", "5prime_non", isomir)) %>% 
    mutate(isomir=ifelse(isomir=="3prime_addition", "3prime_non", isomir)) %>% 
    group_by(isomir) %>% 
    mutate(pct=as.numeric(counts)/sum(as.numeric(counts))*100) %>% 
    mutate(hit=ifelse(hit=="None", "Missed",hit)) %>% 
    ungroup()

list.files("data/round2/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
}) %>% bind_rows() %>% 
    filter(isomir!="NA") %>% 
    mutate(class=paste(hit, class)) %>% 
    ggplot(aes(isomir, pct, fill=class)) + 
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    facet_wrap(~tool) +
    scale_fill_manual(values = c("grey", "orange", "blue", "orange4", "blue4"))  +
    ggsave("figures-local/round2.pdf", width=9, height = 6)

# counts
labs= c("missed","<50%x", "<30%x", "+/-30%x", ">30%x", ">50%x")
simulated = read_tsv("input/round2/simulated_miRNA_BF_counts.tsv")
mirtop = list.files("data/round2/stats/", "counts", full.names = T) %>%
    lapply(., function(fn){
        name = gsub("_counts.tsv","",basename(fn))
        read_tsv(fn) %>% 
            mutate(tool=name) %>% 
            full_join(simulated, by=c("miRNA"="Transcript")) %>%
            select(-Count)
    }) %>% bind_rows() %>% 
    set_names(c("mirna", "unique", "all", "tool")) %>% 
    filter(mirna!="miRNA")

featurecounts = list.files("data/round2/features", ".txt$", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("\\.txt","",basename(fn))
        read_tsv(fn, skip = 1) %>% 
            .[,c(1,7)] %>% 
            set_names("mirna", "fcounts") %>% 
            mutate(tool=name)
    }) %>% bind_rows()

features_multiON=list.files("data/round2/features_multiON", ".txt$", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("\\.txt","",basename(fn))
        read_tsv(fn, skip = 1) %>% 
            .[,c(1,7)] %>% 
            set_names("mirna", "fcountsM") %>% 
            mutate(tool=name)
    }) %>% bind_rows()

features_multiON_fractionON = list.files("data/round2/features_multiON_fractionON", ".txt$", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("\\.txt","",basename(fn))
        read_tsv(fn, skip = 1) %>% 
            .[,c(1,7)] %>% 
            set_names("mirna", "fcountsMF") %>% 
            group_by(mirna) %>% 
            summarise(fcountsMF=sum(fcountsMF)) %>% 
            ungroup() %>% 
            mutate(tool=name)
    }) %>% bind_rows()

full_join(mirtop, featurecounts, by = c("mirna", "tool")) %>% 
    full_join(features_multiON, by = c("mirna", "tool")) %>% 
    full_join(features_multiON_fractionON, by = c("mirna", "tool")) %>% 
    full_join(simulated, by = c("mirna" = "Transcript") ) %>% 
    mutate(unique = ifelse(is.na(unique), 0, unique),
           all = ifelse(is.na(all), 0, all),
           tool = ifelse(is.na(tool), "0missed", tool)) %>% 
    mutate(mirtopU_fc = as.numeric(unique)/Count,
           mirtopA_fc = as.numeric(all)/Count,
           fcountsU_fc = as.numeric(fcounts)/Count,
           fcountsM_fc = as.numeric(fcountsM)/Count,
           fcountsMF_fc = as.numeric(fcountsMF)/Count
           ) %>% # head
    select(mirna, tool,  mirtopU_fc:fcountsMF_fc) %>% 
#    filter(tool=="shortstack", mirna=="hsa-let-7a-5p")
    gather(method, fc, -mirna, -tool) %>%
    filter(!is.na(fc), tool!="0missed") %>% 
    distinct() %>% 
    mutate(fc = cut(fc, c(-1,0.0001,0.5,0.7,1.3,2,100000), 
                    (labs)),
           fc = factor(fc, levels=rev(labs))) %>%  
    
    ggplot(aes(paste(tool, method), fill=fc)) +
    geom_bar() +
    scale_fill_manual(values=c( "#feb24c", "#edf8b1",
                               "#addd8e",
                               "#edf8b1", "#feb24c", "black")) +
    ggthemes::theme_base() +
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) + 
    ggsave("figures-local/round2.counts.pdf", width = 12)
