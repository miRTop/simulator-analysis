library(tidyverse)
# LATEST
###########
# precision
###########
parser = . %>% 
    set_names(c("class", "hit", "isomir", "counts")) %>%
#    mutate(isomir=ifelse(isomir=="", "5prime_non", isomir)) %>% 
#    mutate(isomir=ifelse(isomir=="3prime_addition", "3prime_non", isomir)) %>%
    separate(isomir, into = c("isomir", "type", "size"), sep="-", fill = "right") %>% 
    mutate(size=ifelse(is.na(size), "1", size)) %>% 
    group_by(isomir, type) %>% 
    mutate(pct=as.numeric(counts)/sum(as.numeric(counts))*100) %>% 
    mutate(hit=ifelse(hit=="None", "Missed",hit)) %>% 
    ungroup()

###########
#  round0
###########
r0 = list.files("data/output.S0/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
    }) %>% bind_rows() %>% 
    mutate(round = "cannonical")


###########
#  round1
###########
r1 = list.files("data/output.S1/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
    }) %>% bind_rows() %>% 
    mutate(round = "isomirs")

###########
#  round2
###########

r2 = list.files("data/output.S2/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
    }) %>% bind_rows() %>% 
    mutate(round = "isomirs Comb")


###########
#  round3
###########

r3 = list.files("data/output.S3/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
    }) %>% bind_rows() %>% 
    mutate(round = "sRNA")

###########
#  round4
###########

r4 = list.files("data/output.S4/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
    }) %>% bind_rows() %>% 
    mutate(round = "isomirs Comb + sRNA")


#############################################################################

full = bind_rows(r0,r1,r2,r3,r4) %>% 
#    filter(tool!="bowtie1_default", isomir!="other", tool!="gsnap_m10") %>% 
#    filter(isomir!="NA") %>% 
    mutate(round = factor(round, levels = c("cannonical","isomirs", "isomirs Comb", "isomirs Comb + sRNA"))) %>% 
    mutate(class=paste(hit, class)) %>% 
    mutate(tool_round = paste(tool, round),
           tool_round = factor(tool_round, levels = rev(unique(tool_round)))) %>% 
    mutate(is_comb = ifelse(grepl(":", isomir), "yes", "no"))

full %>% filter(round=="isomirs")  %>%
#    filter(!(isomir %in% c("snp", "indel"))) %>%
    mutate(isomir=paste(isomir,type)) %>% 
    group_by(tool, isomir, size) %>% 
    mutate(total = sum(pct),
           pct = pct / total * 100) %>% 
    ggplot(aes(tool, pct, fill=class)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="bottom") + 
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    facet_grid(size~isomir) +
    scale_fill_manual(values = c("grey", "orange", "blue", "orange4", "blue4")) +
    ggsave("results/S1-accuracy-nocomb-bysize.pdf", width=12, height = 6)

full %>% filter(is_comb=="no", round!="sRNA")  %>% 
    mutate(isomir=paste(isomir,type)) %>% 
    ggplot(aes(tool_round, pct, fill=class)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="bottom") + 
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    coord_flip() +
    facet_wrap(~isomir, nrow=1) +
    scale_fill_manual(values = c("grey", "orange", "blue", "orange4", "blue4")) +
    ggsave("results/accuracy-nocomb.pdf", width=12, height = 6)

full %>% filter(is_comb=="yes") %>% 
    mutate(isomir=paste(isomir,type)) %>% 
    mutate(isomir=gsub(":", "\n",isomir)) %>% 
    ggplot(aes(tool_round, pct, fill=class)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="bottom") + 
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    coord_flip() +
    facet_wrap(~isomir, nrow=1) +
    scale_fill_manual(values = c("grey", "orange", "blue", "orange4", "blue4")) +  
    ggsave("results/accuracy-comb.pdf", width=12, height = 6)

