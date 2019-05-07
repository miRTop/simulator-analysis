library(tidyverse)
# LATEST
###########
# precision
###########
parser = . %>% 
    set_names(c("class", "hit", "isomir", "counts")) %>%
    mutate(isomir=ifelse(isomir=="", "5prime_non", isomir)) %>% 
    mutate(isomir=ifelse(isomir=="3prime_addition", "3prime_non", isomir)) %>% 
    group_by(isomir) %>% 
    mutate(pct=as.numeric(counts)/sum(as.numeric(counts))*100) %>% 
    mutate(hit=ifelse(hit=="None", "Missed",hit)) %>% 
    ungroup()

###########
#  round1
###########
r1 = list.files("data/round1/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
    }) %>% bind_rows() %>% 
    mutate(round = "isomirs")

###########
#  round2
###########

r2 = list.files("data/round2/stats/", "accuracy", full.names = T) %>% 
    lapply(., function(fn){
        name = gsub("_accuracy.tsv","",basename(fn))
        read_tsv(fn, col_names = FALSE) %>% 
            parser %>% mutate(tool=name)
        
    }) %>% bind_rows() %>% 
    mutate(round = "isomirs Comb")

bind_rows(r1,r2) %>% 
    filter(isomir!="NA") %>% 
    mutate(class=paste(hit, class)) %>% 
    ggplot(aes(paste(tool, round), pct, fill=class)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="bottom") + 
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    coord_flip() +
    facet_wrap(~isomir, nrow=1) +
    scale_fill_manual(values = c("grey", "orange", "blue", "orange4", "blue4"))  
ggsave("figures-local/round12.pdf", width=9, height = 6)


bind_rows(r1,r2) %>% 
    filter(isomir!="NA") %>% 
    mutate(class=paste(hit, class)) %>% 
    ggplot(aes(isomir, pct, fill=class)) + 
    geom_bar(stat = "identity") +
    theme(legend.position="bottom") + 
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) +
    facet_wrap(~paste(tool, round), ncol=4) +
    scale_fill_manual(values = c("grey", "orange", "blue", "orange4", "blue4"))  


###########
# counts
###########
labs= c("missed","<50%x", "<30%x", "+/-30%x", ">30%x", ">50%x")

read_files = function(folder){
    mirtop = list.files(file.path(folder, "stats"), "counts", full.names = T) %>%
        lapply(., function(fn){
            name = gsub("_counts.tsv","",basename(fn))
            read_tsv(fn, col_names = c("miRNA", "unqiue", "all")) %>% 
                mutate(tool=name) %>% 
                full_join(simulated, by=c("miRNA"="Transcript")) %>%
                select(-Count)
        }) %>% bind_rows() %>% 
        set_names(c("mirna", "unique", "all", "tool")) %>% 
        filter(mirna!="miRNA")
    
    featurecounts = list.files(file.path(folder, "features"), ".txt$", full.names = T) %>%
        lapply(., function(fn){
            name = gsub("\\.txt","",basename(fn))
            read_tsv(fn, skip = 1) %>%
                .[,c(1,7)] %>%
                set_names("mirna", "fcounts") %>%
                mutate(tool=name)
        }) %>% bind_rows()
    
    features_multiON=list.files(file.path(folder, "features_multiON"), ".txt$", full.names = T) %>%
        lapply(., function(fn){
            name = gsub("\\.txt","",basename(fn))
            read_tsv(fn, skip = 1) %>%
                .[,c(1,7)] %>%
                set_names("mirna", "fcountsM") %>%
                mutate(tool=name)
        }) %>% bind_rows()
    
    full_join(mirtop, featurecounts, by = c("mirna", "tool")) %>% 
        full_join(features_multiON, by = c("mirna", "tool")) 
    
}

plot_data = function(dt){
    dt %>% 
    mutate(unique = ifelse(is.na(unique), 0, unique),
           all = ifelse(is.na(all), 0, all),
           tool = ifelse(is.na(tool), "0missed", tool)) %>% 
        mutate(mirtopU_fc = as.numeric(unique)/Count,
               mirtopM_fc = as.numeric(all)/Count,
               fcountsU_fc = as.numeric(fcounts)/Count,
               internal_fc = as.numeric(internalc)/Count,
        ) %>% # head
        select(mirna, tool,  mirtopU_fc:internal_fc) %>% 
        gather(method, fc, -mirna, -tool) %>%
        filter(!is.na(fc), tool!="0missed") %>%
        filter(!(tool=="manatee" & grepl("mirtop", method))) %>% 
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
        theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5))
}

###########
#  round1
###########
simulated = read_tsv("input/round1/simulated_miRNA_BF_counts.tsv")

r1 = read_files("data/round1/")

own=bind_rows(
    read_tsv("data/round1/own/manatee.tsv") %>% 
        .[,c(3,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "manatee"),
    read_tsv("data/round1/own/shortstack.txt") %>%  
        .[,c(2,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "shortstack")
) %>% group_by(mirna, tool) %>% 
    summarise(internalc=sum(internalc)) %>% 
    ungroup()

full_join(r1,own, by = c("mirna", "tool")) %>% 
    full_join(simulated, by = c("mirna" = "Transcript") ) %>% 
    plot_data() + 
    ggsave("figures-local/round1.counts.pdf", width = 12)

###########
#  round2
###########
simulated = read_tsv("input/round2/simulated_miRNA_BF_counts.tsv")
r2 = read_files("data/round2/")

own=bind_rows(
    read_tsv("data/round2/own/shortstack.txt") %>%  
        .[,c(2,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "shortstack")
) %>% group_by(mirna, tool) %>% 
    summarise(internalc=sum(internalc)) %>% 
    ungroup()

#full_join(r1,own, by = c("mirna", "tool")) %>% 
    full_join(r2,simulated, by = c("mirna" = "Transcript") ) %>% 
        mutate(unique = ifelse(is.na(unique), 0, unique),
               all = ifelse(is.na(all), 0, all),
               tool = ifelse(is.na(tool), "0missed", tool)) %>% 
        mutate(mirtopU_fc = as.numeric(unique)/Count,
               mirtopM_fc = as.numeric(all)/Count,
               fcountsU_fc = as.numeric(fcounts)/Count,
#               internal_fc = as.numeric(internalc)/Count,
        ) %>% # head
        select(mirna, tool,  mirtopU_fc:fcountsU_fc) %>% 
        gather(method, fc, -mirna, -tool) %>%
        filter(!is.na(fc), tool!="0missed") %>%
        filter(!(tool=="manatee" & grepl("mirtop", method))) %>% 
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
