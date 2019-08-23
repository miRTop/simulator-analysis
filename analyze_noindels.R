library(tidyverse)

###########
# counts
###########
labs= c("missed","<50%x", "30-50%x", "10-30%x","+/-10%x", "10-30%x", "30-50%x", ">50%x")

read_files = function(folder){
    # mirtop = list.files(file.path(folder, "stats"), "counts", full.names = T) %>%
    #     lapply(., function(fn){
    #         name = gsub("_counts.tsv","",basename(fn))
    #         read_tsv(fn, col_names = c("miRNA", "unqiue", "all")) %>% 
    #             mutate(tool=name) %>% 
    #             full_join(simulated, by=c("miRNA"="Transcript")) %>%
    #             select(-Count) %>% 
    #             set_names(c("mirna", "mirtopU", "mirtopM", "tool")) %>% 
    #             filter(mirna!="miRNA")
    #     }) %>% bind_rows() 

    featurecounts = list.files(file.path(folder, "noindels","features"), ".txt$", full.names = T) %>%
        lapply(., function(fn){
            name = gsub("\\.txt","",basename(fn))
            read_tsv(fn, skip = 1) %>%
                .[,c(1,7)] %>%
                set_names("mirna", "fcountsU") %>%
                mutate(tool=name)
        }) %>% bind_rows()
    
    features_multiON=list.files(file.path(folder, "noindels","features_multiON"), ".txt$", full.names = T) %>%
        lapply(., function(fn){
            name = gsub("\\.txt","",basename(fn))
            read_tsv(fn, skip = 1) %>%
                .[,c(1,7)] %>%
                set_names("mirna", "fcountsM") %>%
                mutate(tool=name)
        }) %>% bind_rows()
    
    # if (nrow(mirtop) == 0 ){
        df = full_join(featurecounts, features_multiON, by = c("mirna", "tool")) 
    # }else{
        # df = full_join(mirtop, featurecounts, by = c("mirna", "tool")) %>% 
            # full_join(features_multiON, by = c("mirna", "tool")) 
        
    # }
    df %>% gather(method, counts, -mirna, -tool) %>% 
        filter(counts>0)
}

###########
#  round1
###########
simulated = read_tsv("input/round1/simulated_miRNA_BF_counts.tsv")
universe = as.character(simulated$Transcript)
r1 = read_files("data/round1/")


own1=bind_rows(
    read_tsv("data/round1/stats/seqbuster_counts.tsv", col_names = F) %>% 
        group_by(X1) %>% 
        summarise(internalc=sum(X3)) %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "seqbuster"),
    read_tsv("data/round1/stats/razers3_counts.tsv", col_names = F) %>% 
        group_by(X1) %>% 
        summarise(internalc=sum(X3)) %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "razers3"),
    read_tsv("data/round1/own/srnabench.grouped") %>%  
        .[,c(1,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "srnabench"),
    
    read_tsv("data/round1/own/excert/readCounts_miRNAmature_sense.txt") %>% 
        separate_rows(ReferenceID, sep="\\|") %>% 
        separate(ReferenceID, into = "mirna", extra = "drop", sep=":") %>% 
        .[,c(1,3)] %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "excert"),
    read_tsv("data/round1/own/manatee.tsv") %>% 
        .[,c(3,4)] %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "manatee"),
    read_tsv("data/round1/own/shortstack.txt") %>%  
        .[,c(2,4)] %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "shortstack")
) %>% group_by(mirna, tool) %>% 
    summarise(internalc=sum(internalc)) %>% 
    ungroup() %>% 
    full_join(simulated, by=c("mirna"="Transcript")) %>%
    select(-Count) %>% 
    mutate(internalc = ifelse(is.na(internalc), 0, internalc)) %>% 
    gather(method, counts, -mirna, -tool)

r1 = bind_rows(r1,own1) %>%
    full_join(simulated, by = c("mirna" = "Transcript") )
#     plot_data() + 
#     ggsave("figures-local/round1.counts.pdf", width = 12)

###########
#  round2
###########
simulated = read_tsv("input/round2/simulated_miRNA_BF_counts.tsv")
r2 = read_files("data/round2/")

own2=bind_rows(
    read_tsv("data/round2/stats/seqbuster_counts.tsv", col_names = F) %>% 
        group_by(X1) %>% 
        summarise(internalc=sum(X3)) %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "seqbuster"),
    read_tsv("data/round2/stats/razers3_counts.tsv", col_names = F) %>% 
        group_by(X1) %>% 
        summarise(internalc=sum(X3)) %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "razers3"),
    read_tsv("data/round2/own/srnabench.grouped") %>%  
        .[,c(1,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "srnabench"),
    read_tsv("data/miRNAs_expressed_all_samples_1556854694_Headers.csv") %>%   
        .[,c(1,5)] %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "mirdeep"),
    read_tsv("data/round2/own/excert/readCounts_miRNAmature_sense.txt") %>% 
        separate_rows(ReferenceID, sep="\\|") %>% 
        separate(ReferenceID, into = "mirna", extra = "drop", sep=":") %>% 
        .[,c(1,3)] %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "excert"),
    read_tsv("data/round2/own/manatee.tsv") %>% 
        .[,c(1,3)] %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "manatee"),
    read_tsv("data/round2/own/shortstack.txt") %>%  
        .[,c(2,4)] %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "shortstack")
) %>% group_by(mirna, tool) %>% 
    summarise(internalc=sum(internalc)) %>% 
    ungroup() %>% 
    full_join(simulated, by=c("mirna"="Transcript")) %>%
    select(-Count) %>% 
    mutate(internalc = ifelse(is.na(internalc), 0, internalc)) %>% 
    gather(method, counts, -mirna, -tool)

r2 = bind_rows(r2,own2) %>%
    full_join(simulated, by = c("mirna" = "Transcript") )

###########
#  round3
###########

simulated = read_tsv("input/round2/simulated_miRNA_BF_counts.tsv")
r3 = read_files("data/round3/")


own3=bind_rows(    
    read_tsv("data/round3/own/seqbuster.tsv",
             col_names = T) %>% 
        group_by(miRNA) %>% 
        summarise(internalc=sum(simulated_miRNA_BF)) %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "seqbuster"),
    read_tsv("data/round3/stats/razers3_counts.tsv", col_names = F) %>% 
        group_by(X1) %>% 
        summarise(internalc=sum(X3)) %>% 
        set_names("mirna", "internalc") %>% 
        filter(mirna %in% universe) %>% 
        mutate(tool = "razers3"),
    read_tsv("data/round3/own/srnabench.grouped") %>%  
        .[,c(1,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "srnabench"),
    read_tsv("data/miRNAs_expressed_all_samples_1556854694_Headers.csv") %>%   
        .[,c(1,6)] %>% 
        set_names("mirna", "internalc") %>% 
       filter(mirna %in% universe) %>% 
       mutate(tool = "mirdeep"),
    read_tsv("data/round3/own/excert/readCounts_miRNAmature_sense.txt") %>% 
       separate_rows(ReferenceID, sep="\\|") %>% 
       separate(ReferenceID, into = "mirna", extra = "drop", sep=":") %>% 
       .[,c(1,3)] %>% 
        set_names("mirna", "internalc") %>% 
       filter(mirna %in% universe) %>% 
       mutate(tool = "excert"),
    read_tsv("data/round3/own/manatee.tsv") %>% 
        .[,c(1,3)] %>% 
        set_names("mirna", "internalc") %>% 
       filter(mirna %in% universe) %>% 
       mutate(tool = "manatee"),
    read_tsv("data/round3/own/shortstack.txt") %>%  
        .[,c(2,4)] %>% 
        set_names("mirna", "internalc") %>% 
       filter(mirna %in% universe) %>% 
       mutate(tool = "shortstack")
) %>% group_by(mirna, tool) %>% 
    summarise(internalc=sum(internalc)) %>% 
    ungroup() %>% 
    full_join(simulated, by=c("mirna"="Transcript")) %>%
    select(-Count) %>% 
    mutate(internalc = ifelse(is.na(internalc), 0, internalc)) %>% 
    gather(method, counts, -mirna, -tool)

r3 = bind_rows(r3,own3) %>%
    full_join(simulated, by = c("mirna" = "Transcript") )

###########
#  round5
###########

simulated = read_tsv("input/round5/transcripts.tsv") %>% 
    filter(Transcipt %in% universe) %>% 
    rename(Transcript=Transcipt)
r5 = read_files("data/round5/")

own5=bind_rows(
    read_tsv("data/round5/own/seqbuster.mirna") %>% 
        group_by(mir) %>% 
        summarise(internalc=sum(freq)) %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "seqbuster"),
    read_tsv("data/round5/own/razers3.tsv") %>% 
        group_by(miRNA) %>% 
        mutate(freq=1) %>% 
        summarise(internalc=sum(freq)) %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "razers3"),
    read_tsv("data/round5/own/mirdeep.csv") %>%  
        .[,c(1,5)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "mirdeep"),
    read_tsv("data/round5/own/srnabench.grouped") %>%  
        .[,c(1,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "srnabench"),
    read_tsv("data/round5/own/excert/readCounts_miRNAmature_sense.txt") %>% 
        separate_rows(ReferenceID, sep="\\|") %>% 
        separate(ReferenceID, into = "mirna", extra = "drop", sep=":") %>% 
        .[,c(1,3)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "excert"),
    read_tsv("data/round5/own/manatee.tsv") %>% 
        .[,c(1,3)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "manatee"),
    read_tsv("data/round5/own/shortstack.txt") %>%  
        .[,c(2,4)] %>% 
        set_names("mirna", "internalc") %>% 
        right_join(simulated[,"Transcript"],
                   by=c("mirna"="Transcript")) %>% 
        mutate(tool = "shortstack")
) %>% group_by(mirna, tool) %>% 
    summarise(internalc=sum(internalc)) %>% 
    ungroup() %>% 
    full_join(simulated, by=c("mirna"="Transcript")) %>%
    select(-Count) %>% 
    mutate(internalc = ifelse(is.na(internalc), 0, internalc)) %>% 
    gather(method, counts, -mirna, -tool)

r5 = bind_rows(r5,own5) %>%
    full_join(simulated, by = c("mirna" = "Transcript") )


bind_rows(
    r1 %>% mutate(round="S1"),
    r2 %>% mutate(round="S2"),
    r3 %>% mutate(round="S3"),
    r5 %>% mutate(round="S5")
) %>% 
 #    r1 %>% mutate(round="S1") %>% 
    mutate(
        fc = as.numeric(counts)/Count,
    ) %>%
    filter(
           tool!="bowtie1_default",
           tool!="gsnap_m10") %>%
    #filter(!(grepl("mirtop", method) & tool!="razers3")) %>%
    mutate(method = gsub("_fc", "", method)) %>% 
    distinct() %>% 
    mutate(fc = cut(fc,
                    c(-1,0.0001,0.5,0.7,0.9,1.1,1.3,2,1000000), 
                    labs),
           fc = as.character(fc),
           fc = ifelse(is.na(fc),"extra", fc),
           fc = factor(fc, levels=rev(c(labs, "extra")))) %>%  
    
    ggplot(aes(paste(tool, method), fill=fc)) +
    geom_bar() +
    scale_fill_manual(values=c( "black",
                                "orange3", "yellow3",
                                "#addd8e", "deepskyblue3",
                                "#addd8e","#edf8b1",
                                "#feb24c", "black")) +
    ggthemes::theme_base() +
    coord_flip() +
    facet_wrap(~round, nrow=1, scales = "free_x") + 
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
    ggsave("figures-local/expression.pdf", width = 12, height = 12)


