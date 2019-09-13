library(tidyverse)

###########
# counts
###########
labs= c("missed","<50%x", "<30%x", "<10%x","+/-10%x", ">10%x", ">30%x", ">50%x")

read_files = function(folder, simulated, razers=T, bwa=T, seqbuster=T,
                      srnabench=T, manatee=T, shortstack=T, mirdeep=T){
    universe = as.character(simulated$Transcript)
    
    featurecounts = list.files(file.path(folder, "features"),
                               ".txt$",
                               full.names = T) %>%
        lapply(., function(fn){
            name = gsub("\\.txt","",basename(fn))
            read_tsv(fn, skip = 1) %>%
                .[,c(1,7)] %>%
                set_names("mirna", "fcu") %>%
                full_join(simulated, by=c("mirna"="Transcript")) %>%
                rename(real = Count) %>%
                gather(method, counts, -mirna, -real) %>% 
                mutate(counts = ifelse(is.na(counts), 0, counts)) %>% 
                mutate(tool=name)
        }) %>% bind_rows()
    
    features_multiON=list.files(file.path(folder, "features_multiON"),
                                ".txt$",
                                full.names = T) %>%
        lapply(., function(fn){
            name = gsub("\\.txt","",basename(fn))
            read_tsv(fn, skip = 1) %>%
                .[,c(1,7)] %>%
                set_names("mirna", "fcm") %>%
                full_join(simulated, by=c("mirna"="Transcript")) %>%
                rename(real = Count) %>%
                gather(method, counts, -mirna, -real) %>% 
                mutate(counts = ifelse(is.na(counts), 0, counts)) %>% 
                mutate(tool=name)
        }) %>% bind_rows()
    
    # browser()
    
     df = bind_rows(featurecounts, features_multiON) %>% 
         filter(!(tool %in% c("razers3-pre", "ShortStack", "manatee", "bwa-pre")))
    universe = unique(df$mirna)
    # browser()
    # df = df %>% gather(method, counts, -mirna, -tool) %>% 
    #     filter(counts>0)
    mi_counts = . %>% 
        group_by(mirna) %>% 
        summarise(custom=sum(custom)) %>% 
        ungroup() %>% 
        full_join(simulated, by=c("mirna"="Transcript")) %>%
        rename(real = Count) %>% 
        gather(method, counts, -mirna, -real) %>% 
        mutate(counts = ifelse(is.na(counts), 0, counts))
   # browser()
    if (razers==T){
        df = bind_rows(df, 
                       read_tsv(file.path(folder, "stats/razers3-pre_counts.tsv"),
                                col_names = F) %>% 
                           group_by(X1) %>% 
                           summarise(internalc=sum(X3)) %>% 
                           set_names("mirna", "custom") %>% 
                           mi_counts() %>% 
                           mutate(tool = "razers3-pre")
        )
    }
    if (bwa==T){
        df = bind_rows(df, 
                       read_tsv(file.path(folder, "stats/bwa-pre_counts.tsv"),
                                col_names = F) %>% 
                           group_by(X1) %>% 
                           summarise(internalc=sum(X3)) %>% 
                           set_names("mirna", "custom") %>% 
                           mi_counts() %>% 
                           mutate(tool = "bwa-pre")
        )
    }
    if (seqbuster==T){
        df = bind_rows(df, 
                       read_tsv(file.path(folder, "quant/seqbuster.mirna"),
                                col_names = T) %>% 
                           group_by(mir) %>% 
                           summarise(internalc=sum(freq)) %>% 
                           set_names("mirna", "custom") %>% 
                           mi_counts() %>% 
                           mutate(tool = "seqbuster")
        )
    }
    if (srnabench==T){
        df = bind_rows(df, 
                       read_tsv(file.path(folder, "quant/srnabench.txt")) %>%  
                           .[,c(1,4)] %>% 
                           set_names("mirna", "custom") %>% 
                           mi_counts() %>% 
                           mutate(tool = "srnabench")
        )
    }
    if (manatee==T){
        df = bind_rows(df, 
                       read_tsv(file.path(folder, "quant/Manatee_counts.tsv")) %>% 
                           .[,c(1,3)] %>% 
                           set_names("mirna", "custom") %>% 
                           mi_counts() %>% 
                           filter(mirna %in% universe) %>% 
                           mutate(tool = "manatee") 
        )
    }
    if (shortstack==T){
        df = bind_rows(df, 
                       read_tsv(file.path(folder, "quant/shortstack.txt")) %>% 
                           .[,c(2,4)] %>% 
                           set_names("mirna", "custom") %>% 
                           mi_counts() %>% 
                           filter(mirna %in% universe) %>% 
                           mutate(tool = "shortstack") 
        )
    }
    if (mirdeep==T){
        df = bind_rows(df, 
                       read_tsv(file.path(folder, "quant/miRNAs_expressed_all_samples_miRDeep.csv")) %>% 
                           .[,c(1,4)] %>% 
                           set_names("mirna", "custom") %>% 
                           mi_counts() %>% 
                           mutate(tool = "mirdeep") 
        )
    }
        # read_tsv("data/round1/own/excert/readCounts_miRNAmature_sense.txt") %>% 
        #     separate_rows(ReferenceID, sep="\\|") %>% 
        #     separate(ReferenceID, into = "mirna", extra = "drop", sep=":") %>% 
        #     .[,c(1,3)] %>% 
        #     set_names("mirna", "internalc") %>% 
        #     filter(mirna %in% universe) %>% 
        #     mutate(tool = "excert"),
        # read_tsv("data/round1/own/shortstack.txt") %>%  
        #     .[,c(2,4)] %>% 
        #     set_names("mirna", "internalc") %>% 
        #     filter(mirna %in% universe) %>% 
        #     mutate(tool = "shortstack")
}

###########
#  round0
###########
r0 = read_files("data/output.S0", # output files
                read_tsv("input/simulated_data/S0/S0_counts.tsv")) # simulated

###########
#  round1
###########
r1 =read_files("data/output.S1", # output files
               read_tsv("input/simulated_data/S1/S1_counts.tsv")) # simulated

###########
#  round2
###########
r2 =read_files("data/output.S2", # output files
               read_tsv("input/simulated_data/S2/S2_counts.tsv")) # simulated

###########
#  round3
###########
r3 =read_files("data/output.S3", # output files
               read_tsv("input/simulated_data/S3/S3_counts.tsv") %>% 
                   filter(grepl("hsa-", Transcript))) # simulated

###########
#  round4
###########
r4 =read_files("data/output.S4", # output files
               read_tsv("input/simulated_data/S4/S4_counts.tsv") %>% 
                   filter(grepl("hsa-", Transcript))) # simulated


###########
#  round5
###########
r5 =read_files("data/output.S5", # output files
               read_tsv("input/simulated_data/S5/S5_counts.tsv") %>% 
                   rename(Transcript=Transcipt) %>% 
                   filter(grepl("hsa-", Transcript))) # simulated

#############################################################################
allcounst = unique(paste(r0$tool, r0$method))
levels = c(allcounst[grepl("custom", allcounst)], allcounst[!grepl("custom", allcounst)])
# cols = c( "black",
#           "orange3", "yellow3",
#           "#addd8e", "deepskyblue3",
#           "#addd8e","#edf8b1",
#           "#feb24c", "grey")
cols = c( "black",
          "#de2d26", "#fc9272",
          "#fee0d2", "#99d8c9",
          "#deebf7","#9ecae1",
          "#3182bd", "grey70")

bind_rows(
    r0 %>% mutate(round="S0"),
    r1 %>% mutate(round="S1"),
    r2 %>% mutate(round="S2"),
    r3 %>% mutate(round="S3"),
    r4 %>% mutate(round="S4"),
    r5 %>% mutate(round="S5")
) %>% 
    filter(!(is.na(real) & counts==0)) %>% 
 #    r1 %>% mutate(round="S1") %>% 
    mutate(
        fc = as.numeric(counts)/real,
    ) %>%
    mutate(method = gsub("_fc", "", method)) %>% 
    distinct() %>% 
    mutate(fc = cut(fc,
                    c(-1,0.0001,0.5,0.7,0.9,1.1,1.3,2,1000000), 
                    labs),
           fc = as.character(fc),
           fc = ifelse(is.na(fc),"extra", fc),
           fc = factor(fc, levels=rev(c(labs, "extra")))) %>%  
    mutate(tool_method = factor(paste(tool, method), level=levels)) %>% 
    write_csv("results/expression_accuracy.csv")

read_csv("results/expression_accuracy.csv") %>% 
    ggplot(aes(tool_method, fill=fc)) +
    geom_bar() +
    scale_fill_manual(values=cols) +
    ggthemes::theme_base() +
    coord_flip() +
    facet_wrap(~round, nrow=1, scales = "free_x") + 
    theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)) +
    ggsave("figures-local/expression.pdf", width = 12, height = 12)


