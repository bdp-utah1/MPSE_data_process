#!/usr/bin/env R

library(tidyverse)
library(ontologyIndex)

setwd("~/Documents/phd/yandell/software/MPSE_data_process")


clix_clin <- read_tsv("analysis/set_comparison/clix_vs_clinphen.txt", comment="#")
clix_ctak <- read_tsv("analysis/set_comparison/clix_vs_ctakes.txt", comment="#")
clix_meta <- read_tsv("analysis/set_comparison/clix_vs_metamap.txt", comment="#")
clin_ctak <- read_tsv("analysis/set_comparison/clinphen_vs_ctakes.txt", comment="#")
clin_meta <- read_tsv("analysis/set_comparison/clinphen_vs_metamap.txt", comment="#")
ctak_meta <- read_tsv("analysis/set_comparison/ctakes_vs_metamap.txt", comment="#")


set_comp_wide <- bind_rows(list("clix_clin"=clix_clin,
                                 "clix_ctak"=clix_ctak,
                                 "clix_meta"=clix_meta,
                                 "clin_ctak"=clin_ctak,
                                 "clin_meta"=clin_meta,
                                 "ctak_meta"=ctak_meta),
                            .id="pair") %>% 
  separate_wider_delim(pair, delim="_", names=c("set1","set2"), cols_remove=FALSE) %>% 
  mutate(set1 = factor(set1, levels=c("clix","clin","ctak"), labels=c("CLiX","ClinPhen","cTAKES")),
         set2 = factor(set2, levels=c("clin","ctak","meta"), labels=c("ClinPhen","cTAKES","MetaMap"))) %>% 
  rename(real_jac = "jac",
         real_jac_p = "jac_p",
         real_onto_sim = "onto_sim")

set_comp_long <- set_comp_wide %>% 
  pivot_longer(cols=real_jac:rand_onto_sim,
               names_to="source",
               values_to="value") %>% 
  separate_wider_delim(source, delim="_", names=c("sampling","index"), too_many="merge")


ggplot(data=set_comparison, aes(x=cnt1, y=cnt2)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, linetype=2, alpha=0.5) +
  facet_grid(vars(set2), vars(set1)) +
  scale_x_continuous(name="HPO Term Count",
                     limits=c(14,300)) +
  scale_y_continuous(name="HPO Term Count",
                     limits=c(13,300)) +
  theme_linedraw()


set_comp_long %>% 
  filter(index=="jac_p") %>% 
  ggplot(aes(x=value, fill=sampling)) +
  geom_histogram(binwidth=0.02) + 
  facet_grid(vars(set2), vars(set1)) +
  scale_x_continuous(name="Jaccard Index (+parents)",
                     breaks=seq(0,1,0.2),
                     labels=seq(0,1,0.2),
                     limits=c(0,1)) +
  theme_linedraw()


set_comp_long %>% 
  filter(index=="onto_sim") %>% 
  ggplot(aes(x=value, fill=sampling)) +
  geom_histogram(binwidth=0.02) + 
  facet_grid(vars(set2), vars(set1)) +
  scale_x_continuous(name="Semantic Similarity",
                     breaks=seq(0,1,0.2),
                     labels=seq(0,1,0.2),
                     limits=c(0,1)) +
  theme_linedraw()




clix_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_clinithink_stats.txt")
clin_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_clinphen_stats.txt")
ctak_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_ctakes_stats.txt")
meta_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_metamap_stats.txt")


tool_stats <- bind_rows(list("clix"=clix_stats,
               "clin"=clin_stats,
               "ctak"=ctak_stats,
               "meta"=meta_stats),
          .id="tool") %>% 
  mutate(tool = factor(tool, levels=c("clix","clin","ctak","meta"), labels=c("CLiX","ClinPhen","cTAKES","MetaMap")))


tool_stats %>% 
  group_by(tool) %>% 
  summarise(term_cnt = mean(term_count),
            info_cont = mean(mean_info_cont),
            total_info_cont = mean(total_info_cont))


ggplot(tool_stats, aes(x=tool, y=term_count, fill=tool)) +
  geom_violin(alpha=0.5) + 
  geom_dotplot(binaxis="y",
               stackdir="center",
               dotsize=0.4) +
  scale_x_discrete(name="NLP Tool") +
  scale_y_continuous(name="HPO Term Count",
                     breaks=seq(0,300,50),
                     labels=seq(0,300,50),
                     limits=c(0,300)) +
  theme(legend.position="none")

ggplot(tool_stats, aes(x=tool, y=mean_info_cont, fill=tool)) +
  geom_violin(alpha=0.5) + 
  geom_dotplot(binaxis="y",
               stackdir="center",
               dotsize=0.4) +
  scale_x_discrete(name="NLP Tool") +
  scale_y_continuous(name="Mean Information Content",
                     breaks=seq(3.5,6,0.5),
                     labels=seq(3.5,6,0.5),
                     limits=c(3.5,6)) +
  theme(legend.position="none")

ggplot(tool_stats, aes(x=tool, y=total_info_cont, fill=tool)) +
  geom_violin(alpha=0.5) + 
  geom_dotplot(binaxis="y",
               stackdir="center",
               dotsize=0.4) +
  scale_x_discrete(name="NLP Tool") +
  scale_y_continuous(name="Total Information Content",
                     breaks=seq(0,1200,200),
                     labels=seq(0,1200,200),
                     limits=c(0,1200)) +
  theme(legend.position="none")





setwd("~/Documents/phd/yandell/software/MPSE")
training <- read_tsv("analysis/nlp_training_sets/training_preds_combined.tsv")
testing <- read_tsv("analysis/nlp_testing_sets/testing_preds_combined.tsv")


training %>% 
  select(cohort, codes_source, scr) %>% 
  mutate(cohort = factor(cohort, 
                         levels=c("utah","neoseq"),
                         labels=c("UofU Controls", "NeoSeq Cases")),
         codes_source = factor(codes_source, 
                               levels=c("edw_clix","edw_clinphen","edw_ctakes","edw_metamap"), 
                               labels=c("CLiX","ClinPhen","cTAKES","MetaMap"))) %>%
  ggplot(aes(x=scr, color=cohort)) + 
  geom_density() + 
  facet_wrap(vars(codes_source), nrow=4, ncol=1, scales="free_y") + 
  scale_x_continuous(name="MPSE Score",
                     breaks=seq(-100,300,50),
                     labels=seq(-100,300,50)) + 
  scale_color_manual(name="Cohort", 
                     values=c("UofU Controls"="#00BFC4", 
                              "NeoSeq Cases"="#F8766D")) +
  theme_bw()


testing %>% 
  select(cohort, codes_source, scr) %>% 
  mutate(cohort = factor(cohort, 
                         levels=c("utah","neoseq"),
                         labels=c("UofU Controls", "NeoSeq Cases")),
         codes_source = factor(codes_source, 
                               levels=c("edw_clix","edw_clinphen","edw_ctakes","edw_metamap"), 
                               labels=c("CLiX","ClinPhen","cTAKES","MetaMap"))) %>% 
  ggplot(aes(x=scr, color=cohort)) + 
  geom_density() + 
  facet_wrap(vars(codes_source), nrow=4, ncol=1, scales="free_y") + 
  scale_x_continuous(name="MPSE Score",
                     breaks=seq(-75,150,25),
                     labels=seq(-75,150,25)) + 
  scale_color_manual(name="Cohort", 
                     values=c("UofU Controls"="#00BFC4", 
                              "NeoSeq Cases"="#F8766D")) +
  theme_bw()

