#!/usr/bin/env R

library(tidyverse)
library(ontologyIndex)

setwd("~/Documents/phd/yandell/software/MPSE")


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
