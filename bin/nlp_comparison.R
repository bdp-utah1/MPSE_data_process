#!/usr/bin/env R

library(tidyverse)
library(ontologyIndex)
library(pROC)
library(plotROC)

setwd("~/Documents/phd/yandell/software/MPSE_data_process")

hpo <- get_ontology("../MPSE/docs/hp.obo")




orig_rady_train <- read_tsv("tests/test_null/test_preds.tsv") %>% 
  select(seq_status, pos_proba, scr)
null_rady_train <- read_tsv("tests/test_null/weighted_null_preds.tsv") %>% 
  select(seq_status, pos_proba, scr)
df <- bind_rows(list(orig=orig_rady_train, null=null_rady_train), .id="cohort") %>% 
  mutate(cohort=factor(cohort, levels=c("orig","null")))

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Plot ROC curves for original HPO terms vs randomized HPO terms using Rady training dataset
ggplot(df, aes(d=seq_status, m=scr, color=cohort)) + 
  geom_roc() + 
  style_roc() + 
  geom_segment(aes(x=0, y=0, xend=1, yend=1), color="black", linetype=2, linewidth=0.4) + 
  scale_color_discrete(guide = guide_legend(title="Term Type"),
                       labels = c("Original Terms",
                                  "Random Terms"))

auc(roc(orig_rady_train$seq_status, orig_rady_train$pos_proba))
auc(roc(null_rady_train$seq_status, null_rady_train$pos_proba))
## -------------------------------------------------------------------------- ##




clix_clin <- read_tsv("analysis/set_comparison/clix_vs_clinphen.txt", comment="#")
clix_ctak <- read_tsv("analysis/set_comparison/clix_vs_ctakes.txt", comment="#")
clix_meta <- read_tsv("analysis/set_comparison/clix_vs_metamap.txt", comment="#")
clix_medl <- read_tsv("analysis/set_comparison/clix_vs_medlee.txt", comment="#")
clin_ctak <- read_tsv("analysis/set_comparison/clinphen_vs_ctakes.txt", comment="#")
clin_meta <- read_tsv("analysis/set_comparison/clinphen_vs_metamap.txt", comment="#")
clin_medl <- read_tsv("analysis/set_comparison/clinphen_vs_medlee.txt", comment="#")
ctak_meta <- read_tsv("analysis/set_comparison/ctakes_vs_metamap.txt", comment="#")
ctak_medl <- read_tsv("analysis/set_comparison/ctakes_vs_medlee.txt", comment="#")
meta_medl <- read_tsv("analysis/set_comparison/metamap_vs_medlee.txt", comment="#")

manu_clix <- read_tsv("analysis/set_comparison/manual_vs_clix.txt", comment="#")
manu_clin <- read_tsv("analysis/set_comparison/manual_vs_clinphen.txt", comment="#")
manu_ctak <- read_tsv("analysis/set_comparison/manual_vs_ctakes.txt", comment="#")
manu_meta <- read_tsv("analysis/set_comparison/manual_vs_metamap.txt", comment="#")
manu_medl <- read_tsv("analysis/set_comparison/manual_vs_medlee.txt", comment="#")

set_comp_wide <- bind_rows(list("clix_clin"=clix_clin,
                                "clix_ctak"=clix_ctak,
                                "clix_meta"=clix_meta,
                                "clix_medl"=clix_medl,
                                "clin_ctak"=clin_ctak,
                                "clin_meta"=clin_meta,
                                "clin_medl"=clin_medl,
                                "ctak_meta"=ctak_meta,
                                "ctak_medl"=ctak_medl,
                                "meta_medl"=meta_medl,
                                "manu_clix"=manu_clix,
                                "manu_clin"=manu_clin,
                                "manu_ctak"=manu_ctak,
                                "manu_meta"=manu_meta,
                                "manu_medl"=manu_medl),
                            .id="pair") %>% 
  separate_wider_delim(pair, delim="_", names=c("set1","set2"), cols_remove=FALSE) %>% 
  mutate(set1 = factor(set1, levels=c("manu","clix","clin","ctak","meta"), labels=c("Manual","CLiX","ClinPhen","cTAKES","MetaMap")),
         set2 = factor(set2, levels=c("clix","clin","ctak","meta","medl"), labels=c("CLiX","ClinPhen","cTAKES","MetaMap","MedLEE"))) %>% 
  rename(real_jac = "jac",
         real_jac_p = "jac_p",
         real_onto_sim = "onto_sim")

set_comp_long <- set_comp_wide %>% 
  pivot_longer(cols=real_jac:rand_onto_sim,
               names_to="source",
               values_to="value") %>% 
  separate_wider_delim(source, delim="_", names=c("sampling","index"), too_many="merge")

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# summary statistics of Jaccard Similarity Index for randomly sampled HPO sets
# with set counts equal to the real data
set_comp_wide %>% 
  group_by(pair) %>% 
  summarise(min_cnt = min(rand_jac_p),
            q1_cnt = quantile(rand_jac_p, probs=c(0.25)),
            median_cnt = median(rand_jac_p),
            q3_cnt = quantile(rand_jac_p, probs=c(0.75)),
            max_cnt = max(rand_jac_p),
            mean_cnt = mean(rand_jac_p),
            n = n())

# summary statistics of Jaccard Similarity Index between different NLP sets
# '_p" indicates that parent terms have been supplied
set_comp_wide %>% 
  group_by(pair) %>% 
  summarise(min_cnt = min(real_jac_p),
            q1_cnt = quantile(real_jac_p, probs=c(0.25)),
            median_cnt = median(real_jac_p),
            q3_cnt = quantile(real_jac_p, probs=c(0.75)),
            max_cnt = max(real_jac_p),
            mean_cnt = mean(real_jac_p),
            n = n()) %>% 
  mutate(pair = factor(pair, 
                       levels=c("clix_clin","clix_ctak","clix_meta","clix_medl",
                                "clin_ctak","clin_meta","clin_medl","ctak_meta","ctak_medl",
                                "meta_medl","manu_clix","manu_clin","manu_ctak","manu_meta","manu_medl"), 
                       labels=c("CLiX-ClinPhen","CLiX-cTAKES","CLiX-MetaMap","CLiX-MedLEE",
                                "ClinPhen-cTAKES","ClinPhen-MetaMap","ClinPhen-MedLEE","cTAKES-MetaMap","cTAKES-MedLEE",
                                "MetaMap-MedLEE","Manual-CLiX","Manual-ClinPhen","Manual-cTAKES","Manual-MetaMap","Manual-MedLEE"))) %>% 
  arrange(pair)
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# T-tests comparing Jaccard Similarity Index between randomly sampled HPO sets 
# and real data
set_comp_long %>% 
  filter(pair=="clix_clin",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clix_ctak",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clix_meta",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clix_medl",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clin_ctak",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clin_meta",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clin_medl",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="ctak_meta",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="ctak_medl",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="meta_medl",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)

# T-tests comparing Jaccard Similarity Index between randomly sampled HPO sets
# and manual HPO sets
set_comp_long %>% 
  filter(pair=="manu_clix",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_clin",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_ctak",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_meta",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_medl",
         index=="jac_p") %>% 
  t.test(value ~ sampling, data=.)
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Compute z-scores for NLP set Jaccard Similarities using randomly sampled HPO
# sets to estimate the population mean and sd
jac_pop_params <- set_comp_long %>% 
  filter(sampling=="rand",
         index=="jac_p") %>% 
  group_by(pair) %>% 
  summarise(pop_mean=mean(value),
            pop_sd=sd(value))

# z-score summary statistics
set_comp_long %>% 
  filter(sampling=="real",
         index=="jac_p") %>% 
  left_join(jac_pop_params, by="pair") %>% 
  mutate(z_scr=(value - pop_mean)/pop_sd) %>% 
  group_by(pair) %>% 
  summarise(min_z = min(z_scr),
            q1_z = quantile(z_scr, probs=c(0.25)),
            median_z = median(z_scr),
            q3_z = quantile(z_scr, probs=c(0.75)),
            max_z = max(z_scr),
            mean_z = mean(z_scr),
            n = n())
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# summary statistics of Semantic Similarity Index for randomly sampled HPO sets
# with set counts equal to the real data
set_comp_wide %>% 
  group_by(pair) %>% 
  summarise(min_cnt = min(rand_onto_sim),
            q1_cnt = quantile(rand_onto_sim, probs=c(0.25)),
            median_cnt = median(rand_onto_sim),
            q3_cnt = quantile(rand_onto_sim, probs=c(0.75)),
            max_cnt = max(rand_onto_sim),
            mean_cnt = mean(rand_onto_sim),
            n = n())

# summary statistics of Semantic Similarity Index between different NLP sets
set_comp_wide %>% 
  group_by(pair) %>% 
  summarise(min_cnt = min(real_onto_sim),
            q1_cnt = quantile(real_onto_sim, probs=c(0.25)),
            median_cnt = median(real_onto_sim),
            q3_cnt = quantile(real_onto_sim, probs=c(0.75)),
            max_cnt = max(real_onto_sim),
            mean_cnt = mean(real_onto_sim),
            n = n())
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# T-tests comparing Semantic Similarity Index between randomly sampled HPO sets 
# and real data
set_comp_long %>% 
  filter(pair=="clix_clin",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clix_ctak",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clix_meta",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clix_medl",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clin_ctak",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clin_meta",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="clin_medl",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="ctak_meta",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="ctak_medl",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="meta_medl",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)

# T-tests comparing Semantic Similarity Index between randomly sampled HPO sets
# and manual HPO sets
set_comp_long %>% 
  filter(pair=="manu_clix",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_clin",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_ctak",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_meta",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
set_comp_long %>% 
  filter(pair=="manu_medl",
         index=="onto_sim") %>% 
  t.test(value ~ sampling, data=.)
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Compute z-scores for NLP set Jaccard Similarities using randomly sampled HPO
# sets to estimate the population mean and sd
sem_pop_params <- set_comp_long %>% 
  filter(sampling=="rand",
         index=="onto_sim") %>% 
  group_by(pair) %>% 
  summarise(pop_mean=mean(value),
            pop_sd=sd(value))

# z-score summary statistics
set_comp_long %>% 
  filter(sampling=="real",
         index=="onto_sim") %>% 
  left_join(sem_pop_params, by="pair") %>% 
  mutate(z_scr=(value - pop_mean)/pop_sd) %>% 
  group_by(pair) %>% 
  summarise(min_z = min(z_scr),
            q1_z = quantile(z_scr, probs=c(0.25)),
            median_z = median(z_scr),
            q3_z = quantile(z_scr, probs=c(0.75)),
            max_z = max(z_scr),
            mean_z = mean(z_scr),
            n = n())
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# matrix of scatter plots comparing HPO term counts between Manual/NLP sources
ggplot(data=set_comp_wide, aes(x=cnt1, y=cnt2)) +
  geom_point(shape=1) +
  geom_abline(slope=1, intercept=0, linetype=2, alpha=0.5) +
  facet_grid(vars(set2), vars(set1)) +
  scale_x_continuous(name="HPO Term Count",
                     limits=c(3,300)) +
  scale_y_continuous(name="HPO Term Count",
                     limits=c(10,300)) +
  theme_linedraw() +
  theme(axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        strip.text=element_text(size=14))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# matrix of histograms showing Jaccard similarity scores between 
# manual/NLP term sets
set_comp_long %>% 
  filter(index=="jac",
         sampling=="real") %>% 
  ggplot(aes(x=value, fill=sampling)) +
  geom_histogram(binwidth=0.02, fill="#00BFC4") + 
  facet_grid(vars(set2), vars(set1)) +
  scale_x_continuous(name="Jaccard Index",
                     breaks=seq(0,1,0.2),
                     labels=seq(0,1,0.2),
                     limits=c(0,1)) +
  scale_y_continuous(name="Count") +
  theme_linedraw() +
  theme(legend.position="none",
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        strip.text=element_text(size=14))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# matrix of histograms showing Jaccard similarity (+parents) scores between 
# manual/NLP term sets
set_comp_long %>% 
  filter(index=="jac_p") %>% 
  ggplot(aes(x=value, fill=sampling)) +
  geom_histogram(binwidth=0.02) + 
  facet_grid(vars(set2), vars(set1)) +
  scale_x_continuous(name="Jaccard Index (+parents)",
                     breaks=seq(0,1,0.2),
                     labels=seq(0,1,0.2),
                     limits=c(0,1)) +
  scale_y_continuous(name="Count") +
  theme_linedraw() +
  theme(legend.position="none",
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        strip.text=element_text(size=14))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# matrix of histograms showing Semantic similarity scores between manual/NLP 
# term sets
set_comp_long %>% 
  filter(index=="onto_sim") %>% 
  ggplot(aes(x=value, fill=sampling)) +
  geom_histogram(binwidth=0.02) + 
  facet_grid(vars(set2), vars(set1)) +
  scale_x_continuous(name="Semantic Similarity",
                     breaks=seq(0,1,0.2),
                     labels=seq(0,1,0.2),
                     limits=c(0,1)) +
  scale_y_continuous(name="Count") +
  theme_linedraw() + 
  theme(legend.position="none",
        axis.title=element_text(size=16),
        axis.text=element_text(size=12),
        strip.text=element_text(size=14))
## -------------------------------------------------------------------------- ##




rady_mod_coef <- read_tsv("analysis/nlp_model_coefficients/rady/model_coefficients.tsv") %>%
  mutate(ref="rady", coef=log_prob_1 - log_prob_0) %>%
  arrange(coef)
clix_mod_coef <- read_tsv("analysis/nlp_model_coefficients/clix/model_coefficients.tsv") %>%
  mutate(ref="clix", coef=log_prob_1 - log_prob_0) %>%
  arrange(coef)
clin_mod_coef <- read_tsv("analysis/nlp_model_coefficients/clin/model_coefficients.tsv") %>%
  mutate(ref="clin", coef=log_prob_1 - log_prob_0) %>%
  arrange(coef)
ctak_mod_coef <- read_tsv("analysis/nlp_model_coefficients/ctak/model_coefficients.tsv") %>%
  mutate(ref="ctak", coef=log_prob_1 - log_prob_0) %>%
  arrange(coef)
meta_mod_coef <- read_tsv("analysis/nlp_model_coefficients/meta/model_coefficients.tsv") %>%
  mutate(ref="meta", coef=log_prob_1 - log_prob_0) %>%
  arrange(coef)

mod_coef <- bind_rows(rady_mod_coef,
                      clix_mod_coef,
                      clin_mod_coef,
                      ctak_mod_coef,
                      meta_mod_coef) %>% 
  mutate(term_name=map_chr(term, get_term_property, ontology=hpo, property_name="name"))

clix_set <- read_tsv("data/full_NLP_datasets/clinithink.tsv") %>% 
  select(pid,cohort,codes) %>% 
  mutate(source="clix") %>% 
  separate_longer_delim(codes, delim=";") %>% 
  distinct()
clin_set <- read_tsv("data/full_NLP_datasets/clinphen.tsv") %>% 
  select(pid,cohort,codes) %>% 
  mutate(source="clin") %>% 
  separate_longer_delim(codes, delim=";") %>% 
  distinct()
ctak_set <- read_tsv("data/full_NLP_datasets/ctakes.tsv") %>% 
  select(pid,cohort,codes) %>% 
  mutate(source="ctak") %>% 
  separate_longer_delim(codes, delim=";") %>% 
  distinct()
meta_set <- read_tsv("data/full_NLP_datasets/metamaplite.tsv") %>% 
  select(pid,cohort,codes) %>% 
  mutate(source="meta") %>% 
  separate_longer_delim(codes, delim=";") %>% 
  distinct()
nlp_sets <- bind_rows(clix_set,
                      clin_set,
                      ctak_set,
                      meta_set) %>% 
  separate_rows(codes, sep=";") %>% 
  rename(term=codes) %>% 
  left_join(rady_mod_coef, by="term") %>% 
  rename(source=source.x) %>% 
  select(pid, cohort, term, source, coef)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Quantify model coefficient differences between NLP term sets
clix_clin_join <- full_join(clix_set, clin_set, by=c("pid","cohort","codes")) %>% 
  mutate(grouping=if_else(is.na(source.x), "y", if_else(is.na(source.y), "x", "both"))) %>% 
  rename(term=codes) %>% 
  inner_join(rady_mod_coef, by="term") %>% 
  select(pid, cohort, term, grouping, ref, coef) %>% 
  mutate(set_x="clix", set_y="clin")
clix_ctak_join <- full_join(clix_set, ctak_set, by=c("pid","cohort","codes")) %>% 
  mutate(grouping=if_else(is.na(source.x), "y", if_else(is.na(source.y), "x", "both"))) %>% 
  rename(term=codes) %>% 
  inner_join(rady_mod_coef, by="term") %>% 
  select(pid, cohort, term, grouping, ref, coef) %>% 
  mutate(set_x="clix", set_y="ctak")
clix_meta_join <- full_join(clix_set, meta_set, by=c("pid","cohort","codes")) %>% 
  mutate(grouping=if_else(is.na(source.x), "y", if_else(is.na(source.y), "x", "both"))) %>% 
  rename(term=codes) %>% 
  inner_join(rady_mod_coef, by="term") %>% 
  select(pid, cohort, term, grouping, ref, coef) %>% 
  mutate(set_x="clix", set_y="meta")
clin_ctak_join <- full_join(clin_set, ctak_set, by=c("pid","cohort","codes")) %>% 
  mutate(grouping=if_else(is.na(source.x), "y", if_else(is.na(source.y), "x", "both"))) %>% 
  rename(term=codes) %>% 
  inner_join(rady_mod_coef, by="term") %>% 
  select(pid, cohort, term, grouping, ref, coef) %>% 
  mutate(set_x="clin", set_y="ctak")
clin_meta_join <- full_join(clin_set, meta_set, by=c("pid","cohort","codes")) %>% 
  mutate(grouping=if_else(is.na(source.x), "y", if_else(is.na(source.y), "x", "both"))) %>% 
  rename(term=codes) %>% 
  inner_join(rady_mod_coef, by="term") %>% 
  select(pid, cohort, term, grouping, ref, coef) %>% 
  mutate(set_x="clin", set_y="meta")
ctak_meta_join <- full_join(ctak_set, meta_set, by=c("pid","cohort","codes")) %>% 
  mutate(grouping=if_else(is.na(source.x), "y", if_else(is.na(source.y), "x", "both"))) %>% 
  rename(term=codes) %>% 
  inner_join(rady_mod_coef, by="term") %>% 
  select(pid, cohort, term, grouping, ref, coef) %>% 
  mutate(set_x="ctak", set_y="meta")

nlp_coefs_long <- bind_rows(clix_clin_join,
          clix_ctak_join,
          clix_meta_join,
          clin_ctak_join,
          clin_meta_join,
          ctak_meta_join) %>% 
  group_by(pid, cohort, set_x, set_y, grouping, .drop=FALSE) %>% 
  summarise(avg_mc=mean(coef),
            tot_mc=sum(coef),
            term_cnt=n()) %>% 
  ungroup() %>% 
  complete(pid, grouping, nesting(set_x, set_y), fill=list(cohort="utah", avg_mc=0, tot_mc=0, term_cnt=0))

nlp_coefs_wide <- nlp_coefs_long %>% 
  ungroup() %>% 
  group_by(cohort, set_x, set_y, grouping) %>% 
  summarise(avg_mc2=mean(avg_mc),
            tot_mc2=mean(tot_mc),
            avg_term_cnt=mean(term_cnt),
            n=n())

nlp_coefs_wide %>% 
  select(-tot_mc2) %>% 
  pivot_wider(names_from=grouping, values_from=c(avg_mc2, avg_term_cnt))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Pairwise matrix of barcharts showing NLP set unique content comparisons
cnt_order <- clix_clin_join %>% 
  filter(cohort=="neoseq", grouping=="x") %>% 
  group_by(pid) %>% 
  summarise(cnt=n()) %>% 
  arrange(cnt)

bind_rows(clix_clin_join,
          clix_ctak_join,
          clix_meta_join,
          clin_ctak_join,
          clin_meta_join,
          ctak_meta_join) %>% 
  filter(cohort=="neoseq") %>% 
  group_by(pid, set_x, set_y, grouping) %>% 
  summarise(avg_coef=mean(coef),
            cnt=n()) %>% 
  filter(grouping != "both") %>% 
  mutate(cnt_polar=if_else(grouping=="x", -cnt, cnt),
         pid=factor(pid, levels=cnt_order$pid),
         set_x = factor(set_x, levels=c("clix","clin","ctak"), labels=c("CLiX","ClinPhen","cTAKES")),
         set_y = factor(set_y, levels=c("clin","ctak","meta"), labels=c("ClinPhen","cTAKES","MetaMap"))) %>% 
  ggplot(aes(x=cnt_polar, y=pid, fill=grouping)) + 
  geom_bar(stat="identity", position="identity") + 
  facet_grid(vars(set_y), vars(set_x))

coef_order <- clix_clin_join %>% 
  filter(cohort=="neoseq", grouping=="x") %>% 
  group_by(pid) %>% 
  summarise(avg_coef=mean(coef)) %>% 
  arrange(avg_coef)

bind_rows(clix_clin_join,
          clix_ctak_join,
          clix_meta_join,
          clin_ctak_join,
          clin_meta_join,
          ctak_meta_join) %>% 
  filter(cohort=="neoseq") %>% 
  group_by(pid, set_x, set_y, grouping) %>% 
  summarise(avg_coef=mean(coef),
            cnt=n()) %>% 
  filter(grouping != "both") %>% 
  mutate(coef_polar=if_else(grouping=="x", -avg_coef, avg_coef),
         pid=factor(pid, levels=coef_order$pid),
         set_x = factor(set_x, levels=c("clix","clin","ctak"), labels=c("CLiX","ClinPhen","cTAKES")),
         set_y = factor(set_y, levels=c("clin","ctak","meta"), labels=c("ClinPhen","cTAKES","MetaMap"))) %>% 
  ggplot(aes(x=coef_polar, y=pid, fill=grouping)) + 
  geom_bar(stat="identity", position="identity") + 
  facet_grid(vars(set_y), vars(set_x))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Pairwise matrix of density distributions showing NLP set model coefficient comparisons
bind_rows(clix_clin_join,
          clix_ctak_join,
          clix_meta_join,
          clin_ctak_join,
          clin_meta_join,
          ctak_meta_join) %>% 
  filter(cohort=="neoseq", grouping!="both") %>% 
  mutate(set_x = factor(set_x, levels=c("clix","clin","ctak"), labels=c("CLiX","ClinPhen","cTAKES")),
         set_y = factor(set_y, levels=c("clin","ctak","meta"), labels=c("ClinPhen","cTAKES","MetaMap"))) %>% 
  ggplot(aes(x=coef, fill=grouping)) + 
  geom_density(alpha=0.6, linewidth=0.35) + 
  facet_grid(vars(set_y), vars(set_x))
## -------------------------------------------------------------------------- ##




clix_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_clinithink_stats.txt")
clin_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_clinphen_stats.txt")
ctak_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_ctakes_stats.txt")
meta_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_metamap_stats.txt")
medl_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_medlee_stats.txt")
manu_stats <- read_tsv("analysis/set_comparison/NeoSeq_manual_stats.txt")
icdc_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_icd10_stats.txt")
labs_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_labs_stats.txt")
ords_stats <- read_tsv("analysis/set_comparison/NeoSeq_edw_orders_stats.txt")

tool_stats <- bind_rows(list("clix"=clix_stats,
               "clin"=clin_stats,
               "ctak"=ctak_stats,
               "meta"=meta_stats,
               "medl"=medl_stats,
               "manu"=manu_stats,
               "icdc"=icdc_stats,
               "labs"=labs_stats,
               "ords"=ords_stats),
          .id="tool") %>% 
  mutate(tool = factor(tool, 
                       levels=c("clix","clin","ctak","meta","medl","manu","icdc","labs","ords"), 
                       labels=c("CLiX","ClinPhen","cTAKES","MetaMap","MedLEE","Manual","ICD10-CM","Labs","Orders")))

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Manual/NLP set term counts summary statistics
tool_stats %>% 
  group_by(tool) %>% 
  summarise(min_cnt = min(terminal_hpo_cnt + icd_cnt + other_cnt),
            q1_cnt = quantile(terminal_hpo_cnt + icd_cnt + other_cnt, probs=c(0.25)),
            median_cnt = median(terminal_hpo_cnt + icd_cnt + other_cnt),
            q3_cnt = quantile(terminal_hpo_cnt + icd_cnt + other_cnt, probs=c(0.75)),
            max_cnt = max(terminal_hpo_cnt + icd_cnt + other_cnt),
            mean_cnt = mean(terminal_hpo_cnt + icd_cnt + other_cnt),
            n = n())

# Manual/NLP set average information content summary statistics
tool_stats %>% 
  group_by(tool) %>% 
  summarise(min_avg_info = min(mean_info_cont),
            q1_avg_info = quantile(mean_info_cont, probs=c(0.25)),
            median_avg_info = median(mean_info_cont),
            q3_avg_info = quantile(mean_info_cont, probs=c(0.75)),
            max_avg_info = max(mean_info_cont),
            mean_avg_info = mean(mean_info_cont),
            n = n())

# Manual/NLP set total information content summary statistics
tool_stats %>% 
  group_by(tool) %>% 
  summarise(min_total_info = min(total_info_cont),
            q1_total_info = quantile(total_info_cont, probs=c(0.25)),
            median_total_info = median(total_info_cont),
            q3_total_info = quantile(total_info_cont, probs=c(0.75)),
            max_total_info = max(total_info_cont),
            mean_total_info = mean(total_info_cont),
            n = n())

# ANOVAs for term count, average info cont, total info cont
summary(aov(data=tool_stats, (terminal_hpo_cnt + icd_cnt + other_cnt) ~ tool))
summary(aov(data=tool_stats, mean_info_cont ~ tool))
summary(aov(data=tool_stats, total_info_cont ~ tool))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# T-tests comparing term counts between NLP sets
tool_stats %>% 
  filter(tool %in% c("CLiX","ClinPhen")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("CLiX","cTAKES")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("CLiX","MetaMap")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("CLiX","MedLEE")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("ClinPhen","cTAKES")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("ClinPhen","MetaMap")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("ClinPhen","MedLEE")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("cTAKES","MetaMap")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("cTAKES","MedLEE")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
tool_stats %>% 
  filter(tool %in% c("MetaMap","MedLEE")) %>% 
  mutate(tool = fct_drop(tool)) %>% 
  t.test(terminal_hpo_cnt ~ tool, data=.)
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Violin plot showing term counts for NLP sets
tool_stats %>% 
  filter(tool != "Manual") %>% 
  ggplot(aes(x=tool, y=(terminal_hpo_cnt + icd_cnt + other_cnt), fill=tool)) +
  geom_violin(alpha=0.5) + 
  geom_dotplot(binaxis="y",
               stackdir="center",
               dotsize=0.4) +
  scale_x_discrete(name="NLP Tool") +
  scale_y_continuous(name="HPO Term Count",
                     breaks=seq(0,300,50),
                     labels=seq(0,300,50),
                     limits=c(0,300)) +
  theme(legend.position="none",
        axis.title=element_text(size=18),
        axis.text=element_text(size=16))

# Violin plot showing average information contents for Manual/NLP sets
tool_stats %>% 
  filter(tool %in% c("CLiX","ClinPhen","cTAKES","MetaMap","MedLEE","Manual")) %>% 
  ggplot(aes(x=tool, y=mean_info_cont, fill=tool)) +
  geom_violin(alpha=0.5) + 
  geom_dotplot(binaxis="y",
               stackdir="center",
               dotsize=0.4) +
  scale_x_discrete(name="NLP Tool") +
  scale_y_continuous(name="Mean Information Content",
                     breaks=seq(2.0,6,0.5),
                     labels=seq(2.0,6,0.5),
                     limits=c(2.1,6)) +
  theme(legend.position="none",
        axis.title=element_text(size=18),
        axis.text=element_text(size=16))

# Violin plot showing total information contents for NLP sets
tool_stats %>% 
  filter(tool %in% c("CLiX","ClinPhen","cTAKES","MetaMap","MedLEE")) %>% 
  ggplot(aes(x=tool, y=total_info_cont, fill=tool)) +
  geom_violin(alpha=0.5) + 
  geom_dotplot(binaxis="y",
               stackdir="center",
               dotsize=0.4) +
  scale_x_discrete(name="NLP Tool") +
  scale_y_continuous(name="Total Information Content",
                     breaks=seq(0,1200,200),
                     labels=seq(0,1200,200),
                     limits=c(0,1200)) +
  theme(legend.position="none",
        axis.title=element_text(size=18),
        axis.text=element_text(size=16))
## -------------------------------------------------------------------------- ##




combined_meds_coef <- read_tsv("analysis/nlp_training_sets/meds/combined_meds_feature_coefficients.tsv")
meds_coef1 <- read_tsv("analysis/nlp_training_sets/meds/meds_feature_coefficients1.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef1 = coef)
meds_coef2 <- read_tsv("analysis/nlp_training_sets/meds/meds_feature_coefficients2.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef2 = coef)
meds_coef3 <- read_tsv("analysis/nlp_training_sets/meds/meds_feature_coefficients3.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef3 = coef)
meds_coef4 <- read_tsv("analysis/nlp_training_sets/meds/meds_feature_coefficients4.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef4 = coef)
meds_coef5 <- read_tsv("analysis/nlp_training_sets/meds/meds_feature_coefficients5.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef5 = coef)

combined_labs_coef <- read_tsv("analysis/nlp_training_sets/labs/combined_labs_feature_coefficients.tsv")
labs_coef1 <- read_tsv("analysis/nlp_training_sets/labs/labs_feature_coefficients1.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef1 = coef)
labs_coef2 <- read_tsv("analysis/nlp_training_sets/labs/labs_feature_coefficients2.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef2 = coef)
labs_coef3 <- read_tsv("analysis/nlp_training_sets/labs/labs_feature_coefficients3.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef3 = coef)
labs_coef4 <- read_tsv("analysis/nlp_training_sets/labs/labs_feature_coefficients4.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef4 = coef)
labs_coef5 <- read_tsv("analysis/nlp_training_sets/labs/labs_feature_coefficients5.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef5 = coef)

combined_ords_coef <- read_tsv("analysis/nlp_training_sets/ords/combined_ords_feature_coefficients.tsv")
ords_coef1 <- read_tsv("analysis/nlp_training_sets/ords/ords_feature_coefficients1.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef1 = coef)
ords_coef2 <- read_tsv("analysis/nlp_training_sets/ords/ords_feature_coefficients2.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef2 = coef)
ords_coef3 <- read_tsv("analysis/nlp_training_sets/ords/ords_feature_coefficients3.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef3 = coef)
ords_coef4 <- read_tsv("analysis/nlp_training_sets/ords/ords_feature_coefficients4.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef4 = coef)
ords_coef5 <- read_tsv("analysis/nlp_training_sets/ords/ords_feature_coefficients5.tsv") %>% 
  select(terms, coef) %>% 
  rename(coef5 = coef)

combined_icdc_coef <- read_tsv("analysis/nlp_training_sets/icdc/combined_icdc_feature_coefficients.tsv")
icdc_coef1 <- read_tsv("analysis/nlp_training_sets/icdc/icdc_feature_coefficients1.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef1 = coef)
icdc_coef2 <- read_tsv("analysis/nlp_training_sets/icdc/icdc_feature_coefficients2.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef2 = coef)
icdc_coef3 <- read_tsv("analysis/nlp_training_sets/icdc/icdc_feature_coefficients3.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef3 = coef)
icdc_coef4 <- read_tsv("analysis/nlp_training_sets/icdc/icdc_feature_coefficients4.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef4 = coef)
icdc_coef5 <- read_tsv("analysis/nlp_training_sets/icdc/icdc_feature_coefficients5.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef5 = coef)

combined_clix_coef <- read_tsv("analysis/nlp_training_sets/clix/combined_clix_feature_coefficients.tsv")
clix_coef1 <- read_tsv("analysis/nlp_training_sets/clix/clix_feature_coefficients1.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef1 = coef)
clix_coef2 <- read_tsv("analysis/nlp_training_sets/clix/clix_feature_coefficients2.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef2 = coef)
clix_coef3 <- read_tsv("analysis/nlp_training_sets/clix/clix_feature_coefficients3.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef3 = coef)
clix_coef4 <- read_tsv("analysis/nlp_training_sets/clix/clix_feature_coefficients4.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef4 = coef)
clix_coef5 <- read_tsv("analysis/nlp_training_sets/clix/clix_feature_coefficients5.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef5 = coef)

combined_clin_coef <- read_tsv("analysis/nlp_training_sets/clin/combined_clin_feature_coefficients.tsv")
clin_coef1 <- read_tsv("analysis/nlp_training_sets/clin/clin_feature_coefficients1.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef1 = coef)
clin_coef2 <- read_tsv("analysis/nlp_training_sets/clin/clin_feature_coefficients2.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef2 = coef)
clin_coef3 <- read_tsv("analysis/nlp_training_sets/clin/clin_feature_coefficients3.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef3 = coef)
clin_coef4 <- read_tsv("analysis/nlp_training_sets/clin/clin_feature_coefficients4.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef4 = coef)
clin_coef5 <- read_tsv("analysis/nlp_training_sets/clin/clin_feature_coefficients5.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef5 = coef)

combined_ctak_coef <- read_tsv("analysis/nlp_training_sets/ctak/combined_ctak_feature_coefficients.tsv")
ctak_coef1 <- read_tsv("analysis/nlp_training_sets/ctak/ctak_feature_coefficients1.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef1 = coef)
ctak_coef2 <- read_tsv("analysis/nlp_training_sets/ctak/ctak_feature_coefficients2.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef2 = coef)
ctak_coef3 <- read_tsv("analysis/nlp_training_sets/ctak/ctak_feature_coefficients3.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef3 = coef)
ctak_coef4 <- read_tsv("analysis/nlp_training_sets/ctak/ctak_feature_coefficients4.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef4 = coef)
ctak_coef5 <- read_tsv("analysis/nlp_training_sets/ctak/ctak_feature_coefficients5.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef5 = coef)

combined_meta_coef <- read_tsv("analysis/nlp_training_sets/meta/combined_meta_feature_coefficients.tsv")
meta_coef1 <- read_tsv("analysis/nlp_training_sets/meta/meta_feature_coefficients1.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef1 = coef)
meta_coef2 <- read_tsv("analysis/nlp_training_sets/meta/meta_feature_coefficients2.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef2 = coef)
meta_coef3 <- read_tsv("analysis/nlp_training_sets/meta/meta_feature_coefficients3.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef3 = coef)
meta_coef4 <- read_tsv("analysis/nlp_training_sets/meta/meta_feature_coefficients4.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef4 = coef)
meta_coef5 <- read_tsv("analysis/nlp_training_sets/meta/meta_feature_coefficients5.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef5 = coef)

combined_medl_coef <- read_tsv("analysis/nlp_training_sets/medl/combined_medl_feature_coefficients.tsv")
medl_coef1 <- read_tsv("analysis/nlp_training_sets/medl/medl_feature_coefficients1.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef1 = coef)
medl_coef2 <- read_tsv("analysis/nlp_training_sets/medl/medl_feature_coefficients2.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef2 = coef)
medl_coef3 <- read_tsv("analysis/nlp_training_sets/medl/medl_feature_coefficients3.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef3 = coef)
medl_coef4 <- read_tsv("analysis/nlp_training_sets/medl/medl_feature_coefficients4.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef4 = coef)
medl_coef5 <- read_tsv("analysis/nlp_training_sets/medl/medl_feature_coefficients5.tsv") %>% 
  select(codes, coef) %>% 
  rename(coef5 = coef)

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Model coefficients for various datasets
merged_meds_coef <- list(combined_meds_coef, meds_coef1, meds_coef2, meds_coef3, meds_coef4, meds_coef5) %>% 
  reduce(full_join, by="terms") %>% 
  arrange(desc(coef))

merged_labs_coef <- list(combined_labs_coef, labs_coef1, labs_coef2, labs_coef3, labs_coef4, labs_coef5) %>% 
  reduce(full_join, by="terms") %>% 
  arrange(desc(coef))

merged_ords_coef <- list(combined_ords_coef, ords_coef1, ords_coef2, ords_coef3, ords_coef4, ords_coef5) %>% 
  reduce(full_join, by="terms") %>% 
  arrange(desc(coef))

merged_icdc_coef <- list(combined_icdc_coef, icdc_coef1, icdc_coef2, icdc_coef3, icdc_coef4, icdc_coef5) %>% 
  reduce(full_join, by="codes") %>% 
  arrange(desc(coef))

merged_clix_coef <- list(combined_clix_coef, clix_coef1, clix_coef2, clix_coef3, clix_coef4, clix_coef5) %>% 
  reduce(full_join, by="codes") %>% 
  arrange(desc(coef))

merged_clin_coef <- list(combined_clin_coef, clin_coef1, clin_coef2, clin_coef3, clin_coef4, clin_coef5) %>% 
  reduce(full_join, by="codes") %>% 
  arrange(desc(coef))

merged_ctak_coef <- list(combined_ctak_coef, ctak_coef1, ctak_coef2, ctak_coef3, ctak_coef4, ctak_coef5) %>% 
  reduce(full_join, by="codes") %>% 
  arrange(desc(coef))

merged_meta_coef <- list(combined_meta_coef, meta_coef1, meta_coef2, meta_coef3, meta_coef4, meta_coef5) %>% 
  reduce(full_join, by="codes") %>% 
  arrange(desc(coef))

merged_medl_coef <- list(combined_medl_coef, medl_coef1, medl_coef2, medl_coef3, medl_coef4, medl_coef5) %>% 
  reduce(full_join, by="codes") %>% 
  arrange(desc(coef))
## -------------------------------------------------------------------------- ##




training <- read_tsv("analysis/nlp_training_sets/training_preds_combined.tsv",
                     col_types = cols(
                       mrn=col_character(),
                       dwid=col_character()
                     ))
testing <- read_tsv("analysis/nlp_testing_sets/testing_preds_combined.tsv")

permuted_training <- read_tsv("analysis/nlp_training_sets/permuted_seq_status/permuted_training_preds_combined.tsv")
permuted_testing <- read_tsv("analysis/nlp_testing_sets/permuted_seq_status/permuted_testing_preds_combined.tsv")

joined_training <- read_tsv("analysis/nlp_training_sets/joined_sets/joined_training_preds_combined.tsv")

physician_terms <- read_tsv("data/complete_datasets/NeoSeq_physician_HPO.tsv")

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# orthogonal NLP accuracy test using physician-generated "ground truth"
neo_term_sets <- training %>% 
  select(pid:codes) %>% 
  distinct() %>% 
  bind_rows(physician_terms) %>% 
  filter(cohort=="neoseq") %>% 
  filter(codes_source %in% c("edw_clix","edw_clinphen","edw_ctakes","edw_metamap","edw_medlee","physician")) %>% 
  select(pid, codes_source, codes) %>% 
  filter(pid %in% unique(training$pid)) %>% 
  mutate(codes = str_split(codes, pattern=";"),
         codes_p = map(codes, \(x) get_ancestors(hpo, x)))
ground_truth <- neo_term_sets %>% 
  select(-codes_p) %>% 
  pivot_wider(names_from=codes_source, values_from=codes) %>% 
  mutate(phy_cnt = map_int(physician, length),
         clix_hits = map2_int(physician, edw_clix, \(x,y) length(intersect(x,y))) / phy_cnt,
         clin_hits = map2_int(physician, edw_clinphen, \(x,y) length(intersect(x,y))) / phy_cnt,
         ctak_hits = map2_int(physician, edw_ctakes, \(x,y) length(intersect(x,y))) / phy_cnt,
         meta_hits = map2_int(physician, edw_metamap, \(x,y) length(intersect(x,y))) / phy_cnt,
         medl_hits = map2_int(physician, edw_medlee, \(x,y) length(intersect(x,y))) / phy_cnt)
ground_truth_p <- neo_term_sets %>% 
  select(-codes) %>% 
  pivot_wider(names_from=codes_source, values_from=codes_p) %>% 
  mutate(phy_cnt = map_int(physician, length),
         clix_hits = map2_int(physician, edw_clix, \(x,y) length(intersect(x,y))) / phy_cnt,
         clin_hits = map2_int(physician, edw_clinphen, \(x,y) length(intersect(x,y))) / phy_cnt,
         ctak_hits = map2_int(physician, edw_ctakes, \(x,y) length(intersect(x,y))) / phy_cnt,
         meta_hits = map2_int(physician, edw_metamap, \(x,y) length(intersect(x,y))) / phy_cnt,
         medl_hits = map2_int(physician, edw_medlee, \(x,y) length(intersect(x,y))) / phy_cnt)

ground_truth %>% 
  select(clix_hits, clin_hits, ctak_hits, meta_hits, medl_hits) %>% 
  pivot_longer(cols=ends_with("_hits"), names_to="tool", values_to="accuracy") %>% 
  group_by(tool) %>% 
  summarise(min = min(accuracy),
            q1 = quantile(accuracy, probs=c(0.25)),
            median = median(accuracy),
            q3 = quantile(accuracy, probs=c(0.75)),
            max = max(accuracy),
            mean = mean(accuracy),
            n = n())

ground_truth_p %>% 
  select(clix_hits, clin_hits, ctak_hits, meta_hits, medl_hits) %>% 
  pivot_longer(cols=ends_with("_hits"), names_to="tool", values_to="accuracy") %>% 
  group_by(tool) %>% 
  summarise(min = min(accuracy),
            q1 = quantile(accuracy, probs=c(0.25)),
            median = median(accuracy),
            q3 = quantile(accuracy, probs=c(0.75)),
            max = max(accuracy),
            mean = mean(accuracy),
            n = n())

ground_truth_p %>% 
  select(clix_hits, clin_hits, ctak_hits, meta_hits, medl_hits) %>% 
  pivot_longer(cols=ends_with("_hits"), names_to="tool", values_to="accuracy") %>% 
  mutate(tool = factor(tool, levels=c("meta_hits","medl_hits","clix_hits","ctak_hits","clin_hits"))) %>%
  ggplot(aes(x=accuracy, fill=tool)) + 
  geom_density(alpha=0.4) +
  scale_y_continuous(name="Density",
                     minor_breaks=NULL) +
  scale_x_continuous(name="Ground Truth Accuracy",
                     limits=c(0, 1),
                     breaks=seq(0, 1, 0.2),
                     minor_breaks=NULL) +
  scale_fill_manual(values = c("pink","red","yellow","green","blue")) +
  theme_bw()
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# training data AUC for MPSE model trained with CLiX
clix_train <- training %>% 
  filter(codes_source=="edw_clix") %>% 
  select(seq_status, pos_proba)
auc(roc(clix_train$seq_status, clix_train$pos_proba))

# training data AUC for MPSE model trained with ClinPhen
clin_train <- training %>% 
  filter(codes_source=="edw_clinphen") %>% 
  select(seq_status, pos_proba)
auc(roc(clin_train$seq_status, clin_train$pos_proba))

# training data AUC for MPSE model trained with cTAKES
ctak_train <- training %>% 
  filter(codes_source=="edw_ctakes") %>% 
  select(seq_status, pos_proba)
auc(roc(ctak_train$seq_status, ctak_train$pos_proba))

# training data AUC for MPSE model trained with MetaMap
meta_train <- training %>% 
  filter(codes_source=="edw_metamap") %>% 
  select(seq_status, pos_proba)
auc(roc(meta_train$seq_status, meta_train$pos_proba))

# training data AUC for MPSE model trained with MedLEE
medl_train <- training %>% 
  filter(codes_source=="edw_medlee") %>% 
  select(seq_status, pos_proba)
auc(roc(medl_train$seq_status, medl_train$pos_proba))

# training data AUC for MPSE model trained with ICD10 codes
icdc_train <- training %>% 
  filter(codes_source=="edw_icd10") %>% 
  select(seq_status, pos_proba)
auc(roc(icdc_train$seq_status, icdc_train$pos_proba))

# training data AUC for MPSE model trained with labs codes
labs_train <- training %>% 
  filter(codes_source=="edw_labs") %>% 
  select(seq_status, pos_proba)
auc(roc(labs_train$seq_status, labs_train$pos_proba))

# training data AUC for MPSE model trained with meds codes
meds_train <- training %>% 
  filter(codes_source=="edw_meds") %>% 
  select(seq_status, pos_proba)
auc(roc(meds_train$seq_status, meds_train$pos_proba))

# training data AUC for MPSE model trained with orders codes
ords_train <- training %>% 
  filter(codes_source=="edw_orders") %>% 
  select(seq_status, pos_proba)
auc(roc(ords_train$seq_status, ords_train$pos_proba))


# joined training data AUC for MPSE model trained with CLiX AND ICD10 codes
clix_icdc_train <- joined_training %>% 
  filter(codes_source=="clix_icd10") %>% 
  select(seq_status, pos_proba)
auc(roc(clix_icdc_train$seq_status, clix_icdc_train$pos_proba))

# joined training data AUC for MPSE model trained with ClinPhen AND ICD10 codes
clin_icdc_train <- joined_training %>% 
  filter(codes_source=="clin_icd10") %>% 
  select(seq_status, pos_proba)
auc(roc(clin_icdc_train$seq_status, clin_icdc_train$pos_proba))

# joined training data AUC for MPSE model trained with cTAKES AND ICD10 codes
ctak_icdc_train <- joined_training %>% 
  filter(codes_source=="ctak_icd10") %>% 
  select(seq_status, pos_proba)
auc(roc(ctak_icdc_train$seq_status, ctak_icdc_train$pos_proba))

# joined training data AUC for MPSE model trained with MetaMap AND ICD10 codes
meta_icdc_train <- joined_training %>% 
  filter(codes_source=="meta_icd10") %>% 
  select(seq_status, pos_proba)
auc(roc(meta_icdc_train$seq_status, meta_icdc_train$pos_proba))

# joined training data AUC for MPSE model trained with MedLEE AND ICD10 codes
medl_icdc_train <- joined_training %>% 
  filter(codes_source=="medl_icd10") %>% 
  select(seq_status, pos_proba)
auc(roc(medl_icdc_train$seq_status, medl_icdc_train$pos_proba))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# permuted training data AUC for MPSE model trained with CLiX
perm_clix_train <- permuted_training %>% 
  filter(codes_source=="edw_clix") %>% 
  select(seq_status, pos_proba)
auc(roc(perm_clix_train$seq_status, perm_clix_train$pos_proba))

# permuted training data AUC for MPSE model trained with ClinPhen
perm_clin_train <- permuted_training %>% 
  filter(codes_source=="edw_clinphen") %>% 
  select(seq_status, pos_proba)
auc(roc(perm_clin_train$seq_status, perm_clin_train$pos_proba))

# permuted training data AUC for MPSE model trained with cTAKES
perm_ctak_train <- permuted_training %>% 
  filter(codes_source=="edw_ctakes") %>% 
  select(seq_status, pos_proba)
auc(roc(perm_ctak_train$seq_status, perm_ctak_train$pos_proba))

# permuted training data AUC for MPSE model trained with MetaMap
perm_meta_train <- permuted_training %>% 
  filter(codes_source=="edw_metamap") %>% 
  select(seq_status, pos_proba)
auc(roc(perm_meta_train$seq_status, perm_meta_train$pos_proba))

# permuted training data AUC for MPSE model trained with MedLEE
perm_medl_train <- permuted_training %>% 
  filter(codes_source=="edw_medlee") %>% 
  select(seq_status, pos_proba)
auc(roc(perm_medl_train$seq_status, perm_medl_train$pos_proba))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# case/control MPSE score density distributions for UofU training data from
# model trained with different NLP tools
training %>% 
  select(cohort, codes_source, scr) %>% 
  mutate(cohort = factor(cohort, 
                         levels=c("utah","neoseq"),
                         labels=c("UofU Not Sequenced (n=240)", "UofU NeoSeq (n=60)")),
         codes_source = factor(codes_source, 
                               levels=c("edw_clix","edw_clinphen","edw_ctakes","edw_metamap","edw_medlee","edw_icd10","edw_meds","edw_labs","edw_orders"), 
                               labels=c("CLiX","ClinPhen","cTAKES","MetaMap","MedLEE","ICD10-CM","Meds","Labs","Orders"))) %>%
  ggplot(aes(x=scr, fill=cohort)) + 
  geom_density(alpha=0.6, linewidth=0.35) + 
  facet_wrap(vars(codes_source), nrow=5, ncol=2, scales="free_y") + 
  scale_x_continuous(name="MPSE Score",
                     breaks=seq(-100,300,50),
                     labels=seq(-100,300,50)) + 
  scale_y_continuous(name="Density") +
  scale_fill_manual(name="Cohort", 
                    guide = guide_legend(title=NULL),
                    values=c("UofU Not Sequenced (n=240)"="#00BFC4", 
                             "UofU NeoSeq (n=60)"="#F8766D")) +
  theme_bw() +
  theme(legend.text=element_text(size=14),
        legend.position=c(0.75, 0.1),
        axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        strip.text=element_text(size=14))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# AUC values for MPSE models trained with NLP, ICD10, and NLP+ICD10 data
joined_auc <- tribble(
  ~Source, ~group, ~AUROC,
  "meds", "base", 0.822,
  "labs", "base", 0.844,
  "ords", "base", 0.855,
  "icd10", "base", 0.911,
  "clix", "base", 0.909,
  "clix", "join", 0.917,
  "clin", "base", 0.896,
  "clin", "join", 0.913,
  "ctak", "base", 0.891,
  "ctak", "join", 0.907,
  "meta", "base", 0.833,
  "meta", "join", 0.887,
  "medl", "base", 0.827,
) %>% 
  mutate(Source = factor(Source, 
                         levels=c("meds","labs","ords","icd10","clix","clin","ctak","meta","medl"), 
                         labels=c("Meds","Labs","Orders","ICD10-CM", "CLiX", "ClinPhen", "cTAKES", "MetaMap","MedLEE")),
         group = factor(group, levels=c("base","join"), labels=c("Source", "Source+ICD10")))

ggplot(joined_auc, aes(x=Source, y=AUROC, fill=group)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  geom_text(aes(label=AUROC), vjust=1.6, color="black", position=position_dodge(0.9), size=4) + 
  scale_fill_brewer(guide = guide_legend(title=NULL),
                    palette="Paired") + 
  coord_cartesian(ylim=c(0.8,0.95)) +
  theme_minimal() +
  theme(legend.text=element_text(size=12),
        axis.title=element_text(size=16),
        axis.text=element_text(size=14))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# CLiX testing data AUC for MPSE model trained with CLiX
clix_test <- testing %>% 
  filter(codes_source=="edw_clix") %>% 
  select(seq_status, pos_proba)
auc(roc(clix_test$seq_status, clix_test$pos_proba))

# ClinPhen testing data AUC for MPSE model trained with CLiX
clin_test <- testing %>% 
  filter(codes_source=="edw_clinphen") %>% 
  select(seq_status, pos_proba)
auc(roc(clin_test$seq_status, clin_test$pos_proba))

# cTAKES testing data AUC for MPSE model trained with CLiX
ctak_test <- testing %>% 
  filter(codes_source=="edw_ctakes") %>% 
  select(seq_status, pos_proba)
auc(roc(ctak_test$seq_status, ctak_test$pos_proba))

# MetaMap testing data AUC for MPSE model trained with CLiX
meta_test <- testing %>% 
  filter(codes_source=="edw_metamap") %>% 
  select(seq_status, pos_proba)
auc(roc(meta_test$seq_status, meta_test$pos_proba))

# MedLEE testing data AUC for MPSE model trained with CLiX
medl_test <- testing %>% 
  filter(codes_source=="edw_medlee") %>% 
  select(seq_status, pos_proba)
auc(roc(medl_test$seq_status, medl_test$pos_proba))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# case/control MPSE score density distributions for UofU testing data from
# model trained with CLiX only
testing %>% 
  select(cohort, codes_source, scr) %>% 
  mutate(cohort = factor(cohort, 
                         levels=c("utah","neoseq"),
                         labels=c("UofU Not Sequenced (n=240)", "UofU NeoSeq (n=60)")),
         codes_source = factor(codes_source, 
                               levels=c("edw_clix","edw_clinphen","edw_ctakes","edw_metamap","edw_medlee"), 
                               labels=c("CLiX","ClinPhen","cTAKES","MetaMap","MedLEE"))) %>% 
  ggplot(aes(x=scr, fill=cohort)) + 
  geom_density(alpha=0.6, linewidth=0.35) + 
  facet_wrap(vars(codes_source), nrow=5, ncol=1, scales="free_y") + 
  scale_x_continuous(name="MPSE Score",
                     breaks=seq(-75,150,25),
                     labels=seq(-75,150,25)) + 
  scale_y_continuous(name="Density") +
  scale_fill_manual(name="Cohort", 
                    guide = guide_legend(title=NULL),
                    values=c("UofU Not Sequenced (n=240)"="#00BFC4", 
                             "UofU NeoSeq (n=60)"="#F8766D")) +
  theme_bw() +
  theme(legend.text=element_text(size=14),
        legend.position=c(0.75, 0.1),
        axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        strip.text=element_text(size=14))
## -------------------------------------------------------------------------- ##




gem_benchmark <- read_delim("analysis/GEM_benchmark/diagnostic_gene_rank.tsv", delim=" ")

## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# bar chart of diagnostic gene ranks across Manual/NLP tools
gem_benchmark %>% 
  mutate(top1 = if_else(GEM_RANK==1, "y", "n"),
         top2 = if_else(GEM_RANK<=2, "y", "n"),
         top5 = if_else(GEM_RANK<=5, "y", "n"),
         top10 = if_else(GEM_RANK<=10, "y", "n")) %>% 
  pivot_longer(cols=top1:top10,
               names_to="rank_bin",
               values_to="rank_status") %>% 
  mutate(cohort = factor(cohort, 
                         levels=c("clix","clin","ctak","meta","manu"), 
                         labels=c("CLiX","ClinPhen","cTAKES","MetaMap","Manual")),
         rank_bin = factor(rank_bin, 
                           levels=c("top1","top2","top5","top10"), 
                           labels=c("Top 1", "Top 2", "Top 5", "Top 10"))) %>% 
  filter(rank_status=="y") %>% 
  ggplot(aes(x=rank_bin, fill=cohort)) + 
  geom_bar(position=position_dodge()) + 
  scale_x_discrete(name="Diagnostic gene rank") +
  scale_y_continuous(name="Percentage of cases",
                     breaks=seq(0,20,2),
                     labels=c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")) + 
  scale_fill_brewer(name="Tool",
                    palette="Set2") + 
  theme_linedraw() + 
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=15),
        axis.title=element_text(size=18),
        axis.text=element_text(size=16))
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Dot plot of diagnostic gene ranks across Manual/NLP tools
gem_benchmark %>% 
  mutate(cohort = factor(cohort, 
                         levels=c("clix","clin","ctak","meta","manu"), 
                         labels=c("CLiX","ClinPhen","cTAKES","MetaMap","Manual"))) %>% 
  ggplot(aes(x=cohort, y=-GEM_RANK, color=cohort)) + 
  geom_jitter(width=0.1, height=0, shape=1, size=2, stroke=1.2) + 
  geom_abline(slope=0, intercept=0, linetype=3) +
  scale_x_discrete(name=NULL) +
  scale_y_continuous(name="Rank",
                     limits=c(-20,0),
                     breaks=seq(-20,0,5),
                     labels=seq(20,0,-5)) + 
  scale_color_brewer(palette="Set2") + 
  theme_classic()
## -------------------------------------------------------------------------- ##


## -------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------- ##
# Dot plot of diagnostic gene Bayes Factors across Manual/NLP tools
gem_benchmark %>% 
  mutate(cohort = factor(cohort, 
                         levels=c("clix","clin","ctak","meta","manu"), 
                         labels=c("CLiX","ClinPhen","cTAKES","MetaMap","Manual"))) %>% 
  ggplot(aes(x=cohort, y=GEM_SCORE, color=cohort)) + 
  geom_jitter(width=0.1, height=0, shape=1, size=2, stroke=1.2) + 
  geom_abline(slope=0, intercept=0, linetype=3) +
  stat_summary(aes(x=cohort, y=GEM_SCORE, group=cohort), fun="median", geom="crossbar", width=0.4) + 
  scale_x_discrete(name=NULL) +
  scale_y_continuous(name="Bayes Factor",
                     limits=c(0,4),
                     breaks=seq(0,4,1),
                     labels=seq(0,4,1)) + 
  scale_color_brewer(palette="Set2") + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title=element_text(size=18),
        axis.text=element_text(size=16))
## -------------------------------------------------------------------------- ##


