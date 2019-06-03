options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(purrr)
library(stringr)
library(stats)

data <- read.csv("./Cleaned/20181016_AT1_mito_tidy_peptide_stoich_output.csv")

# Condition
data$condition <- NA
data$condition[which(data$sample_id == 10)] <- "KI_6"
data$condition[which(data$sample_id == 11)] <- "KI_6"
data$condition[which(data$sample_id == 12)] <- "KI_6"
data$condition[which(data$sample_id == 13)] <- "KI_6"
data$condition[which(data$sample_id == 20)] <- "WT_3"
data$condition[which(data$sample_id == 21)] <- "WT_3"
data$condition[which(data$sample_id == 22)] <- "WT_3"
data$condition[which(data$sample_id == 23)] <- "WT_3"
data$condition[which(data$sample_id == 30)] <- "TG_3"
data$condition[which(data$sample_id == 31)] <- "TG_3"
data$condition[which(data$sample_id == 32)] <- "TG_3"
data$condition[which(data$sample_id == 33)] <- "TG_3"
data$condition[which(data$sample_id == 40)] <- "WT_6"
data$condition[which(data$sample_id == 41)] <- "WT_6"
data$condition[which(data$sample_id == 42)] <- "WT_6"
data$condition[which(data$sample_id == 43)] <- "WT_6"

data$subcell_fraction <- "mito"

### Preparing the data for ANOVA analysis
ms <- data %>% select(condition, sample_id, PG.ProteinGroups, PG.Genes, EG.StrippedSequence,
                      EG.ModifiedSequence, k_site, stoich_corr)
ms$condition_num <- 0
ms$condition_num[which(ms$condition == "KI_6")] <- 2
ms$condition_num[which(ms$condition == "WT_6")] <- 1
ms$condition_num[which(ms$condition == "TG_3")] <- 4
ms$condition_num[which(ms$condition == "WT_3")] <- 3

ms$condition_num <- as.numeric(ms$condition_num)

ms$name <- paste(ms$PG.ProteinGroups, ms$k_site, ms$EG.ModifiedSequence, sep = "_")

mss <- split(ms, paste(ms$PG.ProteinGroups, ms$k_site, ms$EG.ModifiedSequence, sep = "_"))

# Looking for sites that change at any time point
tmpfn <- function(dat) {
  dat <- dat %>% filter(condition_num == 1 | condition_num == 2)
  if(length(unique(dat$condition_num)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(condition_num), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval_KI <- anova(fit)[,"Pr(>F)"][1]
  else
    pval_KI <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval_KI)
}
(KI_sig <- mss %>%
    map(tmpfn) %>%
    bind_rows(.id = "name") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval_KI)) %>%
  filter(pval_KI <= .1)


###########################################################################
tmpfn <- function(dat) {
  dat <- dat %>% filter(condition_num == 3 | condition_num == 4)
  if(length(unique(dat$condition_num)) < 2)
    return(NULL)
  fit <- lm(stoich_corr ~ factor(condition_num), dat)
  rsq <- summary(fit)$r.squared
  if(rsq < 0.99)
    pval_TG <- anova(fit)[,"Pr(>F)"][1]
  else
    pval_TG <- NA
  coefs <- t(coef(fit))
  data.frame(coefs,rsq, pval_TG)
}
(TG_sig <- mss %>%
    map(tmpfn) %>%
    bind_rows(.id = "name") %>%
    mutate_if(is.numeric, function(x) round(x, 4)) %>%
    arrange(pval_TG)) %>%
  filter(pval_TG <= .1)


#########################################################
# Make a wide version
data$bio_rep <- NA
data$bio_rep[which(data$sample_id == 10)] <- 1
data$bio_rep[which(data$sample_id == 11)] <- 2
data$bio_rep[which(data$sample_id == 12)] <- 3
data$bio_rep[which(data$sample_id == 13)] <- 4
data$bio_rep[which(data$sample_id == 20)] <- 1
data$bio_rep[which(data$sample_id == 21)] <- 2
data$bio_rep[which(data$sample_id == 22)] <- 3
data$bio_rep[which(data$sample_id == 23)] <- 4
data$bio_rep[which(data$sample_id == 30)] <- 1
data$bio_rep[which(data$sample_id == 31)] <- 2
data$bio_rep[which(data$sample_id == 32)] <- 3
data$bio_rep[which(data$sample_id == 33)] <- 4
data$bio_rep[which(data$sample_id == 40)] <- 1
data$bio_rep[which(data$sample_id == 41)] <- 2
data$bio_rep[which(data$sample_id == 42)] <- 3
data$bio_rep[which(data$sample_id == 43)] <- 4


data <- data %>% 
  select(condition, sample_id, bio_rep, everything())  %>% 
  arrange(PG.ProteinGroups, k_site, k_count, condition, bio_rep)

summary <- data %>% 
  group_by(condition, PG.ProteinGroups, PG.ProteinDescriptions, PG.Genes, 
           k_site, k_count, EG.StrippedSequence, EG.ModifiedSequence) %>% 
  summarise(stoich_mean = mean(stoich_corr), 
            stoich_sd = sd(stoich_corr),
            stoich_n = n()) %>% 
  arrange(PG.ProteinGroups, k_site, condition) %>% 
  ungroup()

summary$name <- paste(summary$PG.ProteinGroups, summary$k_site, summary$EG.ModifiedSequence, sep = "_")

###########################################################################
# To merge the pvalue to the previous dataframe to see things that are significantly changing
KI_sig <- KI_sig %>% select(name, rsq, pval_KI)
names(KI_sig)[2] <- "rsq_KI"

TG_sig <- TG_sig %>% select(name, rsq, pval_TG)
names(TG_sig)[2] <- "rsq_TG"

data_stats <- merge(summary, KI_sig, by = "name", all = TRUE)
data_stats <- merge(data_stats, TG_sig, by = "name", all = TRUE)

rm(KI_sig, TG_sig, ms, mss, summary);gc()

####
wide_mean <- data_stats %>%
  select(PG.ProteinGroups, k_site, PG.Genes, PG.ProteinDescriptions,
         EG.ModifiedSequence, EG.StrippedSequence, k_count, condition, 
         pval_KI, pval_TG, stoich_mean) %>% 
  spread(condition, stoich_mean)
names(wide_mean)[10] <- "KI6_mean"
names(wide_mean)[11] <- "TG3_mean"
names(wide_mean)[12] <- "WT3_mean"
names(wide_mean)[13] <- "WT6_mean"

wide_sd <- data_stats %>%
  select(PG.ProteinGroups, k_site, PG.Genes, PG.ProteinDescriptions,
         EG.ModifiedSequence, EG.StrippedSequence, k_count, condition, 
         pval_KI, pval_TG, stoich_sd) %>% 
  spread(condition, stoich_sd)
names(wide_sd)[10] <- "KI6_sd"
names(wide_sd)[11] <- "TG3_sd"
names(wide_sd)[12] <- "WT3_sd"
names(wide_sd)[13] <- "WT6_sd"

wide_n <- data_stats %>%
  select(PG.ProteinGroups, k_site, PG.Genes, PG.ProteinDescriptions, 
         EG.ModifiedSequence, EG.StrippedSequence, k_count, condition, 
         pval_KI, pval_TG, stoich_n) %>% 
  spread(condition, stoich_n)
names(wide_n)[10] <- "KI6_n"
names(wide_n)[11] <- "TG3_n"
names(wide_n)[12] <- "WT3_n"
names(wide_n)[13] <- "WT6_n"

wide <- merge(wide_mean, wide_sd, all = TRUE)
wide <- merge(wide, wide_n, all = TRUE)

rm(wide_mean, wide_sd, wide_n); gc()

wide$KI6_WT6 <- wide$KI6_mean - wide$WT6_mean
wide$TG3_WT3 <- wide$TG3_mean - wide$WT3_mean

wide$KI6_WT6 <- round(wide$KI6_WT6, 4)
wide$TG3_WT3 <- round(wide$TG3_WT3, 4)

#write.table(wide, "./Cleaned/20181219_AT1_mito_clean_withstats.csv", sep = ",", row.names = FALSE)
