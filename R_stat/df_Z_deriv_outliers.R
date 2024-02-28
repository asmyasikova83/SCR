# LMM analysis separately on trained and not_trained or both.
# Choose mode to operate the subset you want  

#rm(list = ls())


library(reshape2)
library(data.table)
library(ggplot2)
library(lme4)
library("ggpubr")
library(emmeans)
library(lmerTest)
library(dplyr)
library(tidyverse)
subjects = c("P701", "P702", "P703", "P704", "P705", "P706", "P707", "P708", "P709",
            "P710", "P711", "P712", "P713", "P715", "P714","P716", "P717", "P718",
            "P719", "P720", "P721", "P722", "P723", "P724", "P725", "P726", "P727",
            "P728", "P729", "P730")

path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"
out_path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"

filename <- paste0(path,'df_LMEMRTtimetrainedpostrisk.csv')
df_large = read.csv(filename)
df_large <- as.data.table(df_large)

#set trial_type
mode = 'postrisk'

df_large = df_large[trial_type == as.character(mode)]
unique(df_large$trial_type)

df_large$X = NULL

####outliers###################################
ncol_min <- grep('scr...1....0.9.',colnames(df_large))
ncol_max <- grep('scr..2.8.2.9.',colnames(df_large))
cols <- c(colnames(df_large)[ncol_min:ncol_max])
cols

df_filt_m <- data.table(melt(df_large, id.vars = c('subject', 'real_time'), measure.vars = cols))
Zmeans <- df_filt_m[, mean(value), by=subject]
setnames(Zmeans,'V1','Zmean_2')

Zsds <- df_filt_m[, sd(value), by=subject]
df_filt_m[, sd:=sd(value), by=subject]
setnames(Zsds,'V1','Zsd_2')

df_large <- merge(df_large,Zmeans,by='subject',all.x = TRUE)
df_large <- merge(df_large,Zsds,by='subject',all.x = TRUE)

colnames(df_large)

#apply Z for scr
cols <- colnames(df_large)[grep('scr..',colnames(df_large))]
for (j in cols) set(df_large, j = j, value = (df_large[[j]]-df_large$Zmean_2)/df_large$Zsd_2)

#outliers
for (j in cols){
  df_large[,flagmax:= df_large[[j]] < (mean(df_large[[j]]) + 3*sd(df_large[[j]]))]
  df_large = df_large[flagmax==TRUE]
  df_large[,flagmin:= df_large[[j]] > (mean(df_large[[j]]) - 3*sd(df_large[[j]]))]
  df_large = df_large[flagmin==TRUE]
}

#derivative

ncol_min <- grep('scr...1....0.9.',colnames(df_large))
ncol_max <- grep('scr..2.8.2.9.',colnames(df_large))
cols <- c(colnames(df_large)[ncol_min:ncol_max])
ncol_min <- grep('scr...0.9..0.8.',colnames(df_large))
ncol_max <- grep('scr..2.8.2.9.',colnames(df_large))
cols2 <- c(colnames(df_large)[ncol_min:ncol_max])
df_large_prev = df_large

for (j in cols) {for (k in cols2){df_large[[j]] = (df_large_prev[[j]] - df_large_prev[[k]])/0.1}}

write.csv(df_large,paste0(out_path, "df_adols_ZZ_', mode, '_trained_deriv.csv"))
