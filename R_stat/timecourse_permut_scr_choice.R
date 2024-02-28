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
library(lmPerm)
#install.packages("devtools")  ## if not already installed
#devtools::install_github("mtorchiano/lmPerm")


sterr <- function(x) sd(x)/sqrt(length(x))
path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"
out_path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"
filename <- paste0(path,'df_LMEMRTtimetrainedpostrisk.csv')

df_large = read.csv(filename)
colnames(df_large)
unique(df_large$trial_type)
unique(df_large$subject)
df_large$X = NULL


#prepare for Z
ncol_min <- grep('scr...1....0.9.',colnames(df_large))
ncol_max <- grep('scr..2.8.2.9.',colnames(df_large))
cols <- c(colnames(df_large)[ncol_min:ncol_max])

df_filt_m <- data.table(melt(df_large, id.vars = c('subject', 'real_time'), measure.vars = cols))
Zmeans <- df_filt_m[, mean(value), by=subject]
setnames(Zmeans,'V1','Zmean_2')

Zsds <- df_filt_m[, sd(value), by=subject]
df_filt_m[, sd:=sd(value), by=subject]
setnames(Zsds,'V1','Zsd_2')

df_large <- merge(df_large,Zmeans,by='subject',all.x = TRUE)
df_large <- merge(df_large,Zsds,by='subject',all.x = TRUE)


#apply Z for scr
cols <- colnames(df_large)[grep('scr..',colnames(df_large))]
for (j in cols) set(df_large, j = j, value = (df_large[[j]]-df_large$Zmean_2)/df_large$Zsd_2)

df_large = data.table(df_large)
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
cols
ncol_min <- grep('scr...0.9..0.8.',colnames(df_large))
ncol_max <- grep('scr..2.8.2.9.',colnames(df_large))
cols2 <- c(colnames(df_large)[ncol_min:ncol_max])
cols2
df_large_prev = df_large
df_large_prev = df_large_prev


for (j in cols) {for (k in cols2){df_large[[j]] = (df_large_prev[[j]] - df_large_prev[[k]])/0.1}}

mode = 'risk'
df_large = df_large[trial_type != 'postrisk']
unique(df_large$trial_type)
############################stat################################################

temp <- df_large

######## for trial_type #############
p_vals <- data.table()

temp$subject <- as.factor(temp$subject)
temp$round <- as.factor(temp$round)
temp$trial_type <- as.factor(temp$trial_type)

unique(temp$trial_type)

#Fitting and testing linear models with permutation tests
for (j in cols){
  m <- lmp(get(j) ~ trial_type, data = temp,  perm="Prob", maxIter = 5000 )

  a <- as.data.table(coef(summary(m)))
  
  a$`Pr(Prob)`[-1] <- format(a$`Pr(Prob)`[-1], digits = 3)
  a$`Pr(Prob)`[-1]
  
  an <- data.table()
  an$interval <- j
  an$trial_type <- format(a$`Pr(Prob)`[-1], digits = 3)
  
  p_vals <- rbind(p_vals,an)
  p_vals
}
#fdr correction of permut pvall

p_vals <- p_vals[, trial_type_fdr:=p.adjust(`trial_type`, method = 'fdr')] 

write.csv(p_vals,paste0(out_path, "p_vals_factor_significance_Z_deriv_scr_adols_permut_', mode, '_trained.csv"))