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

path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"
out_path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/pics/"
filename <- paste0(path,'df_LMEMRTtimetrainedpostrisk.csv')

dftrained = read.csv(filename)
dftrained$train = 'trained'
df_large = dftrained
df_large <- as.data.table(df_large)

#df_large = df_large[trial_type != 'risk']
############################stat################################################


####рreрrocessing###################################
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

#average over the time interval
#decision making
ncol_min <- grep('scr...1....0.9.',colnames(df_large))
ncol_max <- grep('scr...0.2..0.1.',colnames(df_large))
cols <- c(colnames(df_large)[ncol_min:ncol_max])

min = '1....0.9.'
max = '0.2..0.1.'

colname <- paste0('resp_','scr..',min, '_','scr..',max)

df_large <- as.data.table(df_large)
df_large[,(colname):=rowMeans(df_large[,.SD,.SDcol=cols])]
colname

scr_postrisk <- df_large[trial_type == 'postrisk']$resp_scr..1....0.9._scr..0.2..0.1.
scr_norisk <- df_large[trial_type == 'norisk']$resp_scr..1....0.9._scr..0.2..0.1.
scr_risk <- df_large[trial_type == 'risk']$resp_scr..1....0.9._scr..0.2..0.1.

scr <- df_large$resp_scr..1....0.9._scr..0.2..0.1.

RT_norisk <- df_large[trial_type == 'norisk']$response_time
RT_postrisk <- df_large[trial_type == 'postrisk']$response_time
RT_risk <- df_large[trial_type == 'risk']$response_time

RT <- df_large$response_time

#combined = as.data.table(cbind(RT_postrisk, scr_postrisk))
#setnames(combined,c("RT_postrisk","scr_postrisk"))

#combined = as.data.table(cbind(RT_norisk, scr_norisk))
#setnames(combined,c("RT_norisk","scr_norisk"))

#combined = as.data.table(cbind(RT_risk, scr_risk))
#setnames(combined,c("RT_risk","scr_risk"))

combined = as.data.table(cbind(RT, scr))
setnames(combined,c("RT","scr"))

#correlations scr risk scr norisk       
g <- ggscatter(combined, x = 'RT', y = 'scr' , 
               add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE,
               cor.coeff.args = list(label.x = 2000,label.y = 11.95, size = 8, label.sep = "\n"),
               cor.method = "spearman",
               xlab = 'RT m1000 m100ms TRAINED', ylab = "scr") 
g

p1<- ggpar(g,
           #ylim = c(.16, .24),
           font.ytickslab = 30,
           font.xtickslab = 27,
           font.main = 25,
           font.submain = 25,
           font.x = 27,
           font.y = 27)
p1

ggsave(filename = paste0(out_path, 'trained','_scr', '_RT','.png'), g, width =  5, height = 5)
