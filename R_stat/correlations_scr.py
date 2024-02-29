# LMM analysis separately on trained and not_trained or both.
# Choose mode to operate the subset you want  

#rm(list = ls())
library(reshape2)
library(data.table)
library(ggplot2)
library("ggpubr")
library(dplyr)

path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"
out_path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"

filename <- paste0(path,'df_LMEMRTtimetrainedpostrisk.csv')
df_large = read.csv(filename)
#trained only
df_large$train = 'trained'
df_large <- as.data.table(df_large)

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
colname

df_large <- as.data.table(df_large)
df_large[,(colname):=rowMeans(df_large[,.SD,.SDcol=cols])]

scr_norisk <- df_large[trial_type == 'norisk']$resp_scr..1....0.9._scr..0.2..0.1
scr_risk <- df_large[trial_type == 'risk']$resp_scr..1....0.9._scr..0.2..0.1
scr_postrisk <- df_large[trial_type == 'postrisk']$resp_scr..1....0.9._scr..0.2..0.1

combined = as.data.table(cbind(scr_norisk, scr_risk))
setnames(combined,c("scr_norisk","scr_risk"))

#correlations scr risk scr norisk       
g <- ggscatter(combined, x = 'scr_norisk', y = 'scr_risk' , 
               add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE,
               cor.coeff.args = list(label.x = 0.2,label.y = 9.7, size = 8, label.sep = "\n"),
               cor.method = "spearman",
               xlab = 'SCR HP TRAINED', ylab = "SCR LP") 
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
ggsave(filename = paste0(out_path, 'trained','_scr', '_corr','.png'), g, width =  5, height = 5)
