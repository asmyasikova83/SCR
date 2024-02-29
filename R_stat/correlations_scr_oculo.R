    # LMM analysis separately on trained and not_trained or both.
    # Choose mode to operate the subset you want  
    
    #rm(list = ls())
    
    # LMM analysis RT
  
    #library(reshape2)
    library(data.table)
    #library(ggplot2)
    #library("ggpubr")
    library(tidyr)
    library(dplyr)
    
    subjects = c("P701", "P702", "P703", "P704", "P705", "P706", "P707", "P708", "P709",
                 "P710", "P711", "P712", "P713", "P715", "P714","P716", "P717", "P718",
                 "P719", "P720", "P721", "P722", "P723", "P724", "P725", "P726", "P727",
                 "P728", "P729", "P730")
    
    path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/Oculo/TEENS/"
    out_path <- paste0(path,'pics/')
    
    #eyetracker
    rt_teens <- read.csv(paste0(path,"resp_merged_teens_big_extended.csv"))
    rt_teens <- as.data.table(rt_teens)
    rt_teens$group <- 'teens'
    
    df_large_group <- rt_teens
    
    #2
    #RT filter (for all intervals)
    df_large_group <- as.data.table(df_large_group)
    df_large_group <- df_large_group[RT>300]
    
    #standardazitng fname
    df_large_group$fname <- gsub ("_1", "", df_large_group$fname)
    df_large_group$fname <- gsub ("p0", "P0", df_large_group$fname)
    df_large_group$trial_type <-gsub ("_", "", df_large_group$trial_type) 
    
    colnames(df_large_group)
    #1
    #average over time interval 
    #decision making
    int_beg <- 'm999'
    int_end <- 'm100'
    
    ncol_min <- grep(paste0('_', int_beg, '_'), colnames(df_large_group))[1]
    ncol_max <- grep(paste0('_', int_end, ''), colnames(df_large_group))[1]
    col <- c(colnames(df_large_group)[ncol_min:ncol_max])
    
    colname <- paste0('resp_',int_beg,'_',int_end)
    colname
    
    df_large_group[,(colname):=rowMeans(df_large_group[,.SD,.SDcol=col], na.rm=TRUE)]
    
    #Z for pupil
    
    ncol_min <- grep('V_m999_m950',colnames(df_large_group))
    ncol_max <- grep('V_2951_3000',colnames(df_large_group))
    col <- c(colnames(df_large_group)[ncol_min:ncol_max])
    
    df_filt_m <- data.table(melt(df_large_group[blink==F], id.vars = c('fname','time'), measure.vars = col))
    Zmeans <- df_filt_m[, mean(value), by=fname]
    setnames(Zmeans,'V1','Zmean_2')
    Zsds <- df_filt_m[, sd(value), by=fname]
    setnames(Zsds,'V1','Zsd_2')
    
    df_large_group <- merge(df_large_group,Zmeans,by='fname',all.x = TRUE)
    df_large_group <- merge(df_large_group,Zsds,by='fname',all.x = TRUE)
    
    #apply Z for pupil
    cols <- colnames(df_large_group)[grep('resp_|V_|BL_|inter_trial',colnames(df_large_group))]
    for (j in cols) set(df_large_group, j = j, value = (df_large_group[[j]]-df_large_group$Zmean_2)/df_large_group$Zsd_2)
    
    df_large_group$trial_type4 <- df_large_group$trial_type
    df_large_group[prev_risk==0 & risk==1 & next_risk==0, trial_type4:= 'risk']
    df_large_group[prev_risk==0 & risk==0 & next_risk==0, trial_type4:= 'norisk']
    df_large_group[prev_risk==1 & risk==0 & next_risk==0, trial_type4:= 'postrisk']
    
    #3
    df_large_group<- df_large_group[blink==F]
    
    #settings for train untrain
    df_large_group[trained == FALSE, train:='learning'][trained== TRUE, train:='trained']
    
    df_large_group <- df_large_group[trial_type4 %in% c('norisk','risk','postrisk' )]
    df_large_group = df_large_group[train=='trained']
    colnames(df_large_group)[colnames(df_large_group) == 'fname'] <- 'subject'
    
    
    ############################scr################################################
    path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/"
    out_path <- "C:/Users/trosh/OneDrive/jobs_Miasnikova/SCR/pics/"
    filename <- paste0(path,'df_LMEMRTtimetrainedpostrisk.csv')
    
    df_large = read.csv(filename)
    df_large$train = 'trained'
    df_large = dftrained
    df_large <- as.data.table(df_large)
    colnames(df_large)[colnames(df_large) == 'trial_type'] <- 'trial_type4'
    colnames(df_large)[colnames(df_large) == 'response_time'] <- 'RT'
    
    #RT filter
    df_large <- df_large[RT>300]
    #1
    #average over time interval
    
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
    df_large = df_large[train=='trained']
    #colname

    #oculo
    q<- df_large_group %>%
      group_by(subject, train, trial_type4)%>%summarise(mean_pupil=mean(resp_m999_m100))
    
    q <- as.data.table(unique(q))
    
    #SCR
    qq<- df_large %>%
      group_by(subject, train, trial_type4)%>%summarise(mean_scr=mean(resp_scr..1....0.9._scr..0.2..0.1.),
                                                        mean_RT=mean(RT))
    
    qq <- as.data.table(unique(qq))
    
    data = merge(q, qq, by = c('subject', 'trial_type4'),all.x = TRUE)
    data = na.omit(data)

    data = data[trial_type4=='risk']
    

#    m <- aov(data$mean_scr ~ mean_pupil*mean_RT, data=data)
#    summary(m)

#correlations scr risk oculo norisk       
g <- ggscatter(data, x = 'mean_pupil', y = 'mean_scr' , 
               add = "reg.line", conf.int = TRUE,
               cor.coef = TRUE,
               cor.coeff.args = list(label.x = -0.25,label.y = 0.62, size = 8, label.sep = "\n"),
               cor.method = "spearman",
               xlab = 'pupil -1000 ms 0ms \n RISK ', ylab = "scr -1000 ms 0ms") 
g

p1<- ggpar(g,
           #ylim = c(-.3, .3),
           font.ytickslab = 30,
           font.xtickslab = 27,
           font.main = 25,
           font.submain = 25,
           font.x = 27,
           font.y = 27)
p1

#ggsave(filename = paste0(out_path, 'trained','_scr', '_RT','.png'), g, width =  5, height = 5)
