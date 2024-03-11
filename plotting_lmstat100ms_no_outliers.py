### For more details see plot_timecourses Jupyter Notebook  ######

import mne
import os.path as op
from matplotlib import pyplot as plt
import numpy as np
#import copy
import pandas as pd
from scipy import stats
import scipy

#subjects = ["P701", "P702", "P705", "P706", "P707", "P708", "P709", "P710", "P711", "P712",
#                             "P714", "P716",  "P718",  "P723", "P726", "P727", "P728", "P729"]
group = 'adols'

trial_type = ['norisk', 'postrisk']

#df_norisk = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/df_adols_no_outliers_2sd_norisk_trained.csv')
#df_risk = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/df_adols_no_outliers_2sd_risk_trained.csv')
df_norisk = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/df_adols_Z_norisk_trained_ed.csv')
df_risk = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/df_adols_Z_postrisk_trained_ed.csv')
#df_norisk = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/df_adols_ZZ_norisk_trained_deriv.csv')
#df_risk = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/df_adols_ZZ_risk_trained_deriv.csv')
#retrieve data from minus 1.0 sec to 2.0 sec

#norisk = np.log10(df_norisk.iloc[:, 3:43])
#risk = np.log10(df_risk.iloc[:, 3:43])
norisk = df_norisk.iloc[:, 3:42]
risk = df_risk.iloc[:, 3:42]

norisk_mean = norisk.mean(axis = 0)
risk_mean = risk.mean(axis = 0)

print(norisk_mean)
#risk_stderr = scipy.stats.sem(risk, axis=0)
#norisk_stderr = scipy.stats.sem(norisk, axis=0)
#norisk_mean = norisk_mean[np.newaxis, :]
#risk_mean = risk_mean[np.newaxis, :]
time = np.arange(-10, 29)
p_mul = 0

comp1_label = 'HP'
comp2_label = 'post-LP'

time_len = 39 # количество временных отчетов
#################### stat #####################################

##################### timecourse ###############################

cond1 = trial_type[0]
cond2 = trial_type[1]
p_mul_min = -0.1
p_mul_max = 0.1
 
#p_mul_min = -0.1
#p_mul_max = 0.7
#p= pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/p_vals_factor_significance_Z_deriv_scr_adols_permut_risk_trained.csv')
#p = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/p_vals_factor_significance_Z_scr_adols_permut_risk_trained.csv')
p = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/p_vals_factor_significance_Z_scr_adols_permut_postrisk_trained.csv')
p_fdr = p.iloc[:, 3]
p_val = p.iloc[:, 2]

print(p_val)
fig, ax = plt.subplots()
ax.set_xlim(time[0], time[-1])
ax.set_ylim(p_mul_min, p_mul_max)
ax.set_xlabel('Time ms', fontsize=40)
ax.set_ylabel('SCR Z', fontsize=40)

title = '%s_vs_%s_%s' %(comp1_label, comp2_label, group)
ax.set_title(title, fontsize = 40)
ax.plot(time, risk_mean, color='r', linewidth=7, label=comp2_label)
ax.plot(time, norisk_mean, color='b', linewidth=7, label=comp1_label)

ax.fill_between(time, y1 = p_mul_min, y2 = p_mul_max, where = (p_fdr < 0.05), facecolor = 'magenta', alpha = 0.99, step = 'pre')
ax.fill_between(time, y1 = p_mul_min, y2 = p_mul_max, where = (p_val < 0.05), facecolor = 'g', alpha = 0.99, step = 'pre')

ax.tick_params(labelsize = 40)
ax.legend(loc='upper left', fontsize = 40)
fig.set_size_inches((20, 10), forward=False)
plt.savefig(f'/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/Z_postrisk_trained_{cond1}_vs_{cond2}_{group}_')
plt.close()
