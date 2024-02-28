### For more details see plot_timecourses Jupyter Notebook  ######

import mne
import os.path as op
from matplotlib import pyplot as plt
import numpy as np
import copy
import pandas as pd
from scipy import stats


subjects = ["P701", "P702", "P705", "P706", "P707", "P708", "P709", "P710", "P711",
                             "P714", "P716",  "P718",  "P723", "P726", "P727", "P728", "P729"]
#subjects = ["P001", "P002", "P003", "P004", "P005",
#                                    "P006", "P007", "P008", "P009", "P010",
#                                   "P011", "P012", "P043", "P044", "P045",
#                                    "P046", "P047", "P048", "P050"]
group = 'adols'
time_len = 4001 # количество временных отчетов

trial_type = ['norisk', 'risk']
contr = np.zeros((len(subjects), 2, 1, time_len))

prefix = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/averaged_epochs_trained/'
for ind, subj in enumerate(subjects):
    evoked1= mne.Evoked('{0}{1}_{2}_evoked_ave.fif'.format(prefix, subj, trial_type[0]))
    print(evoked1)
    evoked2= mne.Evoked('{0}{1}_{2}_evoked_ave.fif'.format(prefix, subj, trial_type[1]))
    contr[ind, 0, :, :] = evoked1.data
    contr[ind, 1, :, :] = evoked2.data
comp1 = contr[:, 0, :, :]
comp2 = contr[:, 1, :, :]

#################### stat #####################################
#comp1 = np.log10(comp1)
#comp2 = np.log10(comp2)
comp1_mean = comp1.mean(axis=0).mean(axis=0)
comp2_mean = comp2.mean(axis=0).mean(axis=0)

#averaging over 100 ms intervals
sums_norisk = np.zeros((40))
sums_risk = np.zeros((40))
for i in range(40):
    for j in range(100):
        sums_norisk[i] += comp1_mean[i*100 + j]
        sums_risk[i] += comp2_mean[i*100 + j]
    sums_norisk[i] /= 100.0
p = pd.read_csv('/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/p_vals_factor_significance_scr_adols_permut_trained100ms.csv')
p_fdr = p.iloc[:, 3]
p_val = p.iloc[:, 2]

##################### timecourse ###############################

cond1 = trial_type[0]
cond2 = trial_type[1]
#set axes lims
#p_mul_min = -0.9
#p_mul_max = -0.8
 

time = np.arange(-1.0,  3.00, 0.1)
title = f'SCR in {group}  {cond1} vs {cond2} trained'
plt.figure()
plt.rcParams['axes.facecolor'] = 'none'
plt.xlim(time[0], time[-1])
#plt.ylim(p_mul_min, p_mul_max)
plt.plot([0, 0.001], [-50, 50], color='k', linewidth=3, linestyle='--', zorder=1)
plt.plot([-50, 50], [0, 0.001], color='k', linewidth=3, linestyle='--', zorder=1)
#FB axis
plt.plot([1, 1.001], [-50, 50], color='k', linewidth=3, linestyle='--', zorder=1)
plt.plot(time, sums_norisk, color='b', linewidth=3, label=cond1)
plt.plot(time, sums_risk, color='r', linewidth=3, label=cond2)
#stat
#plt.fill_between(time, y1 = p_mul_min, y2 = p_mul_max, where = (p_fdr < 0.05), facecolor = 'm', alpha = 0.46, step = 'pre')

#plt.fill_between(time, y1 = p_mul_min, y2 = p_mul_max, where = (p_val < 0.05), facecolor = 'g', alpha = 0.46, step = 'pre')
plt.tick_params(labelsize = 20)
plt.legend(loc='upper right', fontsize = 14)
plt.title(title, fontsize = 20)
plt.savefig(f'/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/av_{cond1}_vs_{cond2}_{group}_trained')
plt.close()
