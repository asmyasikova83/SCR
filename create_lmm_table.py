
import mne
import os
import os.path as op
import numpy as np
import pandas as pd
from mne import set_log_level
from functions import make_subjects_df


subjects = ["P701", "P702", "P703", "P704", "P705", "P706", "P707", "P708", "P709",
          "P710", "P711", "P712", "P713", "P715", "P714","P716", "P717", "P718",
          "P719", "P720", "P721", "P722", "P723", "P724", "P725", "P726", "P727",
          "P728", "P729", "P730"] 

rounds = [1, 2, 3, 4, 5]

trial_type = ['norisk', 'risk', 'postrisk']
#feedback = ['positive', 'negative']

# interval of interest ( +/- 100 ms)
tmin = -1.0
tmax = 3.0
step = 0.1

prefix_data = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/'
os.makedirs('{0}dataframe_for_LMEM_trained'.format(prefix_data), exist_ok = True)

prefix_events = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/Events_teens_adults_trained/'
#####################################################################################################
df = pd.DataFrame()
for subj in subjects:
    for r in rounds:
        for t in trial_type:
            try:
                epoch = mne.read_epochs('{0}epochs_scr_trained/{1}_run{2}_{3}_epo_scr.fif'.format(prefix_data, subj, r, t), preload = True)
                                        
                df_subj = make_subjects_df(prefix_events, epoch, subj, r, t, tmin, tmax, step)
                df = df.append(df_subj)            
            except (OSError, FileNotFoundError):
                print('This file not exist')
df.to_csv('{0}dataframe_for_LMEM_trained/df_LMEMRTtimetrainedpostrisk.csv'.format(prefix_data))
                    
	
	
		

