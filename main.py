import mne
import os
import os.path as op
import numpy as np
from functions import  make_scr_signal


data_path = '/net/server/data/Archive/prob_learn/pultsinak/adolescence/raw_maxfilter_with_events/'
subjects = ["P701", "P702", "P703", "P704", "P705", "P706", "P707", "P708", "P709",
          "P710", "P711", "P712", "P713", "P715", "P714","P716", "P717", "P718",
          "P719", "P720", "P721", "P722", "P723", "P724", "P725", "P726", "P727",
          "P728", "P729", "P730"]

period_start = -1.0
period_end = 3.0

rounds = [1, 2, 3, 4, 5]
#rounds = [2]

trial_type = ['norisk', 'risk', 'postrisk']

#feedback = ['positive', 'negative']

prefix_events = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/Events_teens_adults_trained/'
prefix = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/epochs_scr_trained/'
##############################################################################################################

for subj in subjects:
    for r in rounds:
        for cond in trial_type:
            try:
                epochs_scr = make_scr_signal(prefix_events, subj, r, cond, data_path, period_start, period_end)
                if epochs_scr != 0:
                    epochs_scr.save('{0}/{1}_run{2}_{3}_epo_scr.fif'.format(prefix, subj, r, cond), overwrite=True)
                    print('Saved!', '{0}/{1}_run{2}_{3}_epo_scr.fif'.format(prefix, subj, r, cond))
            except (OSError):
                    print('This file not exist')
