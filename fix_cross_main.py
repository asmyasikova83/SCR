import mne
import os
import os.path as op
import numpy as np
from functions import fixation_cross_events

data_path = '/home/asmyasnikova83/SCR/links/'
#data_path = '/net/server/data/Archive/prob_learn/pultsinak/adolescence/raw_maxfilter_with_events/'
subjects = ["P701", "P702", "P703", "P704", "P705", "P706", "P707", "P708", "P709",
         "P710", "P711", "P712", "P713", "P715", "P714","P716", "P717", "P718",
         "P719", "P720", "P721", "P722", "P723", "P724", "P725", "P726", "P727",
         "P728", "P729", "P730"]
#subjects = ["P726"]
rounds = [1, 2, 3, 4, 5]
#rounds = [1]

trial_type = ['norisk']
prefix_events = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/Events_teens_adults_trained/'
prefix_data = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/'
os.makedirs("{0}/fix_cross_mio_corr/".format(prefix_data), exist_ok = True)

for subj in subjects:
    for r in rounds:
    
        for t in trial_type:
            try:
                event_fixation_cross_norisk = fixation_cross_events(data_path, prefix_events, subj, r)
                print('{2}/fix_cross_mio_corr/{0}_run{1}_norisk_fix_cross.txt'.format(subj, r, prefix_data))
                print('event_fixation_cross_norisk', event_fixation_cross_norisk)
                np.savetxt("{2}/fix_cross_mio_corr/{0}_run{1}_norisk_fix_cross.txt".format(subj, r, prefix_data), event_fixation_cross_norisk, fmt="%s")
                print('Saved!')
            except OSError:
                print('This file not exist')
