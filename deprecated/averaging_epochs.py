import mne
import os
import os.path as op
import numpy as np
import pandas as pd


subjects = ["P701", "P702", "P705", "P706","P707", "P708", "P709", "P710",
          "P711", "P712", "P713", "P714", "P716", "P718",
          "P721", "P722" ,"P723", "P724", "P725", "P726",
          "P727", "P728", "P729", "P730"]
rounds = [1, 2, 3, 4, 5]

#trial_type = ['norisk', 'prerisk', 'risk', 'postrisk']
trial_type = ['norisk', 'risk']


prefix = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/epochs_scr_trained/'
prefix_ave = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/averaged_epochs_trained/'

os.makedirs(prefix_ave, exist_ok = True)

timing  = 4001
for subj in subjects:
    for t in trial_type:
        ################################ epos ################################
        epos = np.empty((0, 1, timing))
        for r in rounds:
            try:
                epo = mne.read_epochs('{0}{1}_run{2}_{3}_epo_scr.fif'.format(prefix, subj, r, t), preload = True)             
                print(epo)
                epos = np.vstack([epos, epo.get_data()])
               
                
            except (OSError):
                
                print('This file not exist')
        print('epos size', epos.shape)
        ###### Шаг 1. Усреднили все epochs внутри испытуемого (между блоками 1 -6) #################
        if epos.size != 0:
            epos_mean = epos.mean(axis = 0)
            epos_mean = epos_mean.reshape(1,  timing) # добавляем ось для save
            info = mne.create_info(ch_names=['MISC001'], sfreq=1000, ch_types='misc', verbose=None)
            evoked  = mne.EvokedArray(epos_mean, info, tmin=-1.0, kind='average', baseline=None, verbose=None)
            print('evoked', evoked)
            # сохраняем данные, усредненные внутри испытуемого. 
            evoked.save('{0}{1}_{2}_evoked_ave.fif'.format(prefix_ave, subj, t))
            print('Saved!')
        else:
            print('No data')


