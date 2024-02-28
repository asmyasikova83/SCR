import mne
import os
import os.path as op
import numpy as np
import pandas as pd
from scipy import stats
import copy
#import statsmodels.stats.multitest as mul


# функция получения эпох
def make_scr_signal(prefix_events, subj, r, cond, data_path, period_start, period_end):
    #read events
    #events_response = np.loadtxt('{4}/{0}_run{1}_{2}_fb_cur_{3}.txt'.format(subj, r, cond, fb, prefix_events), dtype='int')
    events_response = np.loadtxt('{3}{0}_run{1}_{2}.txt'.format(subj, r, cond, prefix_events), dtype='int')
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводи shape к виду (N,3)
    if events_response.shape == (5,):
        events_response = events_response.reshape(1,5)
   
    #data with restored events
    raw_fname = op.join(data_path, '{0}/run{1}_{0}_raw_ica.fif'.format(subj, r))
    
    raw_data = mne.io.Raw(raw_fname, preload=True)
    #SCR channels picks=['MISC001']
	####### ДЛЯ ДАННЫХ ##############
    if events_response.shape != (0,):
        #we do nott need real time and response time
        events_response = events_response[:,0:3]
        # baseline = None, чтобы не вычитался дефолтный бейзлайн
        epochs = mne.Epochs(raw_data, events_response, event_id = None, tmin = period_start, 
		                tmax = period_end, baseline = None, picks = ['MISC001'], detrend = None, preload = True, reject_by_annotation = True,  on_missing='ignore')
    else:
        print('epochs not found')
        epochs = 0

    return (epochs)   

def make_subjects_df(prefix_events, epoch, subj, r, t, tmin, tmax, step):
    import re

    time_intervals = np.arange(tmin, tmax, step)
    list_of_time_intervals = []
    i = 0
    while i < (len(time_intervals) - 1):
        interval = time_intervals[i:i+2]
        list_of_time_intervals.append(interval)
        #print(interval)
        i = i+1
    list_of_scr = []    
    for i in list_of_time_intervals:
        epoch_in_interval = epoch.copy()
        epoch_in_interval = epoch_in_interval.crop(tmin=i[0], tmax=i[1], include_tmax=True)

        mean_epoch = epoch_in_interval.get_data().mean(axis=-1)
        scr = []

        for j in range(len(mean_epoch)):
            a = mean_epoch[j]
            a = float(a)
            scr.append(a)
        list_of_scr.append(scr)
    trial_number = range(len(mean_epoch))

    events_response = np.loadtxt('{3}{0}_run{1}_{2}.txt'.format(subj, r, t, prefix_events), dtype='int')
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводи shape к виду (N,3)
    if events_response.shape == (5,):
        events_response = events_response.reshape(1,5)

    print('events_response', events_response)
    if len(events_response) != 0:
        real_events = events_response[:,3]
        #timing for further Z normalization 
        real_time = []
        for j in range(len(real_events)):
            r = real_events[j]
            real_time.append(r)

        events = events_response[:,4]
        
        #response time
        rt = []
        for j in range(len(events)):
            r = events[j]
            rt.append(r)

    subject = [subj]*len(mean_epoch)
    run = ['run{0}'.format(r)]*len(mean_epoch)
    trial_type = [t]*len(mean_epoch)
    
    df = pd.DataFrame()
    
    df['trial_number'] = trial_number
    
    # scr на интервалах
    for idx, scr_ in enumerate(list_of_scr):
        df['scr %s'%list_of_time_intervals[idx]] = scr_
    
    
    df['scr'] = scr
    df['subject'] = subject
    df['round'] = run
    df['trial_type'] = trial_type
    df['real_time'] = real_time
    df['response_time'] = rt
    
    #check
    print(df)
    return (df)
    
