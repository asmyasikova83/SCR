import mne
import os
import os.path as op
import numpy as np
import pandas as pd
from scipy import stats
import copy
#import statsmodels.stats.multitest as mul

def fixation_cross_events(data_path, prefix_events, subj, r):
    no_risk = np.loadtxt('{2}{0}_run{1}_norisk.txt'.format(subj, r, prefix_events), dtype='int')    
    raw = op.join(data_path, '{0}/{0}_run{1}_raw_tsss_mc_trans.fif'.format(subj, r))        
    # Load data

    raw_d = mne.io.Raw(raw, preload=True)

    events_raw = mne.find_events(raw_d, stim_channel='STI101', output='onset', 
                                     consecutive='increasing', min_duration=0, shortest_event=1, mask=None, 
                                     uint_cast=False, mask_type='and', initial_event=False, verbose=None)
    print('events_raw', events_raw)
    if no_risk.shape == (5,):
        no_risk = no_risk.reshape(1,5)

    if no_risk.shape != 0:
        idxs = []
        for i in range(len(events_raw)):
            for j in range(len(no_risk)):
                if events_raw[i][2] == 1:
                    if (events_raw[i][0] < no_risk[j][0]):
                        a = events_raw[i][0] - no_risk[j][0]
                        idxs.append(a)
        minn = np.sort(np.abs(idxs))
        minn = np.unique(minn)
        print('minn', minn[0:len(no_risk)])
        minn = minn[0:len(no_risk)]
        print('len minn', len(minn))
        full_ev = []
        for k in minn:
            for i in range(len(events_raw)):
                for j in range(len(no_risk)):
                    if events_raw[i][2] == 1:
                        if (events_raw[i][0] < no_risk[j][0]):
                            if (np.abs(events_raw[i][0] - no_risk[j][0]) == k):
                                #print('K', k)
                                #print('diff', np.abs(events_raw[i][0] - no_risk[j][0]))
                                full_ev.append(events_raw[i])
                                #print('full ev [i]', events_raw[i])
                                #print('no risk j [i]', no_risk[j])
        full_ev = full_ev[0:len(no_risk)]
        print('fix cross', full_ev)
        print('no risk', no_risk)
        assert(len(full_ev) == len(no_risk))
        print('lens equal')
        event_fixation_cross_norisk = np.array(full_ev)
    else:
        event_fixation_cross_norisk = []
    return(event_fixation_cross_norisk)

    #result = []
    #prev_item = [0, 0, 0]
    #for item in events:
    #    cur_item = item
    #    if prev_item[0] != cur_item[0]:
    #        result.append(item)
    #        prev_item = item
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

def make_beta_signal(prefix_events, prefix_fix_cross, subj, r, cond, data_path, period_start, period_end, L_freq, H_freq, f_step, baseline, n_cycles, time_bandwidth = 4):
    freqs = np.arange(L_freq, H_freq, f_step)

    events = np.loadtxt('{2}/{0}_run{1}_norisk_fix_cross.txt'.format(subj, r, prefix_fix_cross), dtype='int')
    if events.shape != (1,):
        events = np.matrix(events)
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводим shape к виду (N,3)
    if events.shape == (3,):
        events = events.reshape(1,3)
    events_response = np.loadtxt('{3}/{0}_run{1}_{2}.txt'.format(subj, r, cond, prefix_events), dtype='int')
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводи shape к виду (N,3)
    if events_response.shape == (5,):
        events_response = events_response.reshape(1,5)
    
    if events.shape != (0,) and events_response.shape != (0,):
        events_response =events_response[:,0:3]
    else:
        return 0
    if not events.any():
        return 0
    #data with restored events
    raw_fname = op.join(data_path, '{0}/run{1}_{0}_raw_ica.fif'.format(subj, r))
    
    raw_data = mne.io.Raw(raw_fname, preload=True)
    picks = mne.pick_types(raw_data.info, meg = True, eog = True, misc = True)
    
        #epochs for baseline
     # baseline = None, чтобы не вычитался дефолтный бейзлайн
    epochs = mne.Epochs(raw_data, events, event_id = None, tmin = -1.0, tmax = 1.0, baseline = None, picks = picks, preload = True)
    epochs.resample(300) #resample in accordance with ti+ime window of multitaper (500 ms)
    
        
    freq_show_baseline = mne.time_frequency.tfr_multitaper(epochs, freqs = freqs, n_cycles = n_cycles, time_bandwidth = time_bandwidth, use_fft = False, return_itc = False, average=False).crop(tmin=baseline[0], tmax=baseline[1], include_tmax=True) #frequency of baselin
    #add up all values according to the frequency axis
    b_line = freq_show_baseline.data.sum(axis=-2)
    
    # логарифмируем данные для бейзлайн и получаем бейзлайн в дБ
    b_line = 10*np.log10(b_line)
	    
	# Для бейзлайна меняем оси местами, на первом месте число каналов
    b_line = np.swapaxes(b_line, 0, 1)
        
    # выстраиваем в ряд бейзлайны для каждого из эвентов, как будто они происходили один за другим
    a, b, c = b_line.shape
    b_line_ave = b_line.reshape(a, b * c)
    # Усредняем бейзлайн по всем точкам, получаем одно число (которое будем вычитать из data для каждого канала) (то же самое если вместо смены осей усреднить сначала по времени, а потом по эпохам)
	                        
    b = b_line_ave.mean(axis=-1)
    b_line_new_shape = b[:, np.newaxis, np.newaxis]  
	####### ДЛЯ ДАННЫХ ##############
    
   
    epochs = mne.Epochs(raw_data, events_response, event_id = None, tmin = period_start, 
		                tmax = period_end, baseline = None, picks = picks, preload = True)
		       
    epochs.resample(300) 

    freq_show = mne.time_frequency.tfr_multitaper(epochs, freqs = freqs, n_cycles = n_cycles, time_bandwidth = time_bandwidth, use_fft = False, return_itc = False, average=False)
    temp = freq_show.data.sum(axis=2)
    # логарифмируем данные и получаем данные в дБ
    temp = 10*np.log10(temp)
	    
	####### Для данных так же меняем оси местами
    data = np.swapaxes(temp, 0, 1)
    data = np.swapaxes(data, 1, 2)
	
	
    #Вычитаем из данных в дБ бейзлайн в дБ
    data_dB = data - b_line_new_shape 
    
    # меняем оси обратно   
    data_dB = np.swapaxes(data_dB, 1, 2)
    data_dB = np.swapaxes(data_dB, 0, 1)
    
            
    freq_show.data = data_dB[:, :, np.newaxis, :]
        
    #freq_show.data = freq_show.data[:, :, np.newaxis, :]
        
    #33 is an arbitrary number. We have to set some frequency if we want to save the file
    freq_show.freqs = np.array([33])
        
    #getting rid of the frequency axis	
    freq_show.data = freq_show.data.mean(axis=2) 
        
    epochs_tfr = mne.EpochsArray(freq_show.data, freq_show.info, tmin = period_start, events = events_response)
        
    return (epochs_tfr)

def make_scr_bline_signal(prefix_events, prefix_fix_cross, subj, r, cond, data_path, period_start, period_end, L_freq, H_freq, f_step, baseline, n_cycles, time_bandwidth = 4):
    freqs = np.arange(L_freq, H_freq, f_step)

    events = np.loadtxt('{2}/{0}_run{1}_norisk_fix_cross.txt'.format(subj, r, prefix_fix_cross), dtype='int')
    if events.shape != (1,):
        events = np.matrix(events)
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводим shape к виду (N,3)
    if events.shape == (3,):
        events = events.reshape(1,3)
        #events = np.unique(events.reshape)
    events_response = np.loadtxt('{3}/{0}_run{1}_{2}.txt'.format(subj, r, cond, prefix_events), dtype='int')
    # если только одна метка, т.е. одна эпоха, то выдается ошибка, поэтому приводи shape к виду (N,3)
    if events_response.shape == (5,):
        events_response = events_response.reshape(1,5)
    
    if events.shape != (0,) and events_response.shape != (0,):
        events_response =events_response[:,0:3]
    else:
        return 0
    if not events.any():
        return 0
    #if not events_response.all():
    #    return 0
    #data with restored events
    raw_fname = op.join(data_path, '{0}/run{1}_{0}_raw_ica.fif'.format(subj, r))
    
    raw_data = mne.io.Raw(raw_fname, preload=True)
    #picks = mne.pick_types(raw_data.info, meg = True, eog = True, misc = True)
    picks=['MISC001']
    
        #epochs for baseline
     # baseline = None, чтобы не вычитался дефолтный бейзлайн
    epochs = mne.Epochs(raw_data, events, event_id = None, tmin = -1.0, tmax = 1.0, baseline = None, picks = picks, preload = True)
    epochs.resample(300) #resample in accordance with ti+ime window of multitaper (500 ms)
    
        
	####### ДЛЯ ДАННЫХ ##############
    
   
    epochs_data = mne.Epochs(raw_data, events_response, event_id = None, tmin = period_start, 
		                tmax = period_end, baseline = None, picks = picks, preload = True)
		       
    epochs_data.resample(300) 

    	
    temp = np.swapaxes(epochs_data.get_data(), 0, 1)
    data_bline = np.swapaxes(temp, 0, 1)
    temp = np.swapaxes(epochs.get_data(), 1, 2)
    data_data = np.swapaxes(temp, 1, 2)
    #Вычитаем из данных в дБ бейзлайн в дБ
    data_blined = data_bline - data_data  
    
    # меняем оси обратно   
    data_blined = np.swapaxes(data_blined, 1, 2)
    data_blined = np.swapaxes(data_dB, 0, 1)
    
            
    data_blined = data_blined[:, :, np.newaxis, :]
        
    info = mne.create_info(ch_names=['MISC001'], sfreq=300, ch_types='misc', verbose=None)
    epochs_scr = mne.EpochsArray(data_blined, info, tmin = period_start, events = events_response)
    print('epochs_scr', epochs_scr)
    exit()
    return (epochs_scr)

    #SCR channels picks=['MISC001']
	####### ДЛЯ ДАННЫХ ##############
    #if events_response.shape != (0,):
    #    #we do nott need real time and response time
    #    events_response = events_response[:,0:3]
    #    # baseline = None, чтобы не вычитался дефолтный бейзлайн
    #    epochs = mne.Epochs(raw_data, events_response, event_id = None, tmin = period_start, 
	#	                tmax = period_end, baseline = None, picks = ['MISC001'], detrend = None, preload = True, reject_by_annotation = True,  on_missing='ignore')

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
    
