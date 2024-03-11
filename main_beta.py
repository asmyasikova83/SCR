import mne
import os
import os.path as op
import numpy as np
from functions import make_beta_signal
#from config import *

L_freq = 16
H_freq = 31
f_step = 2

time_bandwidth = 4 #(by default = 4)
# if delta (1 - 4 Hz) 
#n_cycles = np.array([1, 1, 1, 2]) # уточнить

#for others
freqs = np.arange(L_freq, H_freq, f_step)
n_cycles = freqs//2

period_start = -1.050
period_end = 3.050

baseline = (-0.35, -0.05)

freq_range = 'beta_16_30'

description = 'И для данных и для бейзлайн логарифмирование проводим на самых ранных этапах - сразу после суммирования по частотам.'
rounds = [1, 2, 3, 4, 5]
#rounds = [3]
trial_type = ['norisk', 'risk', 'postrisk']
#trial_type = ['norisk']
feedback = ['positive', 'negative']

#data_path = '/home/asmyasnikova83/SCR/links/'
data_path = '/net/server/data/Archive/prob_learn/pultsinak/adolescence/raw_maxfilter_with_events/'

prefix_events = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/Events_teens_adults_trained'
prefix_data = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR'
prefix_fix_cross = '/net/server/data/Archive/prob_learn/asmyasnikova83/SCR/fix_cross_mio_corr'

os.makedirs('/{0}/{1}_beta'.format(prefix_data, freq_range), exist_ok = True)
os.makedirs('/{0}/{1}_beta/{1}_epo'.format(prefix_data, freq_range), exist_ok = True)
########################## Обязательно делать файл, в котором будет показано какие параметры были заданы, иначе проверить вводные никак нельзя, а это необходимо при возникновении некоторых вопросов ############################################

lines = ["freq_range = {}".format(freq_range), description, "L_freq = {}".format(L_freq), "H_freq = {}, в питоне последнее число не учитывается, т.е. по факту частота (H_freq -1) ".format(H_freq), "f_step = {}".format(f_step), "time_bandwidth = {}".format(time_bandwidth), "period_start = {}".format(period_start), "period_end = {}".format(period_end), "baseline = {}".format(baseline)]

config_name = "/{0}/{1}_beta/{1}_epo/config.txt".format(prefix_data, freq_range)

with open(config_name, "w") as file:
    for  line in lines:
        file.write(line + '\n')


##############################################################################################################

#subjects = ["P716", "P717", "P718",
#            "P719", "P720", "P721", "P722", "P723", "P724", "P725", "P726", "P727",
#            "P728", "P729", "P730"]
subjects = ["P701", "P702", "P703", "P704", "P705", "P706", "P707", "P708", "P709",
            "P710", "P711", "P712", "P713", "P715", "P714","P716", "P717", "P718",
            "P719", "P720", "P721", "P722", "P723", "P724", "P725", "P726", "P727",
            "P728", "P729", "P730"]
for subj in subjects:
    for r in rounds:
        for cond in trial_type:
            try:
                epochs_tfr = make_beta_signal(prefix_events, prefix_fix_cross, subj, r, cond, data_path, period_start, period_end, L_freq, H_freq, f_step, baseline, n_cycles, time_bandwidth = 4)
                if epochs_tfr != 0:
                    epochs_tfr.save('/{0}/{1}_beta/{1}_epo/{2}_run{3}_{4}_epo.fif'.format(prefix_data, freq_range, subj, r, cond), overwrite=True)
            except (OSError):
                print('This file not exist')



