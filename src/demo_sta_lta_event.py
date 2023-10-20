import obspy
from obspy.signal.trigger import classic_sta_lta, plot_trigger, trigger_onset
import matplotlib.pyplot as plt
import numpy as np

file_path = './'

windows = 5    # length of cut window
THRESH_ON = 8   # STA/LTA trigger on (generally, maximum is 80% of the L/S)
THRESH_OFF = 5  # trigger off
SHORT_WINDOW = 0.1
LONG_WINDOW = 1

file = file_path + 'test_data/seis.sac'
input = obspy.read(file)[0]
sampling_interval = input.stats.sac.delta
sampling_frequency = int(1/sampling_interval)
total_point = input.stats.npts  # total points of input data
total_time = total_point*sampling_interval  # time length of input data
total_t = np.arange(0, total_time, sampling_interval)

plt.figure(figsize=(12, 6))
plt.subplot(211)
plt.plot(total_t, input.data, 'k', linewidth=1)
plt.ylabel('Amplitude')

# obspy.signal.trigger.classic_sta_lta
cft = classic_sta_lta(input.data, SHORT_WINDOW/sampling_interval, LONG_WINDOW/sampling_interval)
plt.subplot(212)
plt.plot(total_t, cft, 'y', linewidth=1)

cft1 = trigger_onset(cft, THRESH_ON, THRESH_OFF)    # STA/LTA trigger
for i in range(len(cft1)):
    plt.vlines(cft1[i, 0]*sampling_interval, min(cft), max(cft), colors='r', linestyles='dashed')
    plt.vlines(cft1[i, 1]*sampling_interval, min(cft), max(cft), colors='b', linestyles='dashed')
plt.xlabel('Time (s)')
plt.ylabel('STA/LTA')
plt.savefig('sla_lta.png', dpi=1000, bbox_inches='tight')
 
for i in range(len(cft1)):
    tmp = (cft1[i, 1] - cft1[i, 0])*sampling_interval  # time length for each event
    
    # Take the sum of the time before 1/3*windows of the event and the time after 2/3*windows of the event
    tmp = int((windows - tmp)/3)
    # Cut data
    input_copy = input.copy()
    input_copy.trim(input.stats.starttime + int(cft1[i, 0]*sampling_interval) - tmp, \
                    input.stats.starttime + int(cft1[i, 0]*sampling_interval) - tmp + windows)
    # if (len(input_copy.data) % sampling_frequency) != 0:
    #     input_copy.data = input_copy.data[: -1]
    # if len(input_copy.data) != sampling_frequency*windows:
    #     temp_endtime = input_copy.stats.endtime
    #     input_copy = input.copy()
    #     input_copy.trim(temp_endtime - windows, temp_endtime)
    #     print(len(input_copy.data))
    # input_copy.data = input_copy.data[: -1]
    
    output_path = file_path + 'test_data/seis_new_' + str(i+1) + '.sac'
    input_copy.write(output_path, format='SAC')
print(file + ' has been cut.')