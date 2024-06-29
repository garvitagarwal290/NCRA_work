import numpy as np
import math
import matplotlib.pyplot as plt

# Making simulated time series------------------------------------------------------

P = 500 #in ms
pulse_width = 50 # in ms
pulse_height= 1	
N = 1000000 #total samples in data
tau = 0.2 #ms

time_series = np.ones((N))*.2 #setting noise level to 0.2
num= np.ceil((N* tau)/P) #no. of pulses/cycles. Last cycle might be incomplete
w= int(np.round(pulse_width/tau)) #no. of samples in pulse width

for i in range(int(num)):
	index= int(np.round(P*i/tau)) #index of every cycle
	n=1000 #position of pulse begining in every cycle
	time_series[index + n: index + n + w] = 1
	
#saving the time_series data in binary file time_series.bin
f= open("/home/garvit/Downloads/Prof Yashwant/meeting 12 (24.7.21)/time_series.bin","wb")
#arr = bytearray(time_series)  #(2lines) alternate way for making binary file
#f.write(arr)
time_series.tofile(f)
f.close()

display_window= 15000 # in no.of samples
plt.plot(time_series[:display_window])
plt.show()



