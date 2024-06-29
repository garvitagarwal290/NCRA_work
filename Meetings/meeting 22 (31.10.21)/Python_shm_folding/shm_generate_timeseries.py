import numpy as np
import math
import matplotlib.pyplot as plt
from multiprocessing import shared_memory, resource_tracker

# Generating simulated time series------------------------------------------

P = 500.12345 #in ms
pulse_width = 50.54321 # in ms
pulse_height= 1	
N = 1000000 #total samples in data
tau = 0.234 #ms

shm=shared_memory.SharedMemory(name='shm_folding',create=True, size=N*4)
time_series=np.ndarray(shape=(N,),dtype=np.float32,buffer=shm.buf)

time_series+=0.2
num= np.ceil((N* tau)/P) #no. of pulses/cycles. Last cycle might be incomplete
w= round(pulse_width/tau) #no. of samples in pulse width

n=1000 #position of pulse begining (in terms of no. of samples) in every cycle

for i in range(int(num)):
	index= int(np.round(P*i/tau)) #index of every cycle
	
	time_series[index + n: index + n + w] = 1

display_window= 15000 # in no.of samples

shm.close()




