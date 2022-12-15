import binascii
import numpy as np
import math as m
import matplotlib.pyplot as plt
import csv
import pandas as pd

# command for start receiver vy Hercules AA 55 AA 55 07 03

filename = '2_herc_withoutGen(2).log'
with open(filename, 'rb') as file:
    
    content = file.read()
    
arr = bytearray(content).hex()
file.close()
i_quadr,q_quadr = [],[]
arr_string_q, arr_string_i = [], []
for i in range(8, len(arr)-15, 16):

    arr_string_q.append(''.join(arr[i+6:i+8]) +''.join(arr[i+4:i+6]) +''.join(arr[i+2:i+4]) + ''.join(arr[i:i+2]))
    arr_string_i.append(''.join(arr[i+14:i+16])+''.join(arr[i+12:i+14])+''.join(arr[i+10:i+12]) + ''.join(arr[i+8:i+10]))
    if (arr[i+4:i+6]=='ff'):
        string_q = ''.join(arr[i+2:i+4]) + ''.join(arr[i:i+2])
        dig_q = ~(int(string_q, 16))&(65535)
        q_quadr.append(-dig_q/(8388608 *2))    
    elif (arr[i+4:i+6]=='fe'): 
        string_q = ''.join(arr[i+4:i+6]) +''.join(arr[i+2:i+4]) + ''.join(arr[i:i+2])
        dig_q = (~int(string_q, 16)&131071)
        q_quadr.append(-dig_q/(8388608 *2))   
    else:
        string_q = ''.join(arr[i+4:i+6]) +''.join(arr[i+2:i+4]) + ''.join(arr[i:i+2])
        q_quadr.append(int(string_q, 16)/(8388608 *2))
    if (arr[i+12:i+14]=='ff'):
        string_i = ''.join(arr[i+10:i+12]) + ''.join(arr[i+8:i+10])
        dig_i = ~(int(string_i, 16))&(65535)
        i_quadr.append(-dig_i/(8388608 *2)) 
    elif (arr[i+12:i+14]=='fe'):
        string_i = ''.join(arr[i+12:i+14])+''.join(arr[i+10:i+12]) + ''.join(arr[i+8:i+10])
        dig_i = (~int(string_i, 16)&131071)
        i_quadr.append(-dig_i/(8388608 *2)) 
    else: 
        string_i = ''.join(arr[i+12:i+14])+''.join(arr[i+10:i+12]) + ''.join(arr[i+8:i+10])
        i_quadr.append(int(string_i, 16)/(8388608 *2))


pd_data = pd.DataFrame({"I":i_quadr,"Q":q_quadr, "I_hex":arr_string_i, "Q_hex":arr_string_q})
print(pd_data[4290:4310])
pd_data.to_csv('Quadratures_panda_rez.csv', encoding='cp1251', sep=';')

power = []
sig_fft = np.array(np.random.random(8192), dtype=np.complex64)

for j in range(len(q_quadr)):
    power.append(m.sqrt(q_quadr[j]**2+i_quadr[j]**2))
    sig_fft[j] = complex(i_quadr[j],-(q_quadr[j]))

   
Fs = 2.7e6 # lets say we sampled at 2.7 MHz
# assume x contains your array of IQ samples
N = 8192
power = power[0:N] # we will only take the FFT of the first 1024 samples, see text below

PSD = (np.abs(np.fft.fft(power))/N)**2
print("Масимальная мощность выборки {log_max_power} дБ\n\
индекс максимального элемента {index_max}"\
          .format(log_max_power  = 20.0*np.log10(max(power)), index_max =max(enumerate(power),key=lambda x: x[1])[0]))
PSD_log = 10.0*np.log10(PSD)
PSD_shifted = np.fft.fftshift(PSD_log)
power = power * np.hamming(len(power)) # apply a Hamming window
center_freq = 15e6 # frequency we tuned our SDR to
f = np.arange(Fs/-2.0, Fs/2.0, Fs/N) # start, stop, step.  centered around 0 Hz
f += center_freq # now add center frequency

markers_on = float(max(PSD_shifted.round(6)))
max_index = max(enumerate(PSD_shifted),key=lambda x: x[1])[0]

fig,geeks = plt.subplots()

plt.figure(1)
geeks.plot(f, PSD_shifted,'-g')
geeks.annotate('(%.1f)'%(markers_on), xy =(round(f[max_index], 2), markers_on),
                xytext =(f[max_index]+0.05e6, markers_on), 
                )
geeks.plot(f[max_index],markers_on, marker='o')
geeks.set_xlim(14e6,16e6)
plt.grid()
plt.figure(2)
plt.plot(range(len(q_quadr)), q_quadr,label='q_quadr')
plt.legend()
plt.figure(2)
plt.plot(range(len(i_quadr)), i_quadr, label='i_quadr')
plt.legend()

plt.grid()
plt.show()
