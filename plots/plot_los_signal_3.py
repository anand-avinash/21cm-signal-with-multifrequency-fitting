#To plot los signal for a given los, tau and polystr.
#Plots eor_p_fg, eor_p_noise, eor_p_noise_fg in different panels for all offsets

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import matplotlib as mpl


######################################################
#### Choose (x,y) for plotting los signal
los_array=[0,145]#[124,12]#[202,157]
tau=1000.0
######################################################


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--polystr", type=int, default=-1, help="Exponents in the terms of the fitting polynomial")
args = parser.parse_args()

if args.polystr==-1:
  print("Invalid Input.\nRun `python3 plot_los_signal_3.py --help` for more info")
  sys.exit()

def load_binary_data(filename, dtype=np.float64):
  f = open(filename, "rb") 
  data = f.read() 
  f.close() 
  _data = np.fromstring(data, dtype) 
  if sys.byteorder == 'big':
    _data = _data.byteswap()
  return _data

##############################################################
####Plotting starts here
##############################################################

freq=np.loadtxt("../ska_data/frequency.dat", unpack=True)
freq=[i/1e6 for i in freq]

titles=["eor+fg", "eor+noise", "eor+noise+fg"]

i=los_array[0]
j=los_array[1]

fig = mpl.pyplot.gcf()
fig.set_size_inches(5, 4.0)


####################Plotting the recovered and total field separately


count=311 #Subplot count
for header in ["eor_p_fg", "eor_p_noise", "eor_p_noise_p_fg"]:
  
  axis=plt.subplot(count)
  for offset in [0.1, 1.0, 10.0, 100.0]: #Add here numbers upto 100.0 in factors of 10 starting from 0.1
    
    filename="../results/ska_cube/cube/%s_recovered_21cm_signal_XY_300_%.1f_hrs_Npol_%d_offset_%.6f.bin" % (header, tau, args.polystr, offset)

    data=load_binary_data(filename)

    #ave=np.average(data)
    data=[val-offset for val in data]

    y=[]
    for k in range(21):
      y.append(data[k+21*(j+300*i)])

    #Setting the average to zero
    ave=np.average(y)
    y=[val-ave for val in y]
    
    plt.plot(freq,y,linewidth=1,label=r"offset=%.1f" % offset)
    
  #At this point, offset=100.0
  part=(np.amax(y)-np.amin(y))/5
  yticks=np.arange(np.amin(y)+part,np.amax(y)-part,part)
  #plt.yticks(yticks)
  plt.ylabel(r"$T_{21}$ [K]")
  mpl.rcParams["legend.framealpha"]=0.2
  plt.legend(loc="lower right")
  count+=1

plt.xlabel("Frequency [MHz]")
plt.tight_layout()
plt.subplots_adjust(hspace=0.0)

filepng="ska/los_3/recovered_los_i_%d_j_%d_all_data_all_offset_tau_%.1f_Npol_%d.png" % (i, j, tau,args.polystr)

plt.savefig(filepng)
plt.clf()

    
