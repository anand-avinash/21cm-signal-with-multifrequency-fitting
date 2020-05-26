#To plot los signal for a given los, offset and polystr.
#Plots eor_p_fg, eor_p_noise, eor_p_noise_fg in different panels for all tau

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import matplotlib as mpl
import matplotlib.ticker as mticker


######################################################
#### Choose (x,y) for plotting los signal
los_array=[0,145]#[156,164]#[23,234]
######################################################


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--polystr", type=int, default=-1, help="Exponents in the terms of the fitting polynomial")
parser.add_argument("-o", "--offset", type=float, default=-1, help="Offset added to the original cube")
args = parser.parse_args()

if any(i==-1 for i in [args.polystr,args.offset]):
  print("Invalid Input.\nRun `python3 plot_los_signal_2.py --help` for more info")
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
fig.set_size_inches(5, 6.0)


####################Plotting the recovered and total field separately

for fcount in [0,1]:
  count=311 #Subplot count
  #axis1=plt.subplot(count)
  for header in ["eor_p_fg", "eor_p_noise", "eor_p_noise_p_fg"]:
    axis=plt.subplot(count)#,sharex=axis1)

    for tau in [1000.0, 100.0, 10.0]:
    
      if fcount==0:
        filename="../results/ska_cube/cube/%s_recovered_21cm_signal_XY_300_%.1f_hrs_Npol_%d_offset_%.6f.bin" % (header, tau, args.polystr, args.offset)
      else:
        filename="../ska_data/%s_pre_total_field_300_%.1f_hr-image.bin" %(header, tau)

      data=load_binary_data(filename)

      #if fcount==0:
      #  data=[val-args.offset for val in data]

      print(filename)
      print(np.average(data))

      y=[]
      for k in range(21):
        y.append(data[k+21*(j+300*i)])

      plt.plot(freq,y,linewidth=1,label=r"$\tau=%.1f$ hrs" % tau)
    
    #At this point, tau=10.0
    part=(np.amax(y)-np.amin(y))/5
    yticks=np.arange(np.amin(y)+part,np.amax(y)-part,part)
    plt.yticks(yticks)
    #print(yticks)
    plt.ylabel(r"$T_{21}$ [K]")
    mpl.rcParams["legend.framealpha"]=0.2
    plt.legend(loc="lower right")
    count+=1

  plt.xlabel("Frequency [MHz]")
  plt.tight_layout()
  plt.subplots_adjust(hspace=0.0)
  
  if fcount==0:
    filepng="ska/los_2/recovered_los_i_%d_j_%d_all_data_all_tau_Npol_%d_offset_%.6f.png" % (i, j, args.polystr, args.offset)
  else:
    filepng="ska/los_2/total_los_i_%d_j_%d_all_data_all_tau_Npol_%d_offset_%.6f.png" % (i, j, args.polystr, args.offset)

  plt.savefig(filepng)
  plt.clf()
    
