#To plot los signal for 4 given los for a given offset, polystr and tau

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

######################################################
#### Choose (x,y) for plotting los signal
los_array=[[34,54], [12,183], [234,145], [218,87]]
######################################################


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--polystr", type=int, default=-1, help="Exponents in the terms of the fitting polynomial")
parser.add_argument("-o", "--offset", type=float, default=-1, help="Offset added to the original cube")
parser.add_argument("-t", "--tau", type=float, default=-1, help="Integration time (in hrs)")
args = parser.parse_args()

if any(i==-1 for i in [args.polystr,args.offset,args.tau]):
  print("Invalid Input.\nRun `python3 plot_los_signal.py --help` for more info")
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

#Plan:
#recovered/total[0]=eor_p_fg signal
#recovered/total[1]=eor_p_noise signal
#recovered/total[2]=eor_p_noise_p_fg signal

recovered=[None,None,None]
total=[None,None,None]

k=0
for i in ["eor_p_fg", "eor_p_noise", "eor_p_noise_p_fg"]:
  file1="../results/ska_cube/cube/%s_recovered_21cm_signal_XY_300_%.1f_hrs_Npol_%d_offset_%.6f.bin" % (i, args.tau, args.polystr, args.offset)
  file2="../ska_data/%s_pre_total_field_300_%.1f_hr-image.bin" % (i, args.tau)
  recovered[k]=load_binary_data(file1)
  total[k]=load_binary_data(file2)
  k+=1
fig=plt.figure()
ax=fig.add_subplot(1,1,1)

for los in los_array:

  y1=[]; y2=[]; y3=[]; y4=[]; y5=[]; y6=[]
  i=los[0]
  j=los[1]

  for k in range(21):
    y1.append(recovered[0][k+21*(j+300*i)])
    y2.append(recovered[1][k+21*(j+300*i)])
    y3.append(recovered[2][k+21*(j+300*i)])
    y4.append(total[0][k+21*(j+300*i)])
    y5.append(total[1][k+21*(j+300*i)])
    y6.append(total[2][k+21*(j+300*i)])
  
  ax.plot(freq,y1,"r--",label="eor+fg")
  ax.plot(freq,y2,"b:",label="eor+noise")
  ax.plot(freq,y3,"g-.",label="eor+noise+fg")
  plt.xlabel("Frequency [MHz]")
  plt.ylabel(r"$T_{21}$ [K]")
  plt.tight_layout()
  plt.legend(loc="best")
  fig.savefig("ska/los_1/los_x_%d_y_%d_recovered_all_tau%.1f_polystr%d_offset%.6f.png"%(i,j,args.tau,args.polystr,args.offset))
  plt.cla()

  ax.plot(freq,y4,"r--",label="eor+fg")
  ax.plot(freq,y5,"b:",label="eor+noise")
  ax.plot(freq,y6,"g-.",label="eor+noise+fg")
  plt.xlabel("Frequency [MHz]")
  plt.ylabel(r"$T_{21}$ [K]")
  plt.tight_layout()
  plt.legend(loc="best")
  fig.savefig("ska/los_1/los_x_%d_y_%d_total_all_tau%.1f_polystr%d_offset%.6f.png"%(i,j,args.tau,args.polystr,args.offset))
  plt.cla()

  #Setting average of all the fields to zero
  ave=sum(y1)/len(y1)
  y1=[i-ave for i in y1]

  ave=sum(y2)/len(y2)
  y2=[i-ave for i in y2]

  ave=sum(y3)/len(y3)
  y3=[i-ave for i in y3]

  ave=sum(y4)/len(y4)
  y4=[i-ave for i in y4]

  ave=sum(y5)/len(y5)
  y5=[i-ave for i in y5]

  ave=sum(y6)/len(y6)
  y6=[i-ave for i in y6]

  ax.plot(freq,y1,"r--",label="eor+fg, recovered")
  ax.plot(freq,y4,"b:",label="eor+fg, total")
  plt.xlabel("Frequency [MHz]")
  plt.ylabel(r"$T_{21}$ [K]")
  plt.tight_layout()
  plt.legend(loc="best")
  fig.savefig("ska/los_1/los_x_%d_y_%d_eor_p_fg_tau%.1f_polystr%d_offset%.6f.png"%(i,j,args.tau,args.polystr,args.offset))
  plt.cla()

  ax.plot(freq,y2,"r--",label="eor+noise, recovered")
  ax.plot(freq,y5,"b:",label="eor+noise, total")
  plt.xlabel("Frequency [MHz]")
  plt.ylabel(r"$T_{21}$ [K]")
  plt.tight_layout()
  plt.legend(loc="best")
  fig.savefig("ska/los_1/los_x_%d_y_%d_eor_p_noise_tau%.1f_polystr%d_offset%.6f.png"%(i,j,args.tau,args.polystr,args.offset))
  plt.cla()

  ax.plot(freq,y3,"r--",label="eor+noise+fg, recovered")
  ax.plot(freq,y6,"b:",label="eor+noise+fg, total")
  plt.xlabel("Frequency [MHz]")
  plt.ylabel(r"$T_{21}$ [K]")
  plt.tight_layout()
  plt.legend(loc="best")
  fig.savefig("ska/los_1/los_x_%d_y_%d_eor_p_noise_fg_tau%.1f_polystr%d_offset%.6f.png"%(i,j,args.tau,args.polystr,args.offset))
  plt.cla()