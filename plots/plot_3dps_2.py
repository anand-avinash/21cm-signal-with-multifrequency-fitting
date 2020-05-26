#To plot 3D PS for given polystr, tau and offset
#Plots eor_p_fg, eor_p_noise, eor_p_noise_fg in different panels for both total and recovered signal

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import matplotlib as mpl


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--polystr", type=int, default=-1, help="Exponents in the terms of the fitting polynomial")
parser.add_argument("-o", "--offset", type=float, default=-1, help="Offset added to the original cube")
parser.add_argument("-t", "--tau", type=float, default=-1, help="Integration time (in hrs)")
args = parser.parse_args()

if any(i==-1 for i in [args.polystr,args.offset,args.tau]):
  print("Invalid Input.\nRun `python3 plot_3dps_2.py --help` for more info")
  sys.exit()


##############################################################
####Plotting starts here
##############################################################

titles=["eor+fg", "eor+noise", "eor+noise+fg"]
signal_arr=["Recovered 21cm signal", "Total 21cm signal"]

fig = mpl.pyplot.gcf()
fig.set_size_inches(5, 6.0)


count=311 #Subplot count
for header in ["eor_p_fg", "eor_p_noise", "eor_p_noise_p_fg"]:
  
  sig_count=0
  axis=plt.subplot(count)
  for signal in ["recovered_21cm_field", "total_21cm_field"]:
    
    filename="../results/ska_cube/ps/3d_ps_%s_%s_XY_300_%.1f_hrs_Npol_%d_offset_%.6f.dat" % (header, signal, args.tau, args.polystr, args.offset)

    k, ps = np.loadtxt(filename, usecols=(1,2,) , unpack=True)

    plt.plot(k,ps,linewidth=1,label="%s" % signal_arr[sig_count])
    sig_count+=1
    
  plt.yscale('log')
  plt.xscale('log')
  #At this point, signal=total_21cm_field
  #part=(np.amax(ps)-np.amin(ps))/5
  #yticks=np.arange(np.amin(ps)+part,np.amax(ps)-part,part)
  #plt.yticks(yticks)
  plt.ylabel(r"$\Delta^2(k)$")
  mpl.rcParams["legend.framealpha"]=0.2
  plt.legend(loc="lower right")
  count+=1

plt.xlabel(r"k [$Mpc^{-1}h$]")
plt.tight_layout()
plt.subplots_adjust(hspace=0.0)

filepng="ska/3dps_2/3dps_Npol_%d_tau_%.1f_offset_%.6f.png" % (args.polystr, args.tau, args.offset)

plt.savefig(filepng)
plt.clf()