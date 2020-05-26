#To plot 2D PS for given eor, eor + free-free, eor + free-free + point sources,  eor + total foreground, and eor + foreground + detector noise fields

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--polystr", type=int, default=-1, help="Exponents in the terms of the fitting polynomial")
parser.add_argument("-f", "--frequency", type=int, default=-1, help="Frequency related to the simulation")
parser.add_argument("-s", "--signal", type=str, default=-1, help="Enter original_21cm_field, original_p_ff, original_p_ff_p_point, original_p_fg, total_21cm_field")
parser.add_argument("-t", "--tau", type=float, default=-1, help="Integration time (in hrs)")
args = parser.parse_args()

if any(i==-1 for i in [args.polystr,args.frequency,args.signal,args.tau]):
  print("Invalid Input.\nRun `python3 plot_2dps.py --help` for more info")
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

xfile="../results/additional/ps/k_z_nu_%.0f_XY_128_Z_100_%.1f_hrs_Npol_%d.dat" % (args.frequency, args.tau, args.polystr)
yfile="../results/additional/ps/k_perp_nu_%.0f_XY_128_Z_100_%.1f_hrs_Npol_%d.dat" % (args.frequency, args.tau, args.polystr)

x=np.loadtxt(xfile, unpack=True)
y=np.loadtxt(yfile, unpack=True)


filename="../results/additional/ps/2d_ps_%s_nu_%.0f_XY_128_Z_100_%.1f_hrs_Npol_%d.bin" % (args.signal, args.frequency, args.tau, args.polystr)

data=load_binary_data(filename)
upper=data.max()/10
lower=data.min()
np.clip(data,lower,upper,data)
#1/100 for noise, and total_field

Z=data.reshape((len(x),len(y)))

#plt.figure()
#ax = plt.subplot(111)
im=plt.imshow(Z,extent=(min(y),max(y),max(x),min(x)),interpolation='nearest',vmin=lower, vmax=upper)#,cmap=plt.cm.hot)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)

#plt.colorbar()
plt.colorbar(im,fraction=0.026, pad=0.04, extend="max")
#plt.clim(lower,upper)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$k_{\perp}$ $[Mpc^{-1}h]$")
plt.ylabel(r"$k_{\parallel}$ $[Mpc^{-1}h]$")
plt.tight_layout()

fname="add_ps/2d_%.0f/%s_%.0f_tau_%.1f_Npol_%d.png" % (args.frequency, args.signal, args.frequency, args.tau, args.polystr)
#plt.show()
#fname="ska/2dps_2/noise_tau_%.1f.png" % args.tau
plt.savefig(fname)
