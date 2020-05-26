#To plot 2D PS for a given tau, polystr, offset and header
#plots for given total_signal, total_noise, recovered_signal

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


parser = argparse.ArgumentParser()

parser.add_argument("-p", "--polystr", type=int, default=-1, help="Exponents in the terms of the fitting polynomial")
parser.add_argument("-d", "--header", type=str, default=-1, help="Header for the SKA data")
parser.add_argument("-o", "--offset", type=float, default=-1, help="Offset added to the original cube")
parser.add_argument("-s", "--signal", type=str, default=-1, help="Enter recovered_21cm_signal, noise_signal or total_21cm_signal")
parser.add_argument("-t", "--tau", type=float, default=-1, help="Integration time (in hrs)")
args = parser.parse_args()

if any(i==-1 for i in [args.polystr,args.offset,args.header,args.signal,args.tau]):
  print("Invalid Input.\nRun `python3 plot_2dps_2.py --help` for more info")
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

xfile="../results/ska_cube/ps/%s_k_z_XY_300_%.1f_hrs_Npol_%d_offset_%.6f.dat" % (args.header, args.tau, args.polystr, args.offset)
yfile="../results/ska_cube/ps/%s_k_perp_XY_300_%.1f_hrs_Npol_%d_offset_%.6f.dat" % (args.header, args.tau, args.polystr, args.offset)

x=np.loadtxt(xfile, unpack=True)
y=np.loadtxt(yfile, unpack=True)

#fig = mpl.pyplot.gcf()
#fig.set_size_inches(5, 6.0)
#filename="../results/ska_cube/ps/2d_ps_noise_field_XY_300.bin"

filename="../results/ska_cube/ps/2d_ps_%s_%s_XY_300_%.1f_hrs_Npol_%d_offset_%.6f.bin" % (args.header, args.signal, args.tau, args.polystr, args.offset)

data=load_binary_data(filename)
upper=data.max()/1000
lower=data.min()
np.clip(data,lower,upper,data)

Z=data.reshape((len(x),len(y)))

#plt.figure()
#ax = plt.subplot(111)
im=plt.imshow(Z,extent=(min(y),max(y),max(x),min(x)),interpolation='nearest')#,cmap=plt.cm.hot)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)

#plt.colorbar()
plt.colorbar(im,fraction=0.026, pad=0.04)
plt.clim(lower,upper)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$k_{\perp}$ $[Mpc^{-1}h]$")
plt.ylabel(r"$k_{\parallel}$ $[Mpc^{-1}h]$")
plt.tight_layout()

fname="ska/2dps_2/%s_%s_tau_%.1f_Npol_%d_%.6f.png" % (args.header, args.signal, args.tau, args.polystr, args.offset)
#plt.show()
#fname="ska/2dps_2/noise_tau_%.1f.png" % args.tau
plt.savefig(fname)
