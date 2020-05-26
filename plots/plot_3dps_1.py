#To plot 3D PS for eor, eor + free-free, eor + free-free + point sources,  eor + total foreground, and eor + foreground + detector noise fields
#Plots 50, 100, 150, 200 in different panels for data in additional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


tau=1.0

signal_arr=["eor", "eor + free-free", "eor + free-free + point source", "eor + fg", "Total signal"]#, "Recovered signal"]

fig = mpl.pyplot.gcf()
fig.set_size_inches(6, 7.0)


count=411 #Subplot count
for frequency in [50, 100, 150, 200]:
  
  sig_count=0
  axis=plt.subplot(count)
  for signal in ["original_21cm_field", "original_p_ff", "original_p_ff_p_point", "original_p_fg", "total_21cm_field"]:#, "recovered_21cm_field"]:
    
    filename="../results/additional/ps/3d_ps_%s_nu_%d_XY_128_Z_100_%.1f_hrs_Npol_5310.dat" % (signal, frequency,tau)

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
  #plt.legend(loc="upper center",prop={'size': 8},bbox_to_anchor=(0.5, -0.05), ncol=5)
  count+=1

plt.legend(loc='center left', bbox_to_anchor=(1, 2.))#,prop={'size': 8})
plt.xlabel(r"k [$Mpc^{-1}h$]")
#plt.tight_layout()
plt.subplots_adjust(hspace=0.0)

filepng="add_ps/3d/3dps_Npol_5310_tau_%.1f.png" % tau

#plt.show()
plt.savefig(filepng,bbox_inches="tight")
plt.clf()