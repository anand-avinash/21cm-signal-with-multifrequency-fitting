#Reads the FITS one by one and copy the data to the dat file for each. Also creates a list of frequencies.

import re
import subprocess
from astropy.io import fits

N=300

for header in ["eor_p_fg","eor_p_noise","eor_p_noise_p_fg","noise"]:

  folder="../images/%s" % header

  for int_time in [10.0,100.0,1000.0]:

    freq=[]
    ls1=subprocess.Popen(["ls","-p",folder],stdout=subprocess.PIPE)
    grep1=subprocess.Popen(["grep","-v","/"],stdin=ls1.stdout,stdout=subprocess.PIPE)
    grepout1=grep1.stdout

    for files in grepout1:

      name="%s_(.+)_Hz_t_obsv_%1.1f_hr-image.fits" % (header,int_time)
      fileobject=re.match(name,files.decode("utf-8"))
      if fileobject:
        frequency=int(fileobject.group(1))
        freq.append(frequency)

        dest="%s/%s" % (folder, fileobject.group())
        hdul=fits.open(dest)
        hdul.info()
        data=hdul[0].data
        a=[]
        for i in range(N):
          for j in range(N):
            a.append(data[0,0,i][j])
        hdul.close()

        save_dest="../%s.txt" % fileobject.group()
        f=open(save_dest,"w")
        for i in range(N*N):
          f.write("%0.8lE\n" % a[i])
        f.close()

freq_dest="../frequency.dat"
f=open(freq_dest,"w")
for i in freq:
  f.write("%d\n" % i)
f.close()

subprocess.call(["gcc", "box_copy.c", "-o", "execute"])

for header in ["eor_p_fg","eor_p_noise","eor_p_noise_p_fg","noise"]:
  for int_time in [10.0,100.0,1000.0]:
    subprocess.call(["./execute", header, str(300), str(int_time)])

subprocess.call(["rm ../*.txt execute"], shell=True)