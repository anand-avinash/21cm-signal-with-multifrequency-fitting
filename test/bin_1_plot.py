#An example to demonstrate the plotting of 2D pixel plot from array saved as bin file in row major order
#Refer to bin_1.c

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys


file_test="bindata_1.bin"

def load_binary_data(filename, dtype=np.float64): 
  """ 
  We assume that the data was written 
  with write_binary_data() (little endian). 
  """ 
  f = open(filename, "rb") 
  data = f.read() 
  f.close() 
  _data = np.fromstring(data, dtype) 
  if sys.byteorder == 'big':
    _data = _data.byteswap()
  return _data

data=load_binary_data(file_test)

x=range(1,6,1)
y=range(1,3,1)
plt.xscale("log")
plt.yscale("log")

Z=data.reshape((2,5))
print(Z)
plt.imshow(Z,extent=(min(x),max(x),max(y),min(y)),interpolation='nearest',cmap=plt.cm.hot)
plt.colorbar()
plt.show()
