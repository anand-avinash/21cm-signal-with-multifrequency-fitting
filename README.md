# Statistical Detection of 21 cm Signal from the Epoch of Reionization

Here is the program I wrote for my master's thesis project. This program uses the multi-frequency fitting method to extract the 21 cm signal blindly from a contaminated signal [[1](#ref1),[2](#ref2)]. It also covers an attempt to recover the 21cm signal power spectrum from the contaminated signal cube [[3](#ref3),[4](#ref4)].

## Program description

- Usage:
  1. Modify the compiler flags in `compile.sh` as per your system specs
  2. Run the shell script `compile.sh`:

     ```shell
     sh compile.sh
     ```

  3. Run the executable `execute` with appropriate arguments:

     ```shell
     #To simulate the 21 cm fields
     ./execute [frequency (in MHz)] [XY_pix] [Z_pix] [Int_time (hrs)] [expo_array]
     ```

     OR

     ```shell
     #To read the 21 cm fields from SKA blind challenge data
     ./execute [header] [XY_pix] [Int_time (hrs)] -o [offset] [expo_array]
     ```

- You must have FFTW and GSL installed. See  
  FFTW: <http://www.fftw.org/download.html>  
  GSL: <https://www.gnu.org/software/gsl/>
- Simulation parameters to simulate 21 cm signal and foreground contaminations are stored in `params.h` and `foreg_params.h` respectively. Cosmological parameters to compute the comoving distance between two redshifts are stored in `power_params.h`.

## Notes

- You'll need dimensionless power spectrum in `ps_table/` to simulate the 21 cm field.
- `ska_data/` stores the SKA blind challenge data. See `ska_data/.gitignore` for the list of files that needs to be there.
- `.bin` files need to be generated by combining the datasets from the original `.fits` files of the SKA blind challenge data. To do so, copy the `SKA_Blind_Challenge_Pilot/images` folder to `ska_data/` and run the `ska_data/codes/fit_bin.py`. This will copy the signal data from the `.fits` files to `.bin` in `ska_data/` with the help of `ska_data/codes/box_copy.c`.
- The SKA blind challenge dataset was scaled down and had positive and negative data points as well. The multi-frequency fitting method requires computing the logarithm of the data points to fit the foreground. So in order to apply the foreground removal technique, a constant positive offset was added to all the data points prior to the foreground removal.

## Plotting description

- `plots/plot_2dps_1.py`  
  Plots the 2D power spectrum [[3](#ref3),[4](#ref4)] of simulated eor, eor + free-free, eor + free-free + point sources, and eor + total foreground fields
- `plots/plot_2dps_2.py`  
  Plots the 2D power spectrum of SKA blind challenge data for a given tau, polystr, offset and header
- `plots/plot_3dps_1.py`  
  Plots the 3D power spectrum of simulated eor, eor + free-free, eor + free-free + point sources, and eor + total foreground fields
- `plots/plot_3dps_2.py`  
  Plots the 3D power spectrum of SKA blind challenge data for given polystr, tau and offset. Plots eor_p_fg, eor_p_noise, eor_p_noise_fg in different panels for total and recovered signal together
- `plots/plot_los_signal_1.py`  
  Plots line-of-sight (los) signal for 4 given los for a given offset, polystr and tau
- `plots/plot_los_signal_2.py`  
  Plots los signal for a given los, offset and polystr. Plots eor_p_fg, eor_p_noise, eor_p_noise_fg in different panels for all tau
- `plots/plot_los_signal_3.py`  
  Plots los signal for a given los, tau and polystr. Plots eor_p_fg, eor_p_noise, eor_p_noise_fg in different panels for all offsets

## References

1. <a id="ref1"></a> Li-Ping He. Foreground removal of 21 cm fluctuation with multifrequency fitting. *Research in Astronomy and Astrophysics 9.6 (2009), p. 653*.
2. <a id="ref2"></a> Max Tegmark et al. Foregrounds and Forecasts for the Cosmic Microwave Background. *The Astrophysical Journal 530.1 (2000), p. 133*.
3. <a id="ref3"></a> A. Datta, J. D. Bowman, and C. L. Carilli. Bright Source Subtraction Requirements for Redshifted 21 cm Measurements. *The Astrophysical Journal 724 (Nov. 2010), pp. 526–538*.
4. <a id="ref4"></a> Adrian Liu, Aaron R. Parsons, and Cathryn M. Trott. Epoch of reionization window. II. Statistical methods for foreground wedge reduction. *Phys. Rev. D 90 (2 2014), p. 023019*.
