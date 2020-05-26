# Statistical Detection of 21 cm Signal from the Epoch of Reionization

Here is the program I wrote for my master's thesis project. This program uses the multi-frequency fitting method to extract the 21 cm signal blindly from a contaminated signal. It also covers an attempt to recover the 21cm signal power spectrum from the contaminated signal cube.

## Program description

- Usage:
  1. Modify the compiler flags in `compile.sh` as per your system specs
  2. Run the shell script `compile.sh`:

     ```shell
     sh compile.sh
     ```

  3. Run the executable `execute` with appropriate arguments:

     ```shell
     ./execute [frequency (in MHz)] [XY_pix] [Z_pix] [Int_time (hrs)] [expo_array] #To simulate the 21 cm fields
     ```

     OR

     ```shell
     ./execute [header] [XY_pix] [Int_time (hrs)] -o [offset] [expo_array] #To read the 21 cm fields from SKA blind challenge data
     ```

- You must have FFTW and GSL installed. See  
  FFTW: <http://www.fftw.org/download.html>  
  GSL: <https://www.gnu.org/software/gsl/>
- Simulation parameters to simulate 21 cm signal and foreground contaminations are stored in `params.h` and `foreg_params.h` respectively. Cosmological parameters to compute the comoving distance between two redshifts are stored in `power_params.h`.

## Notes

- You'll need dimensionless power spectrum in `ps_table/` to simulate the 21 cm field.
- The SKA blind challenge dataset was scaled down and had positive and negative data points as well. The multi-frequency fitting method requires us to compute the logarithm of the data points and then attempts to fit the foreground out of it. So in order to apply the foreground removal technique, a constant positive offset was added to all the data points prior to the foreground removal.

## Plotting description

- `plots/plot_2dps_1.py`  
  Plots the 2D power spectrum (Datta, Bowman) of simulated eor, eor + free-free, eor + free-free + point sources, and eor + total foreground fields
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

1. Li-Ping He. Foreground removal of 21 cm fluctuation with multifrequency fitting. *Research in Astronomy and Astrophysics 9.6 (2009), p. 653*.
2. Max Tegmark et al. Foregrounds and Forecasts for the Cosmic Microwave Background. *The Astrophysical Journal 530.1 (2000), p. 133*.
3. A. Datta, J. D. Bowman, and C. L. Carilli. Bright Source Subtraction Requirements for Redshifted 21 cm Measurements. *The Astrophysical Journal 724 (Nov. 2010), pp. 526–538*.
4. Adrian Liu, Aaron R. Parsons, and Cathryn M. Trott. Epoch of reionization window. II. Statistical methods for foreground wedge reduction. *Phys. Rev. D 90 (2 2014), p. 023019*.
