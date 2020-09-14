**two_point.py readme**

command line code for running two point correlation function (Landy-Szalay estimator) with bootstrapped errors
takes data file (fits format) with postitions in and config file inputs as command line arguments, run as:

*python3 two_point.py data.fits config.txt*

returns file with the data for 2 point correlation (angular separation, w_LS, and error in w_LS)


**Config File Description**

parameter | description
----------|------------
ra_col | name of the RA column in the data file (assumed to be in decimal degrees)
dec_col | name of the Dec column in the data file (assumed to be in decimal degrees)
min_sep_arcsec | minimum of range of angular separations to perform two-point correlation for (arcsec)
max_sep_arcsec | maximum of range of angular separations to perform two-point correlation for (arcsec)
N_bins | Number of angular separation bins to use for two-point correlation (bins are in log-space)
N_bs | Number or resamples to use for bootstrap error estimation
