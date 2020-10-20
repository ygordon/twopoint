####perform 2-point correlation function based on data set
####output data file for plotting
###based on command line args and config file

import numpy as np, pandas as pd, argparse, os
from astroML.correlation import bootstrap_two_point
from astropy.table import Table
from astropy import units as u


####################################################################################
####################################################################################
##def functions

def midbins(bins):
    'returns middle value from bin edges for plotting'
    x = []
    for i in range(len(bins)-1):
        x0, x1 = bins[i], bins[i+1]
        midx = x0 + ((x1-x0)/2)
        x.append(midx)
    
    return(np.array(x))


def two_pcf(ra, dec, N_boot, sepbins=np.logspace(-3.5, 1, 30)):
    'determine 2-point correlation (LS estimator) and bootstrap errors'
    
    ###make ra dec into single Nx2 array
    ra, dec, = np.array(ra), np.array(dec)
    skypos = np.vstack((ra, dec)).transpose()
    
    ###perform 2 point and bootstrap errors
    ls_vals, ls_err = bootstrap_two_point(skypos, sepbins, Nbootstrap=N_boot,
                                          method='landy-szalay')
                                 
    ###determin middle of bin for x vals corresponding to LS (y) vals
    x_vals = midbins(sepbins)
    
    return(x_vals, ls_vals, ls_err)


#def main(data_file, racol='RA', deccol='DEC', N_boot=1000, minsep=1,
#         maxsep=36000, n_bins=30):
def main(data_file, config_file):
    'load position data from file, run 2point and write results to file'
    
    ###load config parameters
    cp = pd.read_table(config_file, delim_whitespace=True).replace("'", "", regex=True)
    racol = cp['value'][np.where(cp['parameter']=='ra_col')[0][0]]
    deccol = cp['value'][np.where(cp['parameter']=='dec_col')[0][0]]
    minsep = np.float(cp['value'][np.where(cp['parameter']=='min_sep_arcsec')[0][0]])
    maxsep = np.float(cp['value'][np.where(cp['parameter']=='max_sep_arcsec')[0][0]])
    n_bins = np.int(cp['value'][np.where(cp['parameter']=='N_bins')[0][0]])
    N_boot = np.int(cp['value'][np.where(cp['parameter']=='N_bs')[0][0]])
    
    ###define separation bins (input arcsec -> deg)
    bins = np.logspace(np.log10(minsep/3600), np.log10(maxsep/3600), n_bins)
    
    ###load data and run 2 point
    data = Table.read(data_file, format='fits')
    x, y, yerr = two_pcf(ra=data[racol], dec=data[deccol], sepbins=bins, N_boot=N_boot)
    
    ###convert x to arcsec, create table (add units to x)
    x = x*3600
    
    outtable = Table(np.vstack((x, y, yerr)).transpose(),
                     names=('Separation', 'w_LS', 'E_w_LS'))
    
    outtable['Separation'].unit = u.arcsec
    
    ###write to file - warn if same file name already exists
    ### this will only work first time, intended as a backup to save rerunning from scratch
    outfilename = data_file.split('.')[0] + '_2point_data.fits'
    dlist = os.listdir()
    if outfilename in dlist:
        of2 = outfilename.split('.')[0] + '_2.fits'
        print('WARNING: ' + outfilename + 'already exists, data output as '
              + of2 + ' instead')
        outfilename = of2
    
    outtable.write(outfilename, format='fits')
    
    return


####################################################################################
####################################################################################
##main code

###parse command line arguments

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('data', type=str,
                        help='list of positional data (e.g. catalogue file)')
    parser.add_argument('config', type=str,
                        help='config file')
    args = parser.parse_args()

    main(data_file=args.data, config_file=args.config)

