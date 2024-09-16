import sys

# Note, need to change path to where git repo is cloned
sys.path.append('/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/')
import rvprec
import pandas as pd
import numpy as np
print('Reading rvprec q_files')

# Note, need to change to where q_info is downloaded
qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')
print('Finished')

def get_hpf_rv_precision(teff,jmag,vsini=3.,exptime=1800.):#,qfiles=None,ffiles=None):
    """
    Note:
        TEFF range is 2500K - 5900K
        
    EXAMPLE:
        qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
        ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')
        get_hpf_rv_precision(2500.,11.3,qfiles=qfiles,ffiles=ffiles)
    """
    #if qfiles is None and ffiles is None:
    #    qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
    #    ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')

    # find closest teff
    
    teffs = pd.DataFrame(qfiles).T['teff'].values
    ind = rvprec.find_nearest_idx(teffs,teff)
    teff_using = teffs[ind]
    print('Input Teff={}K, using closest Teff={}K'.format(teff,teff_using))
    
    #Esorted(pd.DataFrame(qfiles).T['teff'].values)
    qdat, fdat = rvprec.get_data(teff_using,5.0,qfiles,ffiles)
    out = rvprec.calc_prec(qdat,fdat,
                           wl1=8100,
                           wl2=12800,
                           resol=55000,
                           vsini_kms=vsini,
                           rad_m=4.,
                           eff=.004,
                           mag=jmag,
                           magtype='johnson,j',
                           exptime_s=exptime,
                           tell=.95,
                           mask_factor=0.5,
                           sampling_pix_per_resel=3.,
                           beam_height_pix=9.,
                           rdn_pix=6.)
    rv_error = out['dv_all']
    print('Overall precision: {:0.3f} m/s'.format(rv_error))
    return rv_error

def get_espresso_rv_precision(teff,mag,vsini=3.,exptime=1800.,verbose=True,magtype='johnson,v'):#,qfiles=None,ffiles=None):
    """
    Note:
        TEFF range is 2500K - 5900K
        
    EXAMPLE:
        qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
        ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')
        get_neid_rv_precision(2500.,11.3,qfiles=qfiles,ffiles=ffiles)
    """
    #if qfiles is None and ffiles is None:
    #    qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
    #    ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')

    # find closest teff
    
    teffs = pd.DataFrame(qfiles).T['teff'].values
    ind = rvprec.find_nearest_idx(teffs,teff)
    teff_using = teffs[ind]
    if verbose: print('Input Teff={}K, using closest Teff={}K'.format(teff,teff_using))
    
    #Esorted(pd.DataFrame(qfiles).T['teff'].values)
    qdat, fdat = rvprec.get_data(teff_using,5.0,qfiles,ffiles)
    out = rvprec.calc_prec(qdat,fdat,
                           wl1=3800,
                           wl2=7880,
                           resol=140000,
                           vsini_kms=vsini,
                           rad_m=8.2/2.,
                           eff=.022,
                           mag=mag,
                           magtype=magtype,
                           exptime_s=exptime,
                           tell=.99,
                           mask_factor=0.1,
                           sampling_pix_per_resel=5.,
                           beam_height_pix=5.,
                           rdn_pix=4.)
    rv_error = out['dv_all']
    if verbose: print('Overall precision: {:0.3f} m/s'.format(rv_error))
    return rv_error

def get_neid_rv_precision(teff,mag,vsini=3.,exptime=1800.,verbose=True,magtype='johnson,v'):#,qfiles=None,ffiles=None):
    """
    Note:
        TEFF range is 2500K - 5900K
        
    EXAMPLE:
        qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
        ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')
        get_neid_rv_precision(2500.,11.3,qfiles=qfiles,ffiles=ffiles)
    """
    #if qfiles is None and ffiles is None:
    #    qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
    #    ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')

    # find closest teff
    
    teffs = pd.DataFrame(qfiles).T['teff'].values
    ind = rvprec.find_nearest_idx(teffs,teff)
    teff_using = teffs[ind]
    if verbose: print('Input Teff={}K, using closest Teff={}K'.format(teff,teff_using))
    
    #Esorted(pd.DataFrame(qfiles).T['teff'].values)
    qdat, fdat = rvprec.get_data(teff_using,5.0,qfiles,ffiles)
    out = rvprec.calc_prec(qdat,fdat,
                           wl1=4000,
                           wl2=9000,
                           resol=90000,
                           vsini_kms=vsini,
                           rad_m=3.5/2.,
                           eff=.012,
                           mag=mag,
                           magtype=magtype,
                           exptime_s=exptime,
                           tell=.99,
                           mask_factor=0.1,
                           sampling_pix_per_resel=5.,
                           beam_height_pix=5.,
                           rdn_pix=4.)
    rv_error = out['dv_all']
    if verbose: print('Overall precision: {:0.3f} m/s'.format(rv_error))
    return rv_error

def get_maroonx_rv_precision(teff,mag,vsini=3.,exptime=1800.,verbose=True,magtype='johnson,v'):#,qfiles=None,ffiles=None):
    """
    Note:
        TEFF range is 2500K - 5900K
        Should double check with the one in MAROON-X HELP, this one was made fairly quickly
        
    EXAMPLE:
        qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
        ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')
        get_neid_rv_precision(2500.,11.3,qfiles=qfiles,ffiles=ffiles)
    """
    #if qfiles is None and ffiles is None:
    #    qfiles = rvprec.index_qs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output/')
    #    ffiles = rvprec.index_fs(directory='/Users/gks/Dropbox/mypylib/notebooks/GIT/rvprec/data/q_info/output_fluxes/')

    # find closest teff
    
    teffs = pd.DataFrame(qfiles).T['teff'].values
    ind = rvprec.find_nearest_idx(teffs,teff)
    teff_using = teffs[ind]
    if verbose: print('Input Teff={}K, using closest Teff={}K'.format(teff,teff_using))
    
    #Esorted(pd.DataFrame(qfiles).T['teff'].values)
    qdat, fdat = rvprec.get_data(teff_using,4.5,qfiles,ffiles)
    out =  rvprec.calc_prec(qdat,fdat,
                       wl1=5000,
                       wl2=9200,
                       resol=85000,
                       vsini_kms=vsini,
                       rad_m=4.0,
                       eff=.06,
                       mag=mag,
                       magtype='johnson,v',
                       exptime_s=exptime,
                       tell=.95,
                       mask_factor=0.5,
                       sampling_pix_per_resel=3.,
                       beam_height_pix=10.,
                       rdn_pix=16.)
    rv_error = out['dv_all']
    if verbose: 
        print('Overall precision: {:0.3f} m/s'.format(rv_error))
    return rv_error
