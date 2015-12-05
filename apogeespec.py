'''
---------------------------------------------------------------------
APOGEESPEC

Basic tools for dealing with APOGEE spectra of emission line stars.

Written by Drew Chojnowski
---------------------------------------------------------------------
'''

import math
import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from scipy.interpolate import interp1d

specdir='../spectra/'
linefile='linelist.dat'
starcatalog='catalog.dat'

'''
APLOAD: LOAD THE APVISIT WAVE,FLUX INTO A SINGLE ARRAY
'''
def apload(file):
    # read the FITS file
    hdulist=fits.open(file)
    wave=hdulist[4].data
    flux=hdulist[1].data

    # get single arrays of wavelength and flux values, ordered properly
    rw=np.array(wave[0]); rw=rw[::-1]; rw=rw.tolist(); rf=flux[0]; rf=rf[::-1]; rf=rf.tolist()
    gw=np.array(wave[1]); gw=gw[::-1]; gw=gw.tolist(); gf=flux[1]; gf=gf[::-1]; gf=gf.tolist()
    bw=np.array(wave[2]); bw=bw[::-1]; bw=bw.tolist(); bf=flux[2]; bf=bf[::-1]; bf=bf.tolist()
    allwave=np.array(bw+gw+rw); allflux=np.array(bf+gf+rf)

    return allwave,allflux

'''
APCLOAD_CHIPS: LOAD THE APVISIT WAVE,FLUX INTO INDIVIDUAL ARRAYS
'''
def apload_chips(file):
    # read the FITS file
    hdulist=fits.open(file)
    wave=hdulist[4].data
    flux=hdulist[1].data

    # get single arrays of wavelength and flux values, ordered properly
    rw=np.array(wave[0]); rw=rw[::-1]; rw=rw.tolist(); rf=flux[0]; rf=rf[::-1]; rf=rf.tolist()
    gw=np.array(wave[1]); gw=gw[::-1]; gw=gw.tolist(); gf=flux[1]; gf=gf[::-1]; gf=gf.tolist()
    bw=np.array(wave[2]); bw=bw[::-1]; bw=bw.tolist(); bf=flux[2]; bf=bf[::-1]; bf=bf.tolist()
    allwave=np.array(bw+gw+rw); allflux=np.array(bf+gf+rf)

    return rw,rf

'''
BROWSER: INTERACTIVE RA/DEC PLOT
'''
def browser():
    startable=ascii.read(starcatalog)
    abeID=np.array(startable['ID'])
    tmassID=np.array(startable['2MASS'])

    spectra=np.array(glob.glob(specdir+'*apV*fits'))
    nspec=len(spectra)
    stars_all=[]
    ra_all=np.zeros(nspec); dec_all=np.zeros(nspec)
    glon_all=np.zeros(nspec); glat_all=np.zeros(nspec)
    snr_all=np.zeros(nspec); hmag_all=np.zeros(nspec)
    for i in range(nspec):
        hdulist=fits.open(spectra[i])
        stars_all.append(hdulist[0].header['objid'])
        ra_all[i]=hdulist[0].header['ra']
        dec_all[i]=hdulist[0].header['dec']
        glon_all[i]=hdulist[0].header['glon']
        glat_all[i]=hdulist[0].header['glat']
        snr_all[i]=hdulist[0].header['snr']
        hmag_all[i]=hdulist[0].header['h']

    stars_all=np.array(stars_all)
    stars,indices=np.unique(stars_all,return_index=True)
    nstars=len(stars)
    ra=ra_all[indices]
    dec=dec_all[indices]
    glon=glon_all[indices]
    glat=glat_all[indices]
    hmag=hmag_all[indices]

    fig, ax = plt.subplots(figsize=(12,7))
    ax.set_title('click on point to show highest S/N spectrum')
    line,=ax.plot(ra,dec,'o',picker=5,markersize=12)  # 5 points tolerance
    for i in range(nstars):
        p=np.where(stars[i]==tmassID)
        ax.text(ra[i]+3,dec[i],abeID[p][0],fontsize=10)

    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.tight_layout() 

    def onpick(event):
        if event.artist != line:
            return True

        N=len(event.ind)
        if not N:
            return True

        for subplotnum, dataind in enumerate(event.ind):
            star=stars[dataind]
            p=np.where(star==stars_all)
            gd=np.where(snr_all[p]==np.max(snr_all[p]))
            specplot=apsplot(spectra[p][gd][0],mark_sky=True)
        return True

    fig.canvas.mpl_connect('pick_event', onpick)

    return

'''
APSPLOT_ON_KEY: KEY PRESS EVENT HANDLER FOR APSPLOT
'''
def apsplot_on_key(event):
    linelist=ascii.read(linefile)
    global key, xdata, ydata
    key=event.key
    xdata=event.xdata; ydata=event.ydata
    xdata_str="%0.3f" % xdata
    ydata_str="%0.3f" % ydata

     # find the nearest line in the linelist
    dist=abs(xdata-linelist['CENT'])
    x=np.where(dist==min(dist))
    z=x[0].astype('int')
    rest=linelist['CENT'][z]
    rest_str="%0.3f" % rest[0]
    dif=rest[0]-xdata
    dif_str="%0.3f" % dif
    closeline=linelist['LABEL'][z][0].tolist()
    if abs(dif) > 8:
        print(key)
        print(xdata_str+r'....No lines within 10 $\AA$')
    else:
        print(key)
        bla=xdata_str+'....nearest line= '+closeline+' '+rest_str+', dif= '+dif_str+')'
        print(bla)

'''
APSPLOT: SPECTRUM PLOTTING PROGRAM
'''
def apsplot(file,mark_sky=False,mark_lines=True,xshift=0.0,do_cont=True,do_cont_chips=False,cont_ord=5):
    linelist=ascii.read(linefile)    # read the linelist

    wave,flux=apload(file)    # get the wavelength and flux arrays
    wave=wave+xshift    # option to add on an xshift, in angstroms

    # get some header values from the FITS file
    hdulist=fits.open(file)
    objid=hdulist[0].header['objid']
    snr=str(hdulist[0].header['snr'])
    mjd=str(hdulist[0].header['mjd5'])
    exptime=str(hdulist[0].header['exptime'])

    # option to normalize the spectrum
    if do_cont is True:
        # option to normalize the chips sepately
        if do_cont_chips is True:
            bw=np.array(bw); bf=np.array(bf)
            coeff=np.polyfit(bw,bf,cont_ord)
            poly=np.poly1d(coeff)
            cont=poly(bw)
            bf=bf/cont
            bf=bf.tolist()

            gw=np.array(gw); gf=np.array(gf)
            coeff=np.polyfit(gw,gf,cont_ord)
            poly=np.poly1d(coeff)
            cont=poly(gw)
            gf=gf/cont
            gf=gf.tolist()

            rw=np.array(rw); rf=np.array(rf)
            coeff=np.polyfit(rw,rf,cont_ord)
            poly=np.poly1d(coeff)
            cont=poly(rw)
            rf=rf/cont
            rf=rf.tolist()

            wave=np.array(bf+gf+rf)
        # else normalize the whole thing at once
        else:
            coeff=np.polyfit(wave,flux,cont_ord)
            poly=np.poly1d(coeff)
            cont=poly(wave)
            flux=flux/cont

    # make the plot
    fig=plt.figure(figsize=(16,8))
    matplotlib.rcParams.update({'font.size': 14, 'font.family':'serif'})
    ax=fig.add_subplot(111)
    ax.plot(wave[np.where((wave>15144)&(wave<15808))],flux[np.where((wave>15144)&(wave<15808))],color='black')
    ax.plot(wave[np.where((wave>15858)&(wave<16432))],flux[np.where((wave>15858)&(wave<16432))],color='black')
    ax.plot(wave[np.where((wave>16473)&(wave<16955))],flux[np.where((wave>16473)&(wave<16955))],color='black')
#    ax.plot(wave,flux,color='black')
    ax.grid(True)
    plt.xlim([15135,16965])
    plt.ylim([0.5,1.5])
    plt.xlabel(r'Observed Wavelength [$\AA$]')
    plt.ylabel(r'Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA$]')
    plt.title(file+'    '+objid+'    S/N='+snr+'    exptime='+exptime+' s')
    plt.tight_layout()

    # option to mark airglow lines
    if mark_sky is True:
        bd=np.where(linelist['LABEL']=='airglow')
        bdlines=linelist[bd]
        for i in range(len(bdlines)):
            pos=bdlines['CENT'][i]+xshift
            sec=np.where(wave>(pos-2))
            wtmp=wave[sec]; ftmp=flux[sec]
            sec=np.where(wtmp<(pos+2))
            wtmp=wtmp[sec]; ftmp=ftmp[sec]
            plt.plot(wtmp,ftmp,color='red')

    # option to mark stellar lines
    if mark_lines is True:
        gd=np.where(linelist['LABEL']!='airglow')
        gdlines=linelist[gd]
        for i in range(len(gdlines)):
            line=gdlines['CENT'][i]
            lab=gdlines['LABEL'][i]
            if lab[0:1]=='H': labcol='blue'
            if lab[0:1]!='H': labcol='green'
            sec=np.where(abs(wave-line)<0.5)
            arrowstart=np.mean(flux[sec])+((max(flux)-min(flux))*0.10)
            arrowlen=(max(flux)-min(flux))*(-0.05)
            arrowheadL=(max(flux)-min(flux))*(0.01)
            sec=np.where(abs(wave-line)<0.5)
            txty=arrowstart+(max(flux)-min(flux))*(0.015)
            if do_cont is True and txty<1: 
                ax.arrow(line,1.09,0,-0.05,head_width=3,head_length=arrowheadL,color=labcol)
                ax.text(line,1.1,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)
            else:
                ax.arrow(line,arrowstart,0,arrowlen,head_width=3,head_length=arrowheadL,color=labcol)
                ax.text(line,txty,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)

    cid=fig.canvas.mpl_connect('key_press_event',apsplot_on_key)

    hdulist.close()
    plt.show()

    return
