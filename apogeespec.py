import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from scipy.interpolate import interp1d

def plotspec(file,linefile='linelist.dat',mark_sky=False,mark_lines=True,xshift=0.0):
    # read the FITS file
    hdulist=fits.open(file)
    wave=hdulist[4].data+xshift
    flux=hdulist[1].data
    bwave=wave[2].tolist()
    gwave=wave[1].tolist()
    rwave=wave[0].tolist()
    allwave=np.array(bwave+gwave+rwave)
    bflux=flux[2].tolist()
    gflux=flux[1].tolist()
    rflux=flux[0].tolist()
    allflux=np.array(bflux+gflux+rflux)
    # get some header values
    objid=hdulist[0].header['objid']
    snr=str(hdulist[0].header['snr'])
    mjd=str(hdulist[0].header['mjd5'])
    exptime=str(hdulist[0].header['exptime'])

    # read the linelist
    linelist=ascii.read(linefile)
    
    # make the plot
    fig=plt.figure(figsize=(18,9))
    matplotlib.rcParams.update({'font.size': 16, 'font.family':'serif'})
    ax=fig.add_subplot(111)
    plt.plot(bwave,bflux,color='black')
    plt.plot(gwave,gflux,color='black')
    plt.plot(rwave,rflux,color='black')
    plt.xlim([15135,16965])
    plt.xlabel(r'Observed Wavelength [$\AA$]')
    plt.ylabel(r'Flux [10$^{-17}$ erg s$^{-1}$ cm$^{2}$ $\AA$]')

    # option to mark airglow lines
    if mark_sky is True:
        bd=np.where(linelist['LABEL']=='airglow')
        bdlines=linelist[bd]
        for i in range(len(bdlines)):
            ax.axvline(bdlines['CENT'][i]+xshift,color='red')

    # option to mark stellar lines
    if mark_lines is True:
        gd=np.where(linelist['LABEL']!='airglow')
        gdlines=linelist[gd]
        for i in range(len(gdlines)):
            line=gdlines['CENT'][i]
            sec=np.where(abs(allwave-line)<0.5)
            arrowstart=np.mean(allflux[sec])+((max(allflux)-min(allflux))*0.10)
            arrowlen=(max(allflux)-min(allflux))*(-0.05)
            arrowheadL=(max(allflux)-min(allflux))*(0.01)
            ax.arrow(line,arrowstart,0,arrowlen,head_width=4,head_length=arrowheadL,color='green')
            txty=np.mean(allflux[sec])+((max(allflux)-min(allflux))*0.11)
            lab=gdlines['LABEL'][i]
            sec=np.where(abs(allwave-line)<0.5)
            poo=arrowstart+(max(allflux)-min(allflux))*(0.05)
            ax.text(line,poo,lab.replace("_"," "),rotation=90,ha='center',va='center',fontsize=12)

    label=file+'    '+objid+'    S/N='+snr+'    exptime='+exptime+' s'
    plt.text(0.5,0.95,label,transform=ax.transAxes,ha='right',fontsize=12)
    plt.tight_layout()

    # define a key press event
    def on_key(event):
        linelist=ascii.read(linefile)
        global key, xdata, ydata
        key=event.key
        xdata=event.xdata
        ydata=event.ydata
        xdata_str="%0.3f" % xdata
        ydata_str="%0.3f" % ydata

        # find the nearest line in the linelist
        dist=abs(xdata-linelist['CENT'])
        x=np.where(dist==min(dist))
        z=x[0].astype('int')
        rest=linelist['CENT'][z]
        rest_str=str(rest[0])
        dif=rest[0]-xdata
        dif_str="%0.3f" % dif
        closeline=linelist['LABEL'][z].tostring()

        if abs(dif) > 8:
            print str(xdata)+r'....No lines within 10 $\AA$'
        else:
            print str(xdata)+'....nearest line= '+closeline+' '+rest_str+', dif= '+dif_str+')'

    cid=fig.canvas.mpl_connect('key_press_event',on_key)

    hdulist.close()

    return

