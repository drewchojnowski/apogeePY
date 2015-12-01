import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import ascii
from scipy.interpolate import interp1d

def splot(file,linefile='linelist.dat',mark_sky=False,mark_lines=True,xshift=0.0,do_cont=True,do_cont_chips=False,cont_ord=5):
    # read the linelist
    linelist=ascii.read(linefile)

    # read the FITS file
    hdulist=fits.open(file)
    wave=hdulist[4].data+xshift
    flux=hdulist[1].data

    # get some header values
    objid=hdulist[0].header['objid']
    snr=str(hdulist[0].header['snr'])
    mjd=str(hdulist[0].header['mjd5'])
    exptime=str(hdulist[0].header['exptime'])

    # get single arrays of wavelength and flux values, ordered properly
    rw=np.array(wave[0]); rw=rw[::-1]; rw=rw.tolist(); rf=flux[0]; rf=rf[::-1]; rf=rf.tolist()
    gw=np.array(wave[1]); gw=gw[::-1]; gw=gw.tolist(); gf=flux[1]; gf=gf[::-1]; gf=gf.tolist()
    bw=np.array(wave[2]); bw=bw[::-1]; bw=bw.tolist(); bf=flux[2]; bf=bf[::-1]; bf=bf.tolist()
    allwave=np.array(bw+gw+rw); allflux=np.array(bf+gf+rf)

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

            allflux=np.array(bf+gf+rf)
        # else normalize the whole thing at once
        else:
            coeff=np.polyfit(allwave,allflux,cont_ord)
            poly=np.poly1d(coeff)
            cont=poly(allwave)
            allflux=allflux/cont

    # make the plot
    fig=plt.figure(figsize=(18,9))
    matplotlib.rcParams.update({'font.size': 16, 'font.family':'serif'})
    ax=fig.add_subplot(111)
    plt.plot(allwave,allflux,color='black')
    plt.xlim([15135,16965])
    plt.xlabel(r'Observed Wavelength [$\AA$]')
    plt.ylabel(r'Flux [10$^{-17}$ erg s$^{-1}$ cm$^{2}$ $\AA$]')
    plt.title(file+'    '+objid+'    S/N='+snr+'    exptime='+exptime+' s')
    plt.tight_layout()

    # option to mark airglow lines
    if mark_sky is True:
        bd=np.where(linelist['LABEL']=='airglow')
        bdlines=linelist[bd]
        for i in range(len(bdlines)):
            pos=bdlines['CENT'][i]+xshift
            sec=np.where(allwave>(pos-2))
            wtmp=allwave[sec]; ftmp=allflux[sec]
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
            sec=np.where(abs(allwave-line)<0.5)
            arrowstart=np.mean(allflux[sec])+((max(allflux)-min(allflux))*0.10)
            arrowlen=(max(allflux)-min(allflux))*(-0.05)
            arrowheadL=(max(allflux)-min(allflux))*(0.01)
            sec=np.where(abs(allwave-line)<0.5)
            txty=arrowstart+(max(allflux)-min(allflux))*(0.015)
            if do_cont is True and txty<1: 
                ax.arrow(line,1.09,0,-0.05,head_width=3,head_length=arrowheadL,color=labcol)
                ax.text(line,1.1,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)
            else:
                ax.arrow(line,arrowstart,0,arrowlen,head_width=3,head_length=arrowheadL,color=labcol)
                ax.text(line,txty,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)

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
        rest_str="%0.3f" % rest[0]
        dif=rest[0]-xdata
        dif_str="%0.3f" % dif
        closeline=linelist['LABEL'][z][0].tolist()

        if abs(dif) > 8:
            print(xdata_str+r'....No lines within 10 $\AA$')
        else:
            bla=xdata_str+'....nearest line= '+closeline+' '+rest_str+', dif= '+dif_str+')'
            print(bla)

    cid=fig.canvas.mpl_connect('key_press_event',on_key)

    hdulist.close()

    return linelist['LABEL']




