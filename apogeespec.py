'''
---------------------------------------------------------------------------
---------------------------------------------------------------------------
APOGEESPEC

Basic tools for dealing with APOGEE spectra of emission line stars.
Written by Drew Chojnowski, 12/2015
---------------------------------------------------------------------------
---------------------------------------------------------------------------
'''
specdir='../spectra/'
linefile='linelist.dat'
starcatalog='catalog.dat'

'''
---------------------------------------------------------------------------
APLOAD: load the apVisit wave & flux into a single array
---------------------------------------------------------------------------
'''
def apload(file):
    import numpy as np
    from astropy.io import fits
    # read the FITS file
    hdulist=fits.open(file)
    wave=hdulist[4].data
    flux=hdulist[1].data

    # get single arrays of wavelength and flux values, ordered properly
    rw=np.array(wave[0]); rw=rw[::-1]; rw=rw.tolist()
    rf=flux[0]; rf=rf[::-1]; rf=rf.tolist()
    gw=np.array(wave[1]); gw=gw[::-1]; gw=gw.tolist()
    gf=flux[1]; gf=gf[::-1]; gf=gf.tolist()
    bw=np.array(wave[2]); bw=bw[::-1]; bw=bw.tolist()
    bf=flux[2]; bf=bf[::-1]; bf=bf.tolist()

    allwave=np.array(bw+gw+rw); allflux=np.array(bf+gf+rf)

    return allwave,allflux

'''
---------------------------------------------------------------------------
BROWSER: interactive plot (button press calls apsplot)
---------------------------------------------------------------------------
'''
def browser(quantities=None):
    import numpy as np
    import matplotlib.pyplot as plt
    import glob
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from astropy.io import ascii
    from astropy.io import fits

    print('gathering browser data... please wait')

    if quantities is None:
        quantities=['RA','DEC']

    startable=ascii.read(starcatalog)
    abeID=np.array(startable['ID'])
    tmassID=np.array(startable['2MASS'])

    spectra=np.array(glob.glob(specdir+'*apV*fits'))
    nspec=len(spectra)
    stars_all=[]
    quant1_all=np.zeros(nspec); quant2_all=np.zeros(nspec)
    snr_all=np.zeros(nspec); hmag_all=np.zeros(nspec)
    for i in range(nspec):
        hdulist=fits.open(spectra[i])
        stars_all.append(hdulist[0].header['objid'])
        if quantities[0]=='J-K':
            quant1_all[i]=hdulist[0].header['J']-hdulist[0].header['K']
        else:
            quant1_all[i]=hdulist[0].header[quantities[0]]
        quant2_all[i]=hdulist[0].header[quantities[1]]
        snr_all[i]=hdulist[0].header['snr']
        hmag_all[i]=hdulist[0].header['h']

    stars_all=np.array(stars_all)
    stars,indices=np.unique(stars_all,return_index=True)
    nstars=len(stars)
    quant1=quant1_all[indices]
    quant2=quant2_all[indices]
    hmag=hmag_all[indices]

    fig, ax = plt.subplots(figsize=(12,7))
    ax.set_title('click on point to show highest S/N spectrum')
    line,=ax.plot(quant1,quant2,'w.',picker=5,markersize=0.1,linewidth=0)  # 5 points tolerance
    plt.xlabel(quantities[0]);  plt.ylabel(quantities[1])
    if quantities is None:
        plt.xlim([-5,365])
        plt.xticks(np.arange(0,361,60))
    scat=plt.scatter(quant1,quant2,c=hmag,s=160,cmap='rainbow',alpha=0.5)

    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", "5%", pad="3%")
    cbar=plt.colorbar(scat,cax=cax)
    cbar.set_label('H mag')

    offset=((max(quant2)-min(quant2))*0.05)
    for i in range(nstars):
        p=np.where(stars[i]==tmassID)
        ax.text(quant1[i],quant2[i]+offset,abeID[p][0],fontsize=9,ha='center')
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



def apcompare(files):
    for i in range(len(files)):
        x=apsplot(files[i])

    return

'''
---------------------------------------------------------------------------
APSPLOT: spectrum plotting program
---------------------------------------------------------------------------
'''
def apsplot(file,mark_sky=False,mark_lines=True,xshift=0.0):
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    from astropy.io import ascii
    from astropy.io import fits
    import os

    linelist=ascii.read(linefile)    # read the linelist

    hdulist=fits.open(file)
    if len(hdulist)>10:
        wave,flux=apload(file)
    else:
        wave=hdulist[1].data
        flux=hdulist[2].data

    # option to add on an xshift, in angstroms
    wave=wave+xshift

    # get some header values from the FITS file
    hdulist=fits.open(file)
    objid=hdulist[0].header['objid']
    snr=str(hdulist[0].header['snr'])
    mjd=str(hdulist[0].header['mjd5'])
    exptime=str(hdulist[0].header['exptime'])

    # make the plot
    fig=plt.figure(figsize=(16,8))
    tmp=os.path.split(file); tmp=tmp[len(tmp)-1]
    fig.canvas.set_window_title('file='+tmp+',     star='+objid+',     S/N='+snr+',     exptime='+exptime+' s')
    matplotlib.rcParams.update({'font.size': 14, 'font.family':'serif'})
    ax=fig.add_subplot(111)
    ax.plot(wave,flux,color='black')
    ax.grid(True)
    plt.xlim([15135,16965])
    if len(hdulist)<10: plt.ylim([0.5,1.5])
    plt.xlabel(r'Observed Wavelength [$\AA$]')
    plt.ylabel(r'Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA$]')
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
            if txty<1: 
                ax.arrow(line,1.09,0,-0.05,head_width=3,head_length=arrowheadL,color=labcol)
                ax.text(line,1.1,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)
            else:
                ax.arrow(line,arrowstart,0,arrowlen,head_width=3,head_length=arrowheadL,color=labcol)
                ax.text(line,txty,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)

    # key press event handler: tells you the nearest spectral line
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
            print(xdata_str+r'....No lines within 10 $\AA$')
        else:
            bla=xdata_str+'....nearest line= '+closeline+' '+rest_str+', dif= '+dif_str+')'
            print(bla)

    winwidth=10.0
    def apsplot_onclick(event):
        # when none of the toolbar buttons is activated and the user clicks in the
        # plot somewhere, compute the median value of the spectrum in a 10angstrom
        # window around the x-coordinate of the clicked point. The y coordinate
        # of the clicked point is not important. Make sure the continuum points
        # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
        toolbar=plt.get_current_fig_manager().toolbar
        if (event.xdata>np.min(wave)) & (event.xdata<np.max(wave)):
            if event.button==1 and toolbar.mode=='':
                window=((event.xdata-winwidth)<=wave) & (wave<=(event.xdata+winwidth))
                y=np.median(flux[window])
                ax.plot(event.xdata,y,'rs',ms=10,picker=5,label='cont_pnt')
            plt.draw()

    def apsplot_onpick(event):
        # when the user right clicks on a continuum point, remove it
        if event.mouseevent.button==3:
            if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
                event.artist.remove()

    fig.canvas.mpl_connect('key_press_event',apsplot_on_key)
    fig.canvas.mpl_connect('button_press_event',apsplot_onclick)
    fig.canvas.mpl_connect('pick_event',apsplot_onpick)

    hdulist.close()
    plt.show()

    return

'''
---------------------------------------------------------------------------
APCONTINUUM: wrapper for continuum normalization (doesn't work)
---------------------------------------------------------------------------
'''
def apcontinuum(file):
    run_normalize_chips(file)
    combine_normalized_chips(file)

    return

'''
---------------------------------------------------------------------------
COMBINE_NORMALIZED_CHIPS: combine normalized chip files into single normalized spectrum
---------------------------------------------------------------------------
'''
def combine_normalized_chips(file):
    import numpy as np
    from astropy.io import fits
    import glob
    import time
    import pyfits
    import os

    # chip gap edges
    gap1_Bedge=15807.0
    gap1_Redge=15859.0
    gap2_Bedge=16431.0
    gap2_Redge=16474.0

    # find the normalized chip files
    chipBfile=glob.glob(file.rsplit('.fits')[0]+'_normB.fits')
    chipGfile=glob.glob(file.rsplit('.fits')[0]+'_normG.fits')
    chipRfile=glob.glob(file.rsplit('.fits')[0]+'_normR.fits')
    if len(chipBfile)==0: print('normalized B chip not found!')
    if len(chipGfile)==0: print('normalized G chip not found!')
    if len(chipRfile)==0: print('normalized R chip not found!')
    chips=[chipBfile,chipGfile,chipRfile]
    if len(chips)==3: print('normalized B,G,R chips found')

    # get single lists of wavelength, flux, and normalized flux
    allw=[]; allf=[]; allfn=[]
    for i in range(0,3):
        tmp=chips[i]
        print(tmp[0])
        data=fits.open(tmp[0])
        w=np.array(data[1].data);    w=w[::-1];   w=w.tolist()
        fn=np.array(data[2].data); fn=fn[::-1]; fn=fn.tolist()
        f=np.array(data[3].data);    f=f[::-1];   f=f.tolist()
        allw=allw+w
        allfn=allfn+fn
        allf=allf+f
        os.remove(tmp[0])

    # convert the lists to array
    allw=np.array(allw)
    allfn=np.array(allfn)
    allf=np.array(allf)

    # set the chip gap normalized flux to 1
    gap1=np.where((allw>=gap1_Bedge) & (allw<=gap1_Redge))
    allfn[gap1]=1.0
    gap2=np.where((allw>=gap2_Bedge) & (allw<=gap2_Redge))
    allfn[gap2]=1.0

    # get the original apVisit header
    orighdu=fits.open(file)
    orighead=orighdu[0].header
    origcards=orighead.cards

    # create the new FITS file, copying the header from the original file
    # and updating it.
    outfile=file.rsplit('.fits')[0]+'_norm.fits'
    hdu=pyfits.PrimaryHDU()
    head=hdu.header                    
    for i in range(len(origcards)):
        head.set(origcards[i][0],origcards[i][1],comment=origcards[i][2])
    head.append('HISTORY','continuum fit by Drew on '+time.strftime("%c"),end=True)
    new_hdu=pyfits.HDUList([hdu])
    new_hdu.append(pyfits.ImageHDU(allw))
    new_hdu.append(pyfits.ImageHDU(allfn))
    new_hdu.append(pyfits.ImageHDU(allf))
    new_hdu.writeto(outfile,clobber=True)
    print('Saved to file: '+outfile)
        
    return 

'''
---------------------------------------------------------------------------
RUN_NORMALIZE_CHIPS: wrapper for normalize_chips
---------------------------------------------------------------------------
'''
def run_normalize_chips(file):
    import numpy as np
    from astropy.io import fits

    hdulist=fits.open(file)
    wave=hdulist[4].data
    flux=hdulist[1].data
    x=normalize_chips(np.array(wave[2]),np.array(flux[2]),file,flabel='normB')
    x=normalize_chips(np.array(wave[1]),np.array(flux[1]),file,flabel='normG')
    x=normalize_chips(np.array(wave[0]),np.array(flux[0]),file,flabel='normR')

    return

'''
---------------------------------------------------------------------------
NORMALIZE_CHIPS: do chip-by-chip continuum normalization
---------------------------------------------------------------------------
'''
def normalize_chips(wave,flux,file,flabel=None,window=10.0):
    from scipy.interpolate import splrep,splev
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits
    import pyfits
    import os

    continuum=None

    hdulist=fits.open(file)
    objid=hdulist[0].header['objid']
    snr=str(hdulist[0].header['snr'])
    mjd=str(hdulist[0].header['mjd5'])
    exptime=str(hdulist[0].header['exptime'])
    utmid=hdulist[0].header['UT-MID']

    fig=plt.figure(figsize=(16,8))
    ax = plt.gca()
    winwidth=window/2.0
    ax.plot(wave,flux,'k-',label='spectrum')
    ax.grid(True)
    ax.set_xlabel(r'Wavelength')
    ax.set_ylabel(r'Flux')
    if flabel=='normB':
        plab='blue detector'; col='blue'
    elif flabel=='normG':
        plab='green detector'; col='green'
    elif flabel=='normR':
        plab='red detector'; col='red'
    else:
        plab=' '; col='black'
    tmp=os.path.split(file); tmp=tmp[len(tmp)-1]
    fig.canvas.set_window_title('NORMALIZE_CHIPS:     file = '+tmp+',     objid = '+objid+',     S/N = '+snr+',      UT-MID = '+utmid)
    plt.title(plab,color=col)
    plt.tight_layout()

    def onclick(event):
        # when none of the toolbar buttons is activated and the user clicks in the
        # plot somewhere, compute the median value of the spectrum in a 10angstrom
        # window around the x-coordinate of the clicked point. The y coordinate
        # of the clicked point is not important. Make sure the continuum points
        # `feel` it when it gets clicked, set the `feel-radius` (picker) to 5 points
        toolbar=plt.get_current_fig_manager().toolbar
        if (event.xdata>np.min(wave)) & (event.xdata<np.max(wave)):
            if event.button==1 and toolbar.mode=='':
                window = ((event.xdata-winwidth)<=wave) & (wave<=(event.xdata+winwidth))
                y = np.median(flux[window])
                ax.plot(event.xdata,y,'rs',ms=10,picker=5,label='cont_pnt')
            plt.draw()

    def onpick(event):
        # when the user right clicks on a continuum point, remove it
        if event.mouseevent.button==3:
            if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
                event.artist.remove()

    def ontype(event):
        # when the user hits enter:
        # 1. Cycle through the artists in the current axes. If it is a continuum
        #    point, remember its coordinates. If it is the fitted continuum from the
        #    previous step, remove it
        # 2. sort the continuum-point-array according to the x-values
        # 3. fit a spline and evaluate it in the wavelength points
        # 4. plot the continuum
        if event.key=='enter':
            cont_pnt_coord = []
            for artist in ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                    cont_pnt_coord.append(artist.get_data())
                elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                    artist.remove()
            cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
            sort_array = np.argsort(cont_pnt_coord[:,0])
            x,y = cont_pnt_coord[sort_array].T
            spline = splrep(x,y,k=2)
            global continuum
            continuum = splev(wave,spline)
            cfunc = lambda w: splev(w, spline)
            plt.plot(wave,continuum,'r-',lw=2,label='continuum')

        # when the user hits 'n' and a spline-continuum is fitted, normalise the
        # spectrum
        elif event.key=='n':
            if continuum is not None:
                ax.cla()
                ax.plot(wave,flux/continuum,'k-',label='normalised')
                ax.set_ylim([0.5,1.5])
                ax.set_title('"w" to write to file; "r" to start over',color=col)

        # when the user hits 'r': clear the axes and plot the original spectrum
        elif event.key=='r':
            continuum = None
            ax.cla()
            ax.plot(wave,flux,'k-')

        # when the user hits 'w': if the normalised spectrum exists, write it to a
        # file.
        elif event.key=='w':
            for artist in ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                    outfile=file.rsplit('.fits')[0]+'_'+flabel+'.fits'

                    ax.cla()
                    ax.plot(wave,flux/continuum,'k-',label='normalised')
                    ax.set_ylim([0.5,1.5])
                    ax.set_title('file written: '+outfile,color=col)
                    data=np.array(artist.get_data())
#                    if flabel is None: flabel='_norm'

                    hdu=pyfits.PrimaryHDU()
                    new_hdu=pyfits.HDUList([hdu])
                    new_hdu.append(pyfits.ImageHDU(wave))
                    new_hdu.append(pyfits.ImageHDU(flux/continuum))
                    new_hdu.append(pyfits.ImageHDU(flux))
                    new_hdu.writeto(outfile,clobber=True)
                    print('Saved to file: '+outfile)
                    break
        plt.draw()

    fig.canvas.mpl_connect('key_press_event',ontype)
    fig.canvas.mpl_connect('button_press_event',onclick)
    fig.canvas.mpl_connect('pick_event',onpick)

    return continuum





