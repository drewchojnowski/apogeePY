'''
---------------------------------------------------------------------------
---------------------------------------------------------------------------
APOGEESPEC

Basic tools for dealing with APOGEE spectra of emission line stars.
Written by Drew Chojnowski, 12/2015

List of procedures:

1. APLOAD: takes an apVisit spectrum and loads the wavelength and flux arrays
           for the blue (left), green (middle), and red (right) detectors
           into single array properly-ordered arrays.

2. GET_STAR_SPECLIST: given a 2MASS or ABE ID, parses the spectra directory
           for associated apVisits and returns a list of them.

3. APCONTINUUM: apVisit continuum removal.
    3a. NORMALIZE_SPECTRUM: generic interactive continuum normalization. User 
               selects continuum points, program fits splines to the points.

    3b. COMBINE_NORMALIZED_CHIPS: finds apVisit*_normB, apVisit*_normG, and
               apVisit*_normR and combines them into a single normalized
               spectrum file.

4. APCONT_STAR: do continuum removal for all spectra pertaining to a star in
               the catalog file (specify by ABE ID).

5. APSPLOT: IRAF_SPLOT-like routine. Plots the spectrum and allows some user
               interaction, e.g. hit 'x' to print to terminal the nearest 
               linelist line to the cursor position.

6. APSPLOTM: IRAF_SPECPLOT-like routine. Plots multiple spectra; user can 
               specify separation.

7. BROWSER: Plots quantities (e.g. RA vs DEC or J-K vs H) pertaining to stars
               in 'catalog.dat'. Clicking a point in the plot calls APSPLOT,
               and plots the highest signal-to-noise spectrum available.

---------------------------------------------------------------------------
---------------------------------------------------------------------------
'''
specdir='../spectra/'
normspecdir='../normspec/'
datadir='data/'
linefile=datadir+'abe_linelist.dat'
airglowfile=datadir+'combine_airglow_linelists.txt'
starcatalog=datadir+'abe_catalog.dat'

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
GET_STAR_SPECLIST: multiple spectrum plotting program
---------------------------------------------------------------------------
'''
def get_star_speclist(star=None,abeid=None,normspec=False):
    import glob
    from pyfits import getheader
    import numpy as np
    from astropy.io import ascii

    startable=ascii.read(starcatalog)
    abeID=np.array(startable['ID'])
    tmassID=np.array(startable['2MASS'])

    if normspec is False: gdir=specdir
    if normspec is True: gdir=normspecdir
    allspectra=np.array(glob.glob(gdir+'*apV*fits'))
    nspec=len(allspectra)

    stars_all=[]; 
    for i in range(nspec):
        head=getheader(allspectra[i],0)
        stars_all.append(head['objid'])

    if abeid is None:
        gd=np.where(stars_all==star)
        spectra=allspectra[gd]
    else:
        p=np.where(abeID==abeid)
        tm=tmassID[p]
        gd=np.where(stars_all==tm)
        spectra=allspectra[gd]

    return spectra


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
    from pyfits import getheader

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
        head=getheader(spectra[i])
        stars_all.append(head['objid'])
        if quantities[0]=='J-K':
            quant1_all[i]=head['J']-head['K']
        else:
            quant1_all[i]=head[quantities[0]]
        quant2_all[i]=head[quantities[1]]
        snr_all[i]=head['snr']
        hmag_all[i]=head['h']

    stars_all=np.array(stars_all)
    stars,indices=np.unique(stars_all,return_index=True)
    nstars=len(stars)
    quant1=quant1_all[indices]
    quant2=quant2_all[indices]
    hmag=hmag_all[indices]

    fig, ax = plt.subplots(figsize=(12,8))
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

    offset=((max(quant2)-min(quant2))*0.025)
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


'''
---------------------------------------------------------------------------
APSPLOTM: multiple spectrum plotting program
---------------------------------------------------------------------------
'''
def apsplotm(star=None,abeid=None,usenorm=True):
    import glob
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    from astropy.io import ascii
    from astropy.io import fits
    from pyfits import getheader

    def options():
        print('\n')
        print('*********************************************************************')
        print('APSPLOTM KEY PRESS OPTIONS:')
        print('*********************************************************************')
        print('     o:  redisplay the options')
        print('     r:  redraw the plot to default')
        print('     x:  print the nearest line in the linelist')
        print(' f1-f9:  change the separation between spectra')
        print('*********************************************************************')
        print('\n')

    startable=ascii.read(starcatalog)
    abeID=np.array(startable['ID'])
    tmassID=np.array(startable['2MASS'])

    if abeid is None:
        pID=star; sID='none'
    else:
        p=np.where(abeID==abeid)
        pID=tmassID[p]; sID=abeid

    test1=get_star_speclist(star=pID)
    test2=get_star_speclist(star=pID,normspec=True)
    if usenorm is True: 
        if len(test1)>len(test2):
            print("some of the apVisits haven't been normalized... defaulting to originals")
            spectra=test1
            usingnorm=False
        else:
            print("all apVisits have been normalized... using them")
            spectra=test2
            usingnorm=True
    else:
        print("usenorm is false, so defaulting to originals")
        spectra=test1
        usingnorm=False
    nspec=len(spectra)

    options()

    global sep
    sep=0.0

    hjd=np.zeros(nspec); hjd_str=[]
    snr=np.zeros(nspec)
    flux_min=np.zeros(nspec); flux_max=np.zeros(nspec)
    flux_lev=np.zeros(nspec)
    pmf=[]
    for i in range(nspec):
        hdulist=fits.open(spectra[i])
        head=hdulist[0].header
        hjd[i]=head['hjd']
        hjd_str.append(str(hjd[i])[:5])
        snr[i]=head['snr']
        if usingnorm is False:
            wave,flux=apload(spectra[i])
        else:
            wave=hdulist[1].data
            flux=hdulist[2].data
        flux_min[i]=np.min(flux)
        flux_max[i]=np.max(flux)
        wave_max=np.max(wave)
        sec=np.where((wave>16950) & (wave<wave_max))
        flux_lev[i]=np.mean(flux[sec])

        tmp=os.path.split(spectra[i])
        tmp=tmp[len(tmp)-1]
        tmp=tmp.rsplit('-')
        pmf.append(tmp[2]+'-'+tmp[3]+'-'+tmp[4][:3])
    pmf=np.array(pmf)
    hjd_str=np.array(hjd_str)
    order=hjd_str.argsort()
    spectra=spectra[order]
    hjd=hjd[order]
    hjd_str=hjd_str[order]
    snr=snr[order]
    pmf=pmf[order]

    maxflux=np.max(flux_max)
    minflux=np.min(flux_min)
    span=(maxflux+minflux)

    cmap=matplotlib.cm.gnuplot
    fig=plt.figure(figsize=(14,8))
    fig.canvas.set_window_title('2MASS ID = '+pID[0]+',    ABE ID = '+sID+',    nvisits = '+str(nspec)+' ')
    fig.subplots_adjust(left=0.08,bottom=0.08,right=0.89,top=0.97)
    matplotlib.rcParams.update({'font.size': 14, 'font.family':'serif'})
    ax=fig.gca()
    ax.grid(True)

    ax.set_xlim([15130,16965])
    ax.set_xlabel(r'Observed Wavelength [$\AA$]')
    ax.set_ylabel(r'Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA$]')

    def plotspectra(norm,sep=0.0,more_sep=None):
        if more_sep is not None: sep=span*0.01*more_sep
        for i in range(nspec):
            hdulist=fits.open(spectra[i])
            if norm is False:
                wave,flux=apload(spectra[i])
            else:
                wave=hdulist[1].data
                flux=hdulist[2].data
            hdulist.close()

            p=ax.plot(wave,flux+sep*i,color=cmap(i/float(nspec)))
            ax.text(np.max(wave)+20,flux_lev[i]+sep*i,pmf[i],color=p[0].get_color(),ha='left',fontsize=12)
            ax.grid(True)
            
            ax.set_xlabel(r'Observed Wavelength [$\AA$]')
            ax.set_ylabel(r'Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA$]')
        return 

    def ontype(event):
        global key, xdata, ydata
        key=event.key; xdata=event.xdata; ydata=event.ydata

        # when the user hits 'x': print the nearest line in linelist
        if event.key=='x':
            linelist=ascii.read(linefile)

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
        elif (len(key)==2) and (key[:1]=='f'):
            ax.cla()
            newsep=float(key[1:])
            xmin,xmax=ax.get_xlim()
            plotspectra(usingnorm,more_sep=newsep)
            ymin,ymax=ax.get_ylim()
            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])
            plt.draw()
        elif key=='r':
            shiftcount=0
            newsep=0.0
            sep=0.0
            ax.cla()
            plotspectra(usingnorm,more_sep=newsep)
            ax.set_xlim([15130,16965])
            plt.draw()
        # when the user hits 'o': redisplay the options
        elif event.key=='o':
            options()

    plotspectra(usingnorm,sep)
#    plt.subplots_adjust(left=0.08,bottom=0.08,right=0.89,top=0.97)

    fig.canvas.mpl_connect('key_press_event',ontype)

#    plt.show()

    return

'''
---------------------------------------------------------------------------
APSPLOT: spectrum plotting program
---------------------------------------------------------------------------
'''
def apsplot(infile,mark_airglow=True,mark_lines=True,xshift=0.0,winwidth=1.0,airglow_width=2.0,airglow_cut=100.0,xr=None):
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    from astropy.io import ascii
    from astropy.io import fits
    import os
    import glob

    def options():
        print('\n')
        print('*********************************************************************')
        print('APSPLOT KEY PRESS OPTIONS:')
        print('*********************************************************************')
        print('           o:  redisplay the options')
        print('           x:  print to terminal the nearest linelist entry')
        print('           n:  display normalized spectrum in new window')
        print('           h:  print the FITS header to terminal')
        print('*********************************************************************')
        print('APSPLOT BUTTON PRESS OPTIONS:')
        print('*********************************************************************')
        print('  left-click:  plot a square symbol at cursor position')
        print(' right-click:  delete the square symbol')
        print('       enter:  write to a text file the square symbol coordinates')
        print('*********************************************************************')
        print('\n')

    options()

    hdulist=fits.open(infile)
    if len(hdulist)>10:
        wave,flux=apload(infile)
    else:
        wave=hdulist[1].data
        flux=hdulist[2].data

    # option to add on an xshift, in angstroms
    wave=wave+xshift

    # get some header values from the FITS file
    hdulist=fits.open(infile)
    head=hdulist[0].header
    objid=head['objid']
    snr=str(head['snr'])
    mjd=str(head['mjd5'])
    exptime=str(head['exptime'])

    # make the plot
    fig=plt.figure(figsize=(16,8))
    tmp=os.path.split(infile); tmp=tmp[len(tmp)-1]
    fig.canvas.set_window_title('file='+tmp+',     star='+objid+',     S/N='+snr+',     exptime='+exptime+' s,     xshift='+str(xshift))
    matplotlib.rcParams.update({'font.size': 14, 'font.family':'serif'})
    ax=fig.add_subplot(111)
    ax.plot(wave,flux,color='black')
    ymin,ymax=ax.get_ylim()
    if (ymax-ymin)<1.2: ax.set_ylim([0.4,1.6])
    ax.grid(True)
    if xr is not None:
        ax.set_xlim(xr)
    else:
        ax.set_xlim([15130,16965])
#    if len(hdulist)<10: plt.ylim([0.5,1.5])
    ax.set_xlabel(r'Observed Wavelength [$\AA$]')
#    ax.set_ylabel(r'Flux [10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]')
    ax.set_ylabel(r'Flux')
    plt.tight_layout()

    # option to mark airglow lines
    if mark_airglow is True:
	airglow=ascii.read(airglowfile)
	gd=np.where(airglow['EMISSION']>airglow_cut)
	airglow=airglow[gd]
        for i in range(len(airglow)):
            pos=airglow['WAVE'][i]+xshift
            sec=np.where(wave>(pos-airglow_width/2.0))
            wtmp=wave[sec]; ftmp=flux[sec]
            sec=np.where(wtmp<(pos+airglow_width/2.0))
            wtmp=wtmp[sec]; ftmp=ftmp[sec]
            ax.plot(wtmp,ftmp,color='red')

    # option to mark stellar lines
    if mark_lines is True:
    	linelist=ascii.read(linefile)
        for i in range(len(linelist)):
            line=linelist['CENT'][i]
            lab=linelist['LABEL'][i]
            if lab[0:1]=='H': labcol='blue'
            if lab[0:1]!='H': labcol='green'
            sec=np.where(abs(wave-line)<0.5)
            arrowstart=np.mean(flux[sec])+((max(flux)-min(flux))*0.10)
            arrowlen=(max(flux)-min(flux))*(-0.05)
            arrowheadL=(max(flux)-min(flux))*(0.01)
            sec=np.where(abs(wave-line)<0.5)
            txty=arrowstart+(max(flux)-min(flux))*(0.015)
            if txty<1: 
                ax.arrow(line,1.09,0,-0.05,head_width=2,head_length=arrowheadL,color=labcol)
                ax.text(line,1.1,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)
            else:
                ax.arrow(line,arrowstart,0,arrowlen,head_width=3,head_length=arrowheadL,color=labcol)
                ax.text(line,txty,lab.replace("_"," "),rotation=90,ha='center',va='bottom',fontsize=9,color=labcol)

    # key press event handler: tells you the nearest spectral line
#    def apsplot_on_key(event):


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
                ax.plot(event.xdata,y,'ys',ms=10,picker=5,label='cont_pnt')
            plt.draw()

    def apsplot_onpick(event):
        # when the user right clicks on a continuum point, remove it
        if event.mouseevent.button==3:
            if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
                event.artist.remove()

    def ontype(event):
        # when the user hits enter: # write a text file with x,y data
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
            data=np.array([x,y])
            np.savetxt('apslot_output.log',data.transpose(),fmt='%.5f')

        # when the user hits 'x': print the nearest line in linelist
        elif event.key=='x':
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

        # when the user hits 'n': show the normalized spectrum if it exists
        elif event.key=='n':
            tmp=os.path.split(infile)
            tmp=tmp[len(tmp)-1]
            tmp1=tmp.rsplit('.fits')
            tmp=tmp1[0]+'_norm.fits'
            normfile=normspecdir+tmp
            search=glob.glob(normfile)
            if len(search)==1:
                apsplot(normfile)
            else:
                print('normalized spectrum does not exist.')

        # when the user hits 'o': redisplay the options
        elif event.key=='o':
            options()

        # when the user hits 'h': print the header to terminal
        elif event.key=='h':
            print('\n')
            print('header of '+infile)
            print(head.cards)
            print('\n')

    fig.canvas.mpl_connect('key_press_event',ontype)
    fig.canvas.mpl_connect('button_press_event',apsplot_onclick)
    fig.canvas.mpl_connect('pick_event',apsplot_onpick)

    hdulist.close()
    plt.show()

#    answer=raw_input('hit "q" to quit: ')
#    if answer=='q': plt.close()
    
    return


'''
---------------------------------------------------------------------------
APCONT_STAR: wrapper for normalize_spectrum & combine_normalized_chips
---------------------------------------------------------------------------
'''
def apcont_star(abeid=None,redo=True):
    import glob
    import os
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib
    from astropy.io import ascii
    from astropy.io import fits
    from pyfits import getheader
    import sys

    startable=ascii.read(starcatalog)
    abeID=np.array(startable['ID'])
    tmassID=np.array(startable['2MASS'])

    if abeid is None:
        pID=star; sID='none'
    else:
        p=np.where(abeID==abeid)
        pID=tmassID[p]; sID=abeid

    test1=get_star_speclist(star=pID)
    test2=get_star_speclist(star=pID,norm=True)
    if len(test1)==len(test2):
        if redo is False:
            print("all of the apVisit have been normalized. Set the redo key to overwrite them.")
            sys.exit()

    spectra=test1

    for i in range(len(spectra)):
        apcontinuum(spectra[i])

    print('All apVisits for ABE-'+abeid+' have been continuum-normalized.')

    return

'''
---------------------------------------------------------------------------
APCONTINUUM: wrapper for normalize_spectrum & combine_normalized_chips
---------------------------------------------------------------------------
'''
def apcontinuum(infile,combine=True):
    import numpy as np
    from astropy.io import fits

    print('-'*80)
    print('APCONTINUUM OPTIONS')
    print('-'*80)
    print(' left-click: mark a continuum point on the continuum')
    print('          c: mark a continuum point at cursor y-position')
    print('right-click: delete a continuum point')
    print('      enter: show continuum fit')
    print('          n: apply continuum fit')
    print('          r: start over')
    print('          w: save continuum fit')
    print('     escape: escape')
    print('-'*80)
    print("ignore the below warning message... I don't know how to suppress it")
    print('-'*80)

    hdulist=fits.open(infile)
    wave=hdulist[4].data
    flux=hdulist[1].data
#    print('-'*80)
#    print('(1) continuum-normalizing the blue chip of '+infile+'...')
    x=normalize_spectrum(np.array(wave[2]),np.array(flux[2]),infile,flabel='normB')
#    print('(2) continuum-normalizing the green chip of '+infile+'...')
    x=normalize_spectrum(np.array(wave[1]),np.array(flux[1]),infile,flabel='normG')
#    print('(3) continuum-normalizing the red chip of '+infile+'...')
    x=normalize_spectrum(np.array(wave[0]),np.array(flux[0]),infile,flabel='normR')
#    print('-'*80)

    if combine is True:
        print('-'*80)
        print('combining the normalized chips of '+infile+'...')
        combine_normalized_chips(infile)
        print('-'*80)
        print('Done')

    return

'''
---------------------------------------------------------------------------
COMBINE_NORMALIZED_CHIPS: combine normalized chips a into single spectrum
---------------------------------------------------------------------------
'''
def combine_normalized_chips(infile):
    import numpy as np
    from astropy.io import fits
    import glob
    import time
    import pyfits
    import os
    import sys

    # chip gap edges
    gap1_Bedge=15807.0
    gap1_Redge=15859.0
    gap2_Bedge=16431.0
    gap2_Redge=16474.0

    # find the normalized chip files
    chipBfile=glob.glob(infile.rsplit('.fits')[0]+'_normB.fits')
    chipGfile=glob.glob(infile.rsplit('.fits')[0]+'_normG.fits')
    chipRfile=glob.glob(infile.rsplit('.fits')[0]+'_normR.fits')
    chips=[chipBfile,chipGfile,chipRfile]
    lenchips=len(chipBfile)+len(chipGfile)+len(chipRfile)
    if lenchips==3:
        print('normalized B,G,R chips found')
    else:
        if len(chipBfile)==0: print('ERROR: normalized B chip not found!')
        if len(chipGfile)==0: print('ERROR: normalized G chip not found!')
        if len(chipRfile)==0: print('ERROR: normalized R chip not found!')
        print('Program execution halted.')
        sys.exit()

    # get single lists of wavelength, flux, and normalized flux
    allw=[]; allf=[]; allfn=[]
    for i in range(0,3):
        tmp=chips[i]
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
    orighdu=fits.open(infile)
    orighead=orighdu[0].header
    barycor=orighead['bc']
    allw=allw-((barycor/299792.458)+1)
    origcards=orighead.cards

    # create the new FITS file, copying the header from the original file
    # and updating it.
    outfile=infile.rsplit('.fits')[0]+'_norm.fits'
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
NORMALIZE_SPECTRUM: do chip-by-chip continuum normalization
---------------------------------------------------------------------------
'''
def normalize_spectrum(wave,flux,infile,flabel=None,window=7.0,splineord=3,contPoints2file=False):
    from scipy.interpolate import splrep,splev
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits
    import pyfits
    import os
    import sys



    continuum=None

    hdulist=fits.open(infile)
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
    tmp=os.path.split(infile); tmp=tmp[len(tmp)-1]
    fig.canvas.set_window_title('NORMALIZE_SPECTRUM:     file = '+tmp+',     objid = '+objid+',     S/N = '+snr+',      UT-MID = '+utmid)
    plt.title(plab,color=col)
    plt.tight_layout()
    fig.show()

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
                ax.grid(True)
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
        if event.key=='c':
            y = event.ydata
            ax.plot(event.xdata,y,'rs',ms=10,picker=5,label='cont_pnt')
            ax.grid(True)
            plt.draw()


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
            spline = splrep(x,y,k=splineord)
            global continuum
            continuum = splev(wave,spline)
            cfunc = lambda w: splev(w, spline)
            plt.plot(wave,continuum,'r-',lw=2,label='continuum')

            if contPoints2file is True:
                outfile=infile.rsplit('.fits')[0]+'_X'+flabel+'.txt'
                np.savetxt(outfile,cont_pnt_coord,fmt='%.5f')

        # when the user hits 'n' and a spline-continuum is fitted, normalise the
        # spectrum
        elif event.key=='n':
            if continuum is not None:
                ax.cla()
                ax.plot(wave,flux/continuum,'k-',label='normalised')
#                ax.set_ylim([0.7,1.3])
                ax.set_title('"w" to write to file; "r" to start over',color=col)
                ax.grid(True)

        # when the user hits 'r': clear the axes and plot the original spectrum
        elif event.key=='r':
            continuum = None
            ax.cla()
            ax.plot(wave,flux,'k-')
            ax.grid(True)

        # when the user hits 'w': if the normalised spectrum exists, write it to a file.
        elif event.key=='w':
            for artist in ax.get_children():
                if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                    outfile=infile.rsplit('.fits')[0]+'_'+flabel+'.fits'

                    ax.cla()
                    ax.plot(wave,flux/continuum,'k-',label='normalised')
#                    ax.set_ylim([0.7,1.3])
                    ax.set_title('*** file written: '+outfile,color=col)
                    ax.grid(True)
                    data=np.array(artist.get_data())

                    hdu=pyfits.PrimaryHDU()
                    new_hdu=pyfits.HDUList([hdu])
                    new_hdu.append(pyfits.ImageHDU(wave))
                    new_hdu.append(pyfits.ImageHDU(flux/continuum))
                    new_hdu.append(pyfits.ImageHDU(flux))
                    new_hdu.writeto(outfile,clobber=True)
#                    print('Saved to file: '+outfile)
                    break
            fig.canvas.stop_event_loop()

        # if the user hits 'escape': exit the program.
        elif event.key=='escape':
            fig.canvas.stop_event_loop()
            sys.exit()
            
        plt.draw()

    def mark() :
        ''' blocking routine to wait for keypress event and return data'''
        # start the blocking, wait indefinitely (-1) for an event
        fig.canvas.start_event_loop(-1)
        # once we get the event, return nothing
        return 0

    fig.canvas.mpl_connect('key_press_event',ontype)
    fig.canvas.mpl_connect('button_press_event',onclick)
    fig.canvas.mpl_connect('pick_event',onpick)

    x=mark()
    
    plt.close()

    return 


