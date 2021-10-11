import datetime
import os
import sys

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)

    if cleanup:
        for ttsuffix in ('.tt0', '.tt1', 'tt2'):
            for suffix in ('pb{tt}', 'weight', 'sumwt{tt}', 'psf{tt}',
                           'model{tt}', 'mask', 'image{tt}', 'residual{tt}',
                           'alpha', ):
                os.system('rm -rf {0}.{1}'.format(myimagebase, suffix).format(tt=ttsuffix))

# Scripit to reduce all cont data.
# Make sure to change values for imsize

#====================================
fld = 'N' # change to 'M' as needed
conf = 'B3' # change to 'B6' as needed
step = 0 # change to 1/2/3. 
# 0 - initial clean
# 1 - p,inf cal + clean
# 2 - p,int cal + clean
# 3 - ap,15s cal + clean
#====================================


fldconf = fld+conf

if fld == 'N':
    mask = 'SgrB2N_regions.crtf'
    
if fld == 'M':
    mask = 'SgrB2M_regions.crtf'

if conf == 'B3':
    imsize = [6272, 6272] #adjust?
    cell = ['0.018arcsec']
    threshold='0.1mJy'

if conf == 'B6':
    imsize = [4480,4480]
    cell =  ['0.01arcsec']
    threshold='3.0mJy'

if fldconf == 'NB3':
    mslist = ['N_cont1.ms','N_cont2.ms']
elif fldconf == 'MB3':
    mslist = ['M_cont1.ms','M_cont2.ms']
elif fldcong == 'NB6':
    mslist = ['N_B6_cont1.ms','N_B6_cont2.ms']
else: mslist = ['M_B6_cont1.ms','M_B6_cont2.ms']


cal = 'cal'+str(step)

if step == 1:
    solint='inf'
    calmode='p'
elif step == 2:
    solint='int'
    calmode='p'
elif step == 3:
    solint='15s'
    calmode='ap'

mslist = [s + '_cal' + str(step) for s in mslist]
print("The folders are: " + str(mslist))

source_spws = (0,1,2,3)

spwtext = ",".join(str(x) for x in source_spws)

suffix='500k0.1mjy'
robust=0.5
niter=500000
field=['0','0']
specmode='mfs'
outframe='LSRK'
deconvolver='mtmfs'
nterms=2
gridder='standard'
weighting='briggs'
pbcor=True
pblimit=0.1
savemodel='modelcolumn'
interactive=False
start=''


imagename = 'sgr_b2.{0}.{1}.cont.r{2}.{3}.{4}'.format(fld, conf, robust, suffix, cal)

if not os.path.exists("{0}.image.tt0.pbcor.fits".format(imagename)):
    print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
    tclean(vis=mslist,
        imagename=imagename,
                   spw=[spwtext, spwtext],
                   field=field,
                   specmode=specmode,
                   start=start,
                   outframe=outframe,
                   threshold=threshold,
                   imsize=imsize,
                   cell=cell,
                   niter=niter,
                   deconvolver=deconvolver,
                   nterms=nterms,
                   gridder=gridder,
                   weighting=weighting,
                   robust=robust,
                   pbcor=pbcor,
                   mask=mask,
                   pblimit=pblimit,
                   savemodel=savemodel,
                   interactive=interactive)
    makefits(imagename)
else:
    print("This file already exists")

if step == 0:
    print('No calibration. Exiting')
    sys.exit()

gaintype='T'
refant='DV09'
combine='spw'
spw=spwcont
minsnr=3.0
minblperant=4

spwmap=[0,0,0,0]
calwt=False
flagbackup=True
applymode='calonly'

tbl1='pcal'+fld+conf+str(step)+'_1'
tbl2='pcal'+fld+conf+str(step)+'_2'

os.system('rm -r tbl1')
os.system('rm -r tbl2')


gaincal(vis=mslist[0],
        caltable=tbl1,
        gaintype=gaintype,
        refant=refant,
        calmode=calmode,
        combine=combine,
        spw=spw,
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant)

gaincal(vis=mslist[1],
        caltable=tbl2,
        gaintype=gaintype,
        refant=refant,
        calmode=calmode,
        combine=combine,
        spw=spw,
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant)


applycal(vis=mslist[0],
         spwmap=spwmap,
         gaintable=[tbl1],
         calwt=calwt,
         flagbackup=flagbackup,
         applymode=applymode)

applycal(vis=mslist[1],
         spwmap=spwmap,
         gaintable=[tbl2],
         calwt=calwt,
         flagbackup=flagbackup,
         applymode=applymode)

if step == 3:
    print('No splitting (last iteration). Exiting')
    sys.exit()

split(vis=mslist[0],
      outputvis=mslist[0]+'_cal'+str(step+1),
      datacolumn='corrected')
split(vis=mslist[1],
      outputvis=mslist[1]+'_cal'+str(step+1),
      datacolumn='corrected')
