import datetime
import os
import sys

def makefits(myimagebase, cleanup=True):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
   # exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
   # exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
   # exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
   # exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
   # exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)

    if cleanup:
        for ttsuffix in ('.tt0', '.tt1', 'tt2'):
            for suffix in ('pb{tt}', 'weight', 'sumwt{tt}', 'psf{tt}',
                           'model{tt}', 'mask', 'image{tt}', 'residual{tt}',
                           'alpha', ):
                os.system('rm -rf {0}.{1}'.format(myimagebase, suffix).format(tt=ttsuffix))

# Scripit to reduce all cont data.

#====================================
fld = 'N' # change to N/M as needed
conf = 'B3' # change to B3/B6 as needed
step = 4 # change to 0/1/2/3. 
# 0 - initial clean
# 1 - p,inf cal + clean
# 2 - p,int cal + clean
# 3 - ap,15s cal + clean
#====================================
print('Field: '+fld)
print('Configuration: '+conf)
print('Step: '+str(step))


fldconf = fld+conf



if conf == 'B3':
    imsize = [6144, 6144]
    cell = ['0.018arcsec']
    if step == 0 or step == 1:
        mask = fld+conf+'_mask.mask'
        print('Using tight mask')
    else:
        mask = 'SgrB2'+fldconf+'_regions.crtf'
        print('Using broad mask')

if conf == 'B6':
    imsize = [4500,4500]
    cell =  ['0.01arcsec']
    mask = 'SgrB2'+fldconf+'_regions_new.crtf'

print('Using mask: '+ mask)

if step == 0 or step == 1:
    if fldconf == 'NB3':
        threshold = '0.15mJy'
    if fldconf == 'MB3':
        threshold = '0.25mJy'
    if fldconf == 'MB6':
        threshold = '1.36mJy'

if step == 2 or step == 3 or step == 4: # can be simplified by saying "threshold = threshold * 2"
    if fldconf == 'NB3':
        threshold = '0.075mJy'
    if fldconf == 'MB3':
        threshold = '0.125mJy'
    if fldconf == 'MB6':
        threshold = '0.68mJy'
if fldconf == 'NB6':
    threshold = '7mJy'
    if step == 3 or step == 4:
        threshold = '1.5mJy' # This is 3 RMS. Hope this works. 2.4mJy worked before.

print('Using threshold value of: '+threshold)

ms1=fld+'_'+conf+'_cont1.ms'
ms2=fld+'_'+conf+'_cont2.ms'

mslist_original = [ms1,ms2]


cal = 'cal'+str(step)

if step == 1:
    solint='inf'
    calmode='p'
elif step == 2:
    solint='int'
    calmode='p'
elif step == 3:
    solint='inf'
    calmode='ap'
elif step == 4: # TEST THIS, WHY NOT?
    solint='15s'
    calmode='ap'

#mslist = [s + '_cal' + str(step-1) for s in mslist_original]
mslist=[ms1+'_cal'+str(step-1),ms2+'_cal'+str(step-1)]

new_mslist=[ms1+'_cal'+str(step),ms2+'_cal'+str(step)]

source_spws = (0,1,2,3)

spwtext = ",".join(str(x) for x in source_spws)

suffix='1m'+threshold
robust=0.5
niter=1000000
field=['0','0']
specmode='mfs'
outframe='LSRK'
deconvolver='mtmfs'
#deconvolver='hogbom'
nterms=2
gridder='standard'
weighting='briggs'
pbcor=True
pblimit=0.1
savemodel='modelcolumn'
interactive=False


imagename = 'sgr_b2.{0}.{1}.cont.r{2}.{3}.{4}'.format(fld, conf, robust, suffix, cal)


if step == 0:
    if not os.path.exists("{0}.image.tt0.pbcor.fits".format(imagename)):
        print("The folders are: " + str(mslist_original))
        print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
        tclean(vis=mslist_original,
               imagename=imagename,
               spw=[spwtext, spwtext],
               field=field,
               specmode=specmode,
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
               interactive=interactive,)
        makefits(imagename)
    
    os.system('rm -r '+new_mslist[0])
    os.system('rm -r '+new_mslist[1])

    #split(vis=mslist_original[0], # DO NOT USE SPLIT HERE!
    #      outputvis=new_mslist[0],
    #      datacolumn='data')
    #split(vis=mslist_original[1],
    #      outputvis=new_mslist[1],
    #      datacolumn='data')
    os.system('cp -r '+mslist_original[0]+" "+new_mslist[0])
    os.system('cp -r '+mslist_original[1]+" "+new_mslist[1])

    print('No calibration. Exiting')
    sys.exit()

gaintype='T'
refant='DV09,DV06'
combine='spw'
spwcont='0,1,2,3'
minsnr=3.0
minblperant=4

spwmap=[0,0,0,0]
calwt=False
flagbackup=True
applymode='calonly'

tbl1='pcal'+fld+conf+str(step)+'_1'
tbl2='pcal'+fld+conf+str(step)+'_2'

os.system('rm -r '+tbl1)
os.system('rm -r '+tbl2)

print("Calibrating folders: " + str(mslist))

gaincal(vis=mslist[0],
        caltable=tbl1,
        gaintype=gaintype,
        refant=refant,
        calmode=calmode,
        combine=combine,
        #spw=spwcont,
        solint=solint,
        minsnr=minsnr,
        minblperant=minblperant)

gaincal(vis=mslist[1],
        caltable=tbl2,
        gaintype=gaintype,
        refant=refant,
        calmode=calmode,
        combine=combine,
        #spw=spwcont,
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

os.system('rm -r '+new_mslist[0])
os.system('rm -r '+new_mslist[1])

split(vis=mslist[0],
      outputvis=new_mslist[0],
      datacolumn='corrected')
split(vis=mslist[1],
      outputvis=new_mslist[1],
      datacolumn='corrected')


if not os.path.exists("{0}.image.tt0.pbcor.fits".format(imagename)):
    print("Imaging folders: " + str(new_mslist))
    print("Imaging {0} at {1}".format(imagename, datetime.datetime.now()))
    tclean(vis=new_mslist,
        imagename=imagename,
                   spw=[spwtext, spwtext],
                   field=field,
                   specmode=specmode,
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
                   interactive=interactive,
                   antenna=antenna)
    makefits(imagename)
else:
    print("This file already exists")
