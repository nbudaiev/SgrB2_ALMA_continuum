from astrodendro import Dendrogram, pp_catalog
from astropy.io import fits
from astropy import units as u
import matplotlib.pyplot as plt
#import regions
import numpy as np
from astropy import coordinates
from astropy import wcs
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
from astropy.table import QTable, Table, Column
import matplotlib.cm as cm
from pyregion.mpl_helper import properties_func_default
from astropy.visualization import simple_norm
from astropy.stats import mad_std
#from regions import DS9Parser #depreciated
import re
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
import warnings
from astropy import visualization
from astropy.io import ascii
# from astropy.table import Table
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp

from gaussfit_catalog import gaussfit_catalog
import pandas as pd
from astropy.table import vstack

import numpy as np
import astropy.constants as c
import astropy.units as u

from radio_beam import Beam

from dendrocat import RadioSource

from astropy.stats import sigma_clip
from regions import Regions
from regions import PixCoord, CirclePixelRegion, CircleAnnulusPixelRegion
from astropy.nddata.utils import Cutout2D

import os

os.chdir('/orange/adamginsburg/sgrb2/NB/the_end')


def read_in(path):
    """
    Simplify reading in the files.
    """
    fh=fits.open(path)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wcs = WCS(fh[0].header)
        hdr = fh[0].header
    datas=fh[0].data #
    return datas, wcs, hdr

def create_regions(datas, num, fldconf):
    size = np.array(datas.shape)
    im_center = size/2
    #center = PixCoord(im_center[0],im_center[1])
    center = PixCoord(im_center[1], im_center[1]) # Using [1] for both values to deal with Band3 cut data.
    middle_row = datas[int(im_center[1])]
    # DANGER!!!!!!!!!
    if fldconf == 'NB3':
        center = PixCoord(im_center[1], 1372.0)
        middle_row = datas[int(1372.0)]
    diameter = len(middle_row[~np.isnan(middle_row)])
    region_rad = (((diameter/2)/int(num)))
    
    region_list = []
    region_list += [CirclePixelRegion(center,region_rad)]
    inner_radius = region_rad
    outer_radius = inner_radius+region_rad
    for i in range(int(num)-1):
        region_list += [CircleAnnulusPixelRegion(center,inner_radius,outer_radius)]
        inner_radius = outer_radius
        outer_radius = inner_radius+region_rad
    #regions.write_ds9(region_list, 'regions_15rings_test.reg')
    return region_list


def set_up_regions(path,wX):
    """
    # I am pretty sure "regionsX" part can be cut from this function
    Modifies a CASA-created DS9 region file to a version that works with Astrodendro package. (and also .contains?)
    Adds 'source_X' text to each region, where X is a source number based on its RA.
    
    Outputs:
    includeX: original DS9 regions
    regionsX: modified version of DS9 regions that works with Astrodendro
    pixelX: same as regionsX, but with pixel coordinates.
    """
    includeX = Regions.read(path)
    strX = Regions.serialize(includeX, format='ds9')
    start = []
    for m in re.finditer('source=1', strX):
        start += [m.end()]
    for i in range(len(start)):
        loc = start[-1-i]
        x = str(i+1)
        strM = strX[:loc] + ' text={source_' + x +'}' + strX[loc:]
    parser = Regions.parse(strX, format='ds9')
    regionsX = parser.serialize(format='ds9')
    #regions.write_ds9(regionsX, 'regions_15_rings_test.reg')
    pixelX = [region.to_pixel(wX) for region in includeX]
    return includeX, regionsX, pixelX



def get_noise_regions(data, RMS, create_regions_out):
    all_regs = create_regions_out
    masks = [reg.to_mask() for reg in all_regs]
    masked_1d = [m.get_values(data) for m in masks]
    rmss = [sigma_clip(m,sigma=7).std() for m in masked_1d]
    
    rmss_adj=rmss[:]
    lowest = np.min(rmss)
    for i in range(len(rmss)):
        if rmss[i] > lowest:
            rmss_adj[i] = lowest
        else:
            break
            
    rmss_adj = np.array(rmss_adj)*RMS
    min_values = []
    [min_values.append(x) for x in rmss_adj if x not in min_values]
    min_values = rmss_adj[:]
    
    return min_values


def run_dendros(datas, min_values, RMS, min_delta_factor, min_npix):
    dendros = []
    min_delta = (min_values[0]/RMS)*5
    dndr = Dendrogram.compute(datas, min_value = min_values[0], min_delta = min_delta, min_npix = min_npix)
    dendros.append(dndr)

    temp_min_values = min_values[1:]
    for i in range(len(temp_min_values)):
        if temp_min_values[i] == temp_min_values[i-1]:
            dendros.append(dndr)
        else:
            min_delta = (temp_min_values[i]/RMS)*min_delta_factor
            dndr = Dendrogram.compute(datas, min_value = temp_min_values[i], min_delta = min_delta, min_npix = min_npix)
            dendros.append(dndr)
    return dendros   


def set_up_catalog(d, wX, hdrX, includeX, extendedX, central=True):
    """
    Creates a dendrogram catalog.
    Removes sources that are not inside includeX regions.
    Removes sources inside extendedX regions.
    Removes sources that are above level 5 in the dendrogram (removes some unwated 
    detections from extended structures).
    """
    
    my_beam = Beam.from_fits_header(hdrX)  
    scale = wX.proj_plane_pixel_scales()[0]
    
    cat = pp_catalog(d.leaves,metadata={'data_unit':u.Jy / u.beam,'beam_major':my_beam.major,'beam_minor':my_beam.minor,'beam_pa':my_beam.pa, 'wcs': wX, 'spatial_scale':scale},verbose=False)
    
    peaks = []
    peak_x_px = []
    peak_y_px = []
    peak_x_wcs = []
    peak_y_wcs = []
    npix = []
    for i in range(len(d.leaves)):
        (y, x), peak = d.leaves[i].get_peak()
        peak_x_px += [x] # x 
        peak_y_px += [y] # y 
        peaks += [peak]
        npix += [d.leaves[i].get_npix()]
        
        coordsX_peak_wcs = PixCoord(x,y).to_sky(wX)
        peak_x_wcs += [coordsX_peak_wcs.ra.deg]
        peak_y_wcs += [coordsX_peak_wcs.dec.deg]
        
        
    cat.add_column(peaks,name='peak',index=4)
    cat.add_column(peak_x_px, name = 'x_peak_px', index=5)
    cat.add_column(peak_y_px, name = 'y_peak_px', index=6)
    cat.add_column(peak_x_wcs, name = 'x_peak', index=7)
    cat.add_column(peak_y_wcs, name = 'y_peak', index=8)
    cat.add_column(npix,name='npix',index=9)

    print('Total detections: '+str(len(d.leaves)))
    purgeX = []
    for i in range(len(d)):
        if d[i].level > 5: 
            purgeX += [d[i].idx]
    coordsX_wcs = coordinates.SkyCoord(cat['x_cen'], cat['y_cen'], unit=(u.deg,u.deg), frame=wX.celestial.wcs.radesys.lower())
    containsX = np.zeros(len(cat),dtype=bool)

    coordsX_pix = PixCoord.from_sky(coordsX_wcs, wX)
    containsX[includeX.contains(coordsX_pix)] = 1

    print('Total detections inside input regions: '+str(np.sum(containsX)))
    cat_cut_temp = cat[containsX]
    cat_coords_temp = coordinates.SkyCoord(cat_cut_temp['x_cen'], cat_cut_temp['y_cen'], unit=(u.deg,u.deg), frame=wX.celestial.wcs.radesys.lower())
    
    ####
    containsX1 = np.zeros(len(cat_cut_temp),dtype=bool)
    for X in extendedX:
        containsX1[X.contains(cat_coords_temp, wX)] = 1

    print('Total detections within extended structures (removed): '+str(np.sum(containsX1)))
    cat_cut = cat_cut_temp[~containsX1]
    
    ####
    
    if central:
        print('Central region. No trimming performed')
        cat_coords = coordinates.SkyCoord(cat_cut['x_cen'], cat_cut['y_cen'], unit=(u.deg,u.deg), frame=wX.celestial.wcs.radesys.lower())
        return cat_coords, cat_cut, cat_cut['_idx']
    else:
        cat_final= cat_cut[:0].copy()
        for i in range(len(cat_cut)):    
            if cat_cut['_idx'][i] not in purgeX:
                cat_final.add_row(cat_cut[i])
        print('Total detections after triming above level 5: '+str(len(cat_final)))
        cat_coords = coordinates.SkyCoord(cat_final['x_cen'], cat_final['y_cen'], unit=(u.deg,u.deg), frame=wX.celestial.wcs.radesys.lower())
        return cat_coords, cat_final, cat_final['_idx']
    
    
# This creates list instead of astropy tables
def catalog_from_rings(dendros, wcs_, hdr, all_regs, extended):
    cat_coords = []
    cat_final = []
    IDs = []
    concentric_ID = []
    for i in range(len(dendros)):
        if i == 0:
            cat_coords, cat_final, IDs = set_up_catalog(dendros[i],wcs_,hdr, all_regs[i], extended)
            #wcs_coords_, wcs_ = set_up_catalog(dendros[i],wcs,hdr, all_regs[i], extended)
            #return wcs_coords_, wcs_
        ########################
            #cat_coords, cat_final, IDs = set_up_catalog(dendros[i],wcs_,hdr, all_regs[i])
            concentric_ID = [i] * len(cat_final)
        else:
            x,y,z = set_up_catalog(dendros[i],wcs_,hdr,all_regs[i], extended)
            CID = [i] * len(y)
            try:
                cat_coords = vstack([cat_coords,x])
                cat_final = vstack([cat_final,y])
                IDs = vstack([IDs])
                concentric_ID = vstack([concentric_ID,CID])
            except:
                pass
    return cat_coords, cat_final, IDs, concentric_ID

def add_pix_coords(cat, wcs):
    xwcs = np.array(cat['x_cen'])
    ywcs = np.array(cat['y_cen'])
    
    c_wcs = SkyCoord(xwcs, ywcs, unit="deg", frame="icrs")
    c_pix = wcs.world_to_pixel(c_wcs)
    
    cat.add_column(c_pix[0],name='x_cen_px')
    cat.add_column(c_pix[1],name='y_cen_px')

def get_radius(datas, cat, wcs_, fldconf):
    # maybe add physical units
    xc, yc = datas.shape[1]/2, datas.shape[1]/2 # careful with datas variable here
    if fldconf == 'NB3':
        xc, yc = datas.shape[1]/2, 1372.0
    x = cat['x_cen_px']
    y = cat['y_cen_px']
    radius_px = np.sqrt((x - xc)**2 + (y - yc)**2)

    # This is really messy and probably doesn't work correctly
    temp_sky_coords = SkyCoord(cat['x_cen'], cat['y_cen'], unit=(u.deg,u.deg), distance=8.5*u.kpc, frame='icrs')
    c1= temp_sky_coords
    xc_wcs =  wcs_.pixel_to_world(xc, yc).ra
    yc_wcs =  wcs_.pixel_to_world(xc, yc).dec
    center_sky_coords = SkyCoord(xc_wcs, yc_wcs, distance = 8.5*u.kpc, frame = 'icrs')
    c2 = center_sky_coords
    sep = c1.separation_3d(c2).to(u.AU)
    
    return radius_px, sep

def add_columns(datas, cat, wcs_, concentric_ID, fldconf):
    cat.add_column(concentric_ID['col0'],name='c_ID', index = 1)
    add_pix_coords(cat, wcs_)
    radius_px, radius_AU = get_radius(datas, cat, wcs_, fldconf)
    cat.add_column(radius_px, name='radial_d_px')
    cat.add_column(radius_AU, name='radial_d_AU')
    if not '_name' in cat.colnames:
        fwhm_maj = cat['major_sigma']*np.sqrt(8*np.log(2))
        fwhm_min = cat['minor_sigma']*np.sqrt(8*np.log(2))
        cat.add_column(fwhm_maj,index=11,name='major_fwhm')
        cat.add_column(fwhm_min,index=12,name='minor_fwhm')
        #cat.add_column(cat['_idx'],index=1,name='_name')
        cat.add_column(0.0,name='rejected')
        numbered = list(range(1, len(cat)+1))
        cat.add_column(numbered, name = '_name', index = 2)

def regions_gaussfit(cat_final, fldconf):
    """
    Creates a DS9 region file from the catalog. The centers and sizes of each region
    are based on ['x_cen'], ['y_cen'], and ['radius'] columns in the catalog.
    """
    centerX = SkyCoord(cat_final['x_cen'][::-1], cat_final['y_cen'][::-1],unit='deg')
    # There are some problems with displaying these regions in DS9.
    # Use 1st option to display the regions in DS9.
    # Use 2nd option to have proper radii for later use.
    #radiusX = Angle(cat_final['radius'][::-1].value*.0000055, unit='deg')
    #radiusX = Angle(cat_final['radius'][::-1].value*0.02, unit='arcsec')
    
    # fixed stuff?
    radiusX = Angle(cat_final['radius'][::-1].value, unit='deg') #### IMPORTANT!!! ONLY RADIUS IS BACKWARDS, X_CEN AND Y_CEN ARE NOT???
    cat_regionsX = []
    for i in range(len(centerX)):
        cat_regionsX += [CircleSkyRegion(centerX[i], radiusX[i])]
    #str_cat_regionsX = regions.ds9_objects_to_string(cat_regionsX, coordsys='icrs') # regions 0.5
    str_cat_regionsX = Regions(cat_regionsX).serialize(format='ds9')
    start = []
    for m in re.finditer('\)', str_cat_regionsX):
        start += [m.end()]
    str1 = ' # select=1 highlite=1 fixed=0 edit=1 move=1 delete=1 source=1 color=#2EE6D6 dashlist=8 3 width=2 dash=0 font="helvetica 10 normal roman" text={source_'
    str2 = '}'
    for i in range(len(start)):
        loc = start[-1-i]
        x = str(i+1)
        str_cat_regionsX = str_cat_regionsX[:loc] + str1 + x + str2 + str_cat_regionsX[loc:]
    #parser = DS9Parser(str_cat_regionsX) #regions 0.5
    parser = Regions.parse(str_cat_regionsX, format='ds9')
    #str_cat_regionsX = parser.shapes.to_regions()
    #regions.write_ds9(str_cat_regionsX, ''.reg")
    #return str_cat_regionsX
    return parser
    
    
def catalog_for_configuration(fldconf, RMS = 5.0, num_regions = 7, min_delta_factor = 0.5, min_npix = 3):
    if fldconf == 'NB3':
        #path = 'sgr_b2.N.B3.cont.r0.5.1m0.075mJy.cal3.image.tt0.pbcor.fits'
        path = 'sgr_b2.N.B3.cont.r0.5.1m0.075mJy.cal2.image.tt0.pbcor.fits'
        position = (3072, 3922)
        size = (4444, 6144)
    if fldconf == 'MB3':
        path = 'sgr_b2.M.B3.cont.r0.5.1m0.125mJy.cal3.image.tt0.pbcor.fits'
        position = (3072, 2135)
        size = (4270, 6144)
    if fldconf == 'NB6':
        path = 'sgr_b2.N.B6.cont.r0.5.1m1.5mJy.cal4.image.tt0.pbcor.fits'
    if fldconf == 'MB6':
        path = 'sgr_b2.M.B6.cont.r0.5.1m0.68mJy.cal3.image.tt0.pbcor.fits'

    datas, wcs_, hdr = read_in(path)

    if fldconf == 'NB3' or fldconf == 'MB3':
        cutout = Cutout2D(datas, position, size, wcs_)
        datas = cutout.data
        wcs_ = cutout.wcs
    
    regions_rings = create_regions(datas, num_regions, fldconf) # number of rings
    min_values = get_noise_regions(datas, RMS, regions_rings)
    
    extended, regions_, pixel_ = set_up_regions(fldconf+'_extended_structures.reg',wcs_)
    
    print('Running dendrograms')
    dendros = run_dendros(datas, min_values, RMS, min_delta_factor, min_npix)
    print('Finished running dendrograms')
    
    print('Getting catalog')
    cat_coords, cat_final, IDs, concentric_ID = catalog_from_rings(dendros, wcs_, hdr, regions_rings, extended)  # is IDs really needed?
    print('Done cataloging')
    
    #newly added:
    if np.min(cat_final['minor_sigma']) == 0:
        print('Spotted "1D" sources')
        print('Length before'+ str(len(cat_final)))
        concentric_ID = concentric_ID[cat_final['minor_sigma']>0]
        cat_final = cat_final[cat_final['minor_sigma']>0]

        print('Length after'+ str(len(cat_final)))

    add_columns(datas, cat_final, wcs_, concentric_ID, fldconf)
    
    # This is to fix the pixel coordinates being recorded based on the cut image
    if fldconf == 'NB3':
        cat_final['y_cen_px'] = cat_final['y_cen_px']+1700
        cat_final['y_peak_px'] = cat_final['y_peak_px']+1700
    

    catalogname='catalog_'+fldconf+'_'+str(int(RMS))+'RMS_'+str(num_regions)+'rings_'+str(min_delta_factor)+'mindelta'+str(min_npix)+'npix'+'.csv'
    ascii.write(cat_final, catalogname, format='csv', overwrite=True)

    #str_cat_regionsX = regions_gaussfit(cat_final, fldconf)
    #regions_name = 'regions_'+fldconf+'_'+str(int(RMS))+'RMS_'+str(num_regions)+'rings_'+str(min_delta_factor)+'mindelta'+str(min_npix)+'npix'+'.reg'
    #str_cat_regionsX.write(regions_name, overwrite=True, format='ds9')
    
import sys
    
catalog_for_configuration(fldconf = sys.argv[1], RMS = float(sys.argv[2]), num_regions = float(sys.argv[3]), min_delta_factor = float(sys.argv[4]), min_npix = 3)