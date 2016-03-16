from astropy.table import Table
from scipy.spatial import Voronoi
import matplotlib.path as mplPath
import astropy.wcs as wcs
from astropy.io import fits
import astropy.units as u
import numpy as np
from galaxies import Galaxy

m100 = Galaxy('M100')
tbl = Table.read('m100.co10.kkms_props_cprops_withSFR_subsets.fits')
hdr = fits.Header.fromtextfile('m100.header')
w = wcs.WCS(hdr)
ypix,xpix = np.meshgrid(np.linspace(0,799,800),np.linspace(0,799,800))
ra,dec = w.celestial.wcs_pix2world(xpix,ypix,0)
radius = m100.radius(ra=ra,dec=dec)

rclds = m100.radius(ra=tbl[~tbl['Outliers']]['XPOS'],dec=tbl[~tbl['Outliers']]['YPOS'])
maxrad = np.max(rclds.value)
xcloud,ycloud,zcloud = w.wcs_world2pix(tbl['XPOS'],tbl['YPOS'],tbl['VPOS'],0)
vor = Voronoi(zip(xcloud,ycloud))
label = np.zeros((800,800),dtype=np.str)
keys = ['Nuclear','Interarm','Arm']

for idx,region in enumerate(vor.regions):
    verts = vor.vertices[region]
    path = mplPath.Path(verts)
    inout = path.contains_points(zip(xcloud,ycloud))
    regions = path.contains_points(zip(xpix.ravel(),ypix.ravel()))
    thislabel = tbl[inout][keys]
    print(thislabel)
    if -1 not in region:
        for k in keys:
            if thislabel[k]:
                label[regions.reshape(800,800)]=k

for thislabel,thisname in zip(['A','I','N'],['arm','interarm','nuclear']):
    map = (label==thislabel)*(radius.value<9e3)
    map = map.astype(np.float)
    hdu = fits.PrimaryHDU(map,header=hdr)
    hdu.writeto('m100.{0}.fits'.format(thisname),clobber=True)
    pix_solidangle = np.abs(np.linalg.det(w.celestial.pixel_scale_matrix))*u.deg**2
    pix_area = pix_solidangle.to(u.sr)*m100.distance**2/np.cos(m100.inclination)
    region_area = pix_area * map.sum()
    print "Area for {0} region is {1} kpc^2".format(thisname,region_area.to(u.kpc**2*u.sr))
