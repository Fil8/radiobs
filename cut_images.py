import string
from astropy.io import fits, ascii
from astropy import wcs
import numpy as np


def ra2deg(ra_degrees):
    
    ra = string.split(ra_degrees, ':')

    hh = float(ra[0])*15
    mm = (float(ra[1])/60)*15
    ss = (float(ra[2])/3600)*15
    
    return hh+mm+ss


# DMS -> degrees
def dec2deg(dec_degrees):
    dec = string.split(dec_degrees, ':')
    
    hh = abs(float(dec[0]))
    mm = float(dec[1])/60
    ss = float(dec[2])/3600
    
    return hh+mm+ss

xmin_ra='03:25:11.242'
ra_min=ra2deg(xmin_ra)
ymin_dec='-37:41:57.45'
dec_min=-dec2deg(ymin_dec)

xmax_ra='03:20:14.275'
ra_max=ra2deg(xmax_ra)

ymax_dec='-36:42:50.0'
dec_max=-dec2deg(ymax_dec)

t=ascii.read('maps_filenames.def')

maps_filenames = t['MAP']

for i in xrange(0,len(maps_filenames)):
	
	print maps_filenames[i]
	f=fits.open(maps_filenames[i])
	d=f[0].data
	h=f[0].header
	print d.shape
	w = wcs.WCS(h)    

	xmin,ymin=w.wcs_world2pix(ra_min,dec_min,0)
	print xmin,ymin
	xmin=int(np.round(xmin,0))
	ymin=int(np.round(ymin,0))
	
	xmax,ymax=w.wcs_world2pix(ra_max,dec_max,0)
	xmax=int(np.round(xmax,0))
	ymax=int(np.round(ymax,0))
	naxis1=xmax-xmin
	naxis2=ymax-ymin
	print naxis1,naxis2
	
	h['NAXIS1']=naxis1
	h['NAXIS2']=naxis2
	h['CRPIX1']=naxis1/2
	h['CRPIX2']=naxis2/2

    
	aaa = string.split(maps_filenames[i], '.fits')
	output=aaa[0]+'_masked.fits'
	fits.writeto(output,d[ymin:ymax,xmin:xmax],h,overwrite=True)
