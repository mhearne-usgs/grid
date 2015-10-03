#!/usr/bin/env python

import rasterio
import numpy as np
from grid2d import Grid2D
from dataset import DataSetException

class GDALGrid(Grid2D):
    def __init__(self,data,geodict):
        """Construct a GMTGrid object.
        :param data:
           2D numpy data array (must match geodict spec)
        :param geodict:
           Dictionary specifying the spatial extent,resolution and shape of the data.
        :returns:
           A GMTGrid object.
        :raises DataSetException:
          When data and geodict dimensions do not match. 
        """
        m,n = data.shape
        if m != geodict['nrows'] or n != geodict['ncols']:
            raise DataSetException('Input geodict does not match shape of input data.')
        self._data = data
        self._geodict = geodict

    @classmethod
    def getFileGeoDict(cls,filename):
        geodict = {}
        with rasterio.drivers():
            with rasterio.open(filename) as src:
                aff = src.affine
                geodict['xmin'] = aff.xoff
                geodict['ymax'] = aff.yoff
                
                geodict['xdim'] = aff.a
                geodict['ydim'] = aff.e
                shp = src.shape
                if len(shp) > 2:
                    raise DataSetException('Cannot support grids with more than one band')
                geodict['nrows'] = src.height
                geodict['ncols'] = src.width
                geodict['xmax'] = geodict['xmin'] + geodict['ncols']*geodict['xdim']
                geodict['ymin'] = geodict['ymax'] - geodict['nrows']*geodict['ydim']

        return geodict

    @classmethod
    def readGDAL(cls,filename,bounds=None):
        geodict = cls.getFileGeoDict(filename)
        data = None
        with rasterio.drivers():
            with rasterio.open(filename) as src:
                if bounds is None:
                    if len(src.shape) > 2:
                        raise DataSetException('Cannot support grids with more than one band')
                    data = src.read()
                else:
                    xmin,xmax,ymin,ymax = bounds
                    gxmin,gxmax,gymin,gymax = (geodict['xmin'],geodict['xmax'],geodict['ymin'],geodict['ymax'])
                    xdim,ydim = (geodict['xdim'],geodict['ydim'])
                    if xmin < xmax:
                        ixmin = int((xmin - gxmin)/xdim)
                        ixmax = int((xmax - gxmin)/xdim)
                        iymin = int((ymin - gymin)/ydim)
                        iymax = int((ymax - gymin)/ydim)
                        window = ((iymin,iymax),(ixmin,ixmax))
                        data = src.read_band(1,window)
                        geodict['xmin'] = xmin
                        geodict['xmax'] = xmax
                        geodict['ymin'] = ymin
                        geodict['ymax'] = ymax
                        m,n = data.shape
                        geodict['nrows'] = m
                        geodict['ncols'] = m
                    else:
                        #cut user's request into two regions - one from the minimum to the
                        #meridian, then another from the meridian to the maximum.
                        (region1,region2) = super(GDALGrid,cls)._createSections((xmin,xmax,ymin,ymax),geodict)
                        (iulx1,iuly1,ilrx1,ilry1) = region1
                        (iulx2,iuly2,ilrx2,ilry2) = region2
                        outcols1 = long(ilrx1-iulx1+1)
                        outcols2 = long(ilrx2-iulx2+1)
                        outcols = long(outcols1+outcols2)
                        outrows = long(ilry1-iuly1+1)
                        window1 = (iuly1,ilry1,iulx1,ilrx1)
                        window2 = (iuly2,ilry2,iulx2,ilrx2)
                        section1 = src.read_band(1,window1)
                        section2 = src.read_band(1,window2)
                        data = np.concatenate((section1,section2),axis=1)
                        xmin = (ulx + iulx1*xdim)
                        ymax = uly - iuly1*ydim
                        xmax = ulx + ilrx2*xdim
                        ymin = bymax - outrows*ydim
                        geodict['xmin'] = xmin
                        geodict['xmax'] = xmax
                        geodict['ymin'] = ymin
                        geodict['ymax'] = ymax
                        geodict['nrows'],geodict['ncols'] = data.shape
        return (data,geodict)

    def save(self,filename,format='EHdr'):
        env = rasterio._drivers.GDALEnv()
        env.start()
        fdict = env.drivers()
        env.stop()
        if format not in fdict.keys():
            raise DataSetException('Format "%s" not supported.  Call GDALGrid().getFormats() for supported list.' % format)
        with rasterio.drivers():
            with rasterio.open(filename) as src:
                pass
    
    @classmethod
    def load(cls,filename,samplegeodict=None,resample=False,method='linear',doPadding=False,padValue=np.nan):
        if samplegeodict is None:
            data,geodict = cls.readGDAL(filename)
        else:
            bounds = (samplegeodict['xmin'],samplegeodict['xmax'],samplegeodict['ymin'],samplegeodict['ymax'])
            #if the user wants resampling, we can't just read the bounds they asked for, but instead
            #go outside those bounds.  if they asked for padding and the input bounds exceed the bounds
            #of the file, then we can pad.  If they *didn't* ask for padding and input exceeds, raise exception.
            if resample:
                PADFACTOR = 2 #how many cells will we buffer out for resampling?
                filegeodict = cls.getFileGeoDict(gmtfilename)
                xdim = filegeodict['xdim']
                ydim = filegeodict['ydim']
                fbounds = (filegeodict['xmin'],filegeodict['xmax'],filegeodict['ymin'],filegeodict['ymax'])
                isOutside = False
                #make a bounding box that is one row/col greater than what the user asked for
                rbounds = (bounds[0]-xdim*PADFACTOR,bounds[1]+xdim*PADFACTOR,bounds[2]-ydim*PADFACTOR,bounds[3]+ydim*PADFACTOR)
                #compare that bounding box to the file bounding box
                if fbounds[0] > rbounds[0] or fbounds[1] < rbounds[1] or fbounds[2] > rbounds[2] or fbounds[3] < rbounds[3]:
                    isOutside = True
                if isOutside:
                    if doPadding==False:
                        raise DataSetException('Cannot resample data given input bounds, unless doPadding is set to True.')
                    else:
                        samplebounds = bounds #we'll just read what's in the file and pad later?
                else:
                    samplebounds = rbounds
            #TODO - fill in with code!
            data = None
            geodict = None
        return cls(data,geodict)


def _test_load():
    pass

def createSample(M,N):
    data = np.arange(0,M*N).reshape(M,N)
    data = data.astype(np.int32) #arange gives int64 by default, not supported by netcdf3
    xvar = np.arange(0.5,0.5+N,1.0)
    yvar = np.arange(0.5,0.5+M,1.0)
    geodict = {'nrows':N,
               'ncols':N,
               'xmin':0.5,
               'xmax':xvar[-1],
               'ymin':0.5,
               'ymax':yvar[-1],
               'xdim':1.0,
               'ydim':1.0}
    gmtgrid = GDALGrid(data,geodict)
    return gmtgrid

def _test_save():
    try:
        gdalgrid = createSample(4,3)
        gdalgrid.save('test.grd',format='netCDF')
        
            
        

