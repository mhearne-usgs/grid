#!/usr/bin/env python

#stdlib imports
import os.path
import sys
from collections import OrderedDict
import warnings

#third party imports
import rasterio
import numpy as np
from grid2d import Grid2D
from dataset import DataSetException,DataSetWarning

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
                    data = np.squeeze(data)
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
                        data = np.squeeze(data)
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
                        section1 = np.squeeze(section1)
                        section2 = src.read_band(1,window2)
                        section2 = np.squeeze(section2)
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
                #Put NaN's back in where nodata value was
                nodata = src.get_nodatavals()[0]
                if (data==nodata).any():
                    data[data == nodata] = np.nan
        return (data,geodict)

    def _getHeader(self):
        hdr = {}
        if sys.byteorder == 'little':
            hdr['BYTEORDER'] = 'LSBFIRST'
        else:
            hdr['BYTEORDER'] = 'MSBFIRST'
        hdr['LAYOUT'] = 'BIL'
        hdr['NROWS'],hdr['NCOLS'] = self._data.shape
        hdr['NBANDS'] = 1
        if self._data.dtype == np.uint8:
            hdr['NBITS'] = 8
            hdr['PIXELTYPE'] = 'UNSIGNEDINT'
        elif self._data.dtype == np.int8:
            hdr['NBITS'] = 8
            hdr['PIXELTYPE'] = 'SIGNEDINT'
        elif self._data.dtype == np.uint16:
            hdr['NBITS'] = 16
            hdr['PIXELTYPE'] = 'UNSIGNEDINT'
        elif self._data.dtype == np.int16:
            hdr['NBITS'] = 16
            hdr['PIXELTYPE'] = 'SIGNEDINT'
        elif self._data.dtype == np.uint32:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'UNSIGNEDINT'
        elif self._data.dtype == np.int32:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'SIGNEDINT'
        elif self._data.dtype == np.float32:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'FLOAT'
        elif self._data.dtype == np.float64:
            hdr['NBITS'] = 32
            hdr['PIXELTYPE'] = 'FLOAT'
        else:
            raise DataSetException('Data type "%s" not supported.' % str(self._data.dtype))
        hdr['BANDROWBYTES'] = hdr['NCOLS']*(hdr['NBITS']/8)
        hdr['TOTALROWBYTES'] = hdr['NCOLS']*(hdr['NBITS']/8)
        hdr['ULXMAP'] = self._geodict['xmin']
        hdr['ULYMAP'] = self._geodict['ymax']
        hdr['XDIM'] = self._geodict['xdim']
        hdr['YDIM'] = self._geodict['ydim']
        #try to have a nice readable NODATA value in the header file
        zmin = np.nanmin(self._data)
        zmax = np.nanmax(self._data)
        if self._data.dtype in [np.int8,np.int16,np.int32]:
            nodata = np.array([-1*int('9'*i) for i in range(3,20)])
            if zmin > nodata[-1]:
                NODATA = nodata[np.where(nodata < zmin)[0][0]]
            else: #otherwise just pick an arbitrary value smaller than our smallest
                NODATA = zmin - 1
        else:
            nodata = np.array([int('9'*i) for i in range(3,20)])
            if zmin < nodata[-1]:
                NODATA = nodata[np.where(nodata > zmin)[0][0]]
            else: #otherwise just pick an arbitrary value smaller than our smallest
                NODATA = zmax + 1
        hdr['NODATA'] = NODATA
        keys = ['BYTEORDER','LAYOUT','NROWS','NCOLS','NBANDS','NBITS','BANDROWBYTES','TOTALROWBYTES','PIXELTYPE',
                'ULXMAP','ULYMAP','XDIM','YDIM','NODATA']
        hdr2 = OrderedDict()
        for key in keys:
            hdr2[key] = hdr[key]
        return hdr2
    
    def save(self,filename,format='EHdr'):
        #http://webhelp.esri.com/arcgisdesktop/9.3/index.cfm?TopicName=BIL,_BIP,_and_BSQ_raster_files
        #http://resources.esri.com/help/9.3/arcgisdesktop/com/gp_toolref/conversion_tools/float_to_raster_conversion_.htm
        supported = ['EHdr']
        if format not in supported:
            raise DataSetException('Only "%s" file formats supported for saving' % str(supported))
        hdr = self._getHeader()
        data = self._data #create a reference to the data - this may be overridden by a downcasted version for doubles
        if self._data.dtype == np.float32:
            data = self._data.astype(np.float32) #so we can find/reset nan values without screwing up original data
            data[np.isnan(data)] = hdr['NODATA']
        elif self._data.dtype == np.float64:
            data = self._data.astype(np.float32)
            data[np.isnan(data)] = hdr['NODATA']
            warnings.warn(DataSetWarning('Down-casting double precision floating point to single precision'))

        data.tofile(filename)
        #write out the header file
        basefile,ext = os.path.splitext(filename)
        hdrfile = basefile+'.hdr'
        f = open(hdrfile,'wt')
        for key,value in hdr.iteritems():
            value = hdr[key]
            f.write('%s  %s\n' % (key,str(value)))
        f.close()
            
    
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

def _test():
    try:
        for dtype in [np.uint8,np.uint16,np.uint32,np.int8,np.int16,np.int32,np.float32,np.float64]:
            print 'Testing saving/loading of data with type %s...' % str(dtype)
            data = np.arange(0,16).reshape(4,4).astype(dtype)
            if dtype in [np.float32,np.float64]:
                data[1,1] = np.nan
            geodict = {'xmin':0.5,'xmax':3.5,'ymin':0.5,'ymax':3.5,'xdim':1.0,'ydim':1.0,'nrows':4,'ncols':4}
            gdalgrid = GDALGrid(data,geodict)
            gdalgrid.save('test.bil')
            gdalgrid2 = GDALGrid.load('test.bil')
            np.testing.assert_almost_equal(gdalgrid2.getData(),gdalgrid.getData())
            print 'Passed saving/loading of data with type %s...' % str(dtype)

    except Exception,obj:
        print 'Failed tests with message: "%s"' % str(obj)
    os.remove('test.bil')
        
        
if __name__ == '__main__':
    _test()
        

