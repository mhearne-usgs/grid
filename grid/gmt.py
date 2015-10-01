#!/usr/bin/env python

#stdlib imports
import struct
import os.path

#third party imports
import numpy as np
from scipy.io import netcdf
from grid2d import Grid2D
from gridbase import Grid
from dataset import DataSetException
import h5py

#TODO:
#1) Figure out how COARDS specifies pixel vs gridline registration (default seems to be gridline)
#2) Write a test for 180 meridian crossing grids

def __createSections(bounds,geodict):
    """Given a grid that goes from 0 to 180 degrees, figure out the two pixel regions that up both sides of the subset
    :param bounds:
       Tuple of (xmin,xmax,ymin,ymax)
    :param geodict:
       Geodict dictionary
    :returns:
      Two tuples of 4 elements each - (iulx,iuly,ilrx,ilry). The first tuple defines the pixel offsets for the left
      side of the subsetted region, and the second tuple defines the pixel offsets for the right side.
    """
    (bxmin,bxmax,bymin,bymax) = bounds
    ulx = geodict['xmin']
    uly = geodict['ymax']
    xdim = geodict['xdim']
    ydim = geodict['ydim']
    ncols = geodict['ncols']
    nrows = geodict['nrows']
    #section 1
    iulx1 = int(np.floor((bxmin - ulx)/xdim))
    iuly1 = int(np.ceil((uly - bymax)/ydim))
    ilrx1 = int(ncols-1)
    ilry1 = int(np.floor((uly - bymin)/ydim))
    #section 2
    iulx2 = 0
    iuly2 = int(np.ceil((uly - bymax)/ydim))
    ilrx2 = int(np.ceil((bxmax - ulx)/xdim))
    ilry2 = int(np.floor((uly - bymin)/ydim))

    region1 = (iulx1,iuly1,ilrx1,ilry1)
    region2 = (iulx2,iuly2,ilrx2,ilry2)
    return(region1,region2)

def createSampleGrid(N):
    """Used for internal testing - create an NxN grid with lower left corner at 0.5,0.5, xdim/ydim = 1.0
    :param N:
       Number of rows and columns in output grid
    :returns:
       GMTGrid object where data values are an NxN array of values from 0 to N-squared minus 1, and geodict
       lower left corner is at 0.5/0.5 and cell dimensions are 1.0.
    """
    data = np.arange(0,np.power(N,2)).reshape(N,N)
    data = data.astype(np.int32) #arange gives int64 by default, not supported by netcdf3
    xvar = np.arange(0.5,0.5+N,1.0)
    yvar = np.arange(0.5,0.5+N,1.0)
    geodict = {'nrows':N,
               'ncols':N,
               'xmin':0.5,
               'xmax':xvar[-1],
               'ymin':0.5,
               'ymax':yvar[-1],
               'xdim':1.0,
               'ydim':1.0}
    gmtgrid = GMTGrid(data,geodict)
    return gmtgrid

class GMTGrid(Grid2D):
    """
    A class that implements a Grid2D object around GMT NetCDF/HDF/native gridded data sets.
    """
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
    def getFileType(cls,grdfile):
        """Get the GMT sub-format (netcdf, hdf, or GMT binary).
        :param grdfile:
          File name of suspected GMT grid file.
        :returns:
          One of 'netcdf' (NetCDF version 3), 'hdf' (NetCDF version 4), 'binary' (so-called GMT native format), or 'unknown'.
        """
        f = open(grdfile,'rb')
        #check to see if it's HDF or CDF
        f.seek(1,0)
        hdfsig = ''.join(struct.unpack('ccc',f.read(3)))
        ftype = 'unknown'
        if hdfsig == 'HDF':
            ftype = 'hdf'
        else:
            f.seek(0,0)
            cdfsig = ''.join(struct.unpack('ccc',f.read(3)))
            if cdfsig == 'CDF':
                ftype = 'netcdf'
            else:
                f.seek(8,0)
                offset = struct.unpack('I',f.read(4))[0]
                if offset == 0 or offset == 1:
                    ftype = 'binary'
                    
        f.close()
        return ftype

    @classmethod
    def getFileGeoDict(cls,filename):
        """Get the spatial extent, resolution, and shape of grid inside NetCDF file.
        :param filename:
           File name of NetCDF file.
        :returns:
          GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
        :raises DataSetException:
          When the file is not detectable as one of the three flavors of GMT grids.
        """
        ftype = cls.getFileType(filename)
        if ftype == 'native':
            geodict = cls.getNativeHeader(filename)
        elif ftype == 'netcdf':
            geodict,xvar,yvar = cls.getNetCDFHeader(filename)
        elif ftype == 'hdf':
            geodict,xvar,yvar = cls.getHDFHeader(filename)
        else:
            raise DataSetException('Unknown file type for file "%s".' % filename)
        return geodict
            
    @classmethod
    def getNativeHeader(cls,fname):
        """Get the header information from a GMT native grid file.
        :param fname:
           File name of GMT native grid
        :returns:
          GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
        """
        f = open(grdfile,'rb')
        f.seek(0,0)
        geodict = {}
        geodict['ncols'] = struct.unpack('I',f.read(4))[0]
        geodict['nrows'] = struct.unpack('I',f.read(4))[0]
        offset = struct.unpack('I',f.read(4))[0]
        geodict['xmin'] = struct.unpack('d',f.read(8))[0]
        geodict['xmax'] = struct.unpack('d',f.read(8))[0]
        geodict['ymin'] = struct.unpack('d',f.read(8))[0]
        geodict['ymax'] = struct.unpack('d',f.read(8))[0]
        zmin = struct.unpack('d',f.read(8))[0]
        zmax = struct.unpack('d',f.read(8))[0]
        geodict['xdim'] = struct.unpack('d',f.read(8))[0]
        geodict['ydim'] = struct.unpack('d',f.read(8))[0]
        f.close()
        return geodict

            
    @classmethod
    def readGMTNative(cls,fname,bounds=None):
        """Read the data and geo-referencing information from a GMT native grid file, subsetting if requested.
        :param fname:
          File name of GMT native grid
        :param bounds:
           Tuple of (xmin,xmax,ymin,ymax)
        :returns:
          Tuple of data (2D numpy array of data, possibly subsetted from file) and geodict (see above).
        :raises NotImplementedError:
          For any bounds not None (we'll get to it eventually!)
        """
        f = open(grdfile,'rb')
        f.seek(0,0)
        geodict = {}
        geodict['ncols'] = struct.unpack('I',f.read(4))[0]
        geodict['nrows'] = struct.unpack('I',f.read(4))[0]
        offset = struct.unpack('I',f.read(4))[0]
        geodict['xmin'] = struct.unpack('d',f.read(8))[0]
        geodict['xmax'] = struct.unpack('d',f.read(8))[0]
        geodict['ymin'] = struct.unpack('d',f.read(8))[0]
        geodict['ymax'] = struct.unpack('d',f.read(8))[0]
        zmin = struct.unpack('d',f.read(8))[0]
        zmax = struct.unpack('d',f.read(8))[0]
        geodict['xdim'] = struct.unpack('d',f.read(8))[0]
        geodict['ydim'] = struct.unpack('d',f.read(8))[0]
        zscale = struct.unpack('d',f.read(8))[0]
        zoffset = struct.unpack('d',f.read(8))[0]
        xunits = f.read(80).strip()
        yunits = f.read(80).strip()
        zunits = f.read(80).strip()
        title = f.read(80).strip()
        command = f.read(320).strip()
        remark = f.read(160).strip()
        #nota bene - the extent specified in a GMT grid is for the edges of the
        #grid, regardless of whether you've specified grid or pixel
        #registration.
        #TODO  - test above assertion!
        geodict['xmin'] = geodict['xmin'] + geodict['xdim']/2.0
        geodict['xmax'] = geodict['xmax'] - geodict['xdim']/2.0
        geodict['ymin'] = geodict['ymin'] + geodict['ydim']/2.0
        geodict['ymax'] = geodict['ymax'] - geodict['ydim']/2.0
        if bounds is None:
            sfmt = '%i%s' % (geodict['ncols']*geodict['nrows'],fmt)
            dwidths = {'i':2,'l':4,'f':4,'d':8}
            dwidth = dwidths[fmt]
            dbytes = f.read(geodict['ncols']*geodict['nrows']*dwidth)
            data = struct.unpack(sfmt,dbytes)
            data = np.array(data).reshape(geodict['nrows'],-1)
            data = (data * zscale) + zoffset
        else:
            raise NotImplementedError('Reading of subset of GMT "native" format not yet supported')
        f.close()
        return (data,geodict)

    @classmethod
    def getNetCDFHeader(cls,filename):
        """Get the header information from a GMT NetCDF3 file.
        :param fname:
           File name of GMT NetCDF3 grid
        :returns:
          GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
        """
        cdf = netcdf.netcdf_file(filename)
        geodict = {}
        xvarname = None
        if 'x' in cdf.variables.keys():
            xvarname = 'x'
            yvarname = 'y'
        elif 'lon' in cdf.variables.keys():
            xvarname = 'lon'
            yvarname = 'lat'
        if xvarname is not None:
            xvar = cdf.variables[xvarname].data.copy()
            yvar = cdf.variables[yvarname].data.copy()
            geodict['ncols'] = len(xvar)
            geodict['nrows'] = len(yvar)
            geodict['xmin'] = xvar.min()
            geodict['xmax'] = xvar.max()
            geodict['ymin'] = yvar.min()
            geodict['ymax'] = yvar.max()
            newx = np.linspace(geodict['xmin'],geodict['xmax'],num=geodict['ncols'])
            newy = np.linspace(geodict['ymin'],geodict['ymax'],num=geodict['nrows'])
            geodict['xdim'] = newx[1]-newx[0]
            geodict['ydim'] = newy[1]-newy[0]
        elif 'x_range' in cdf.variables.keys():
            geodict['xmin'] = cdf.variables['x_range'].data[0]
            geodict['xmax'] = cdf.variables['x_range'].data[1]
            geodict['ymin'] = cdf.variables['y_range'].data[0]
            geodict['ymax'] = cdf.variables['y_range'].data[1]
            geodict['ncols'],geodict['nrows'] = cdf.variables['dimension'].data
            #geodict['xdim'],geodict['ydim'] = cdf.variables['spacing'].data
            xvar = np.linspace(geodict['xmin'],geodict['xmax'],num=ncols)
            yvar = np.linspace(geodict['ymin'],geodict['ymax'],num=nrows)
            geodict['xdim'] = xvar[1] - xvar[0]
            geodict['ydim'] = yvar[1] - yvar[0]
        else:
            raise DataSetException('No support for CDF data file with variables: %s' % str(cdf.variables.keys()))
        return (geodict,xvar,yvar)
    
    @classmethod
    def readNetCDF(cls,filename,bounds=None):
        """Read the data and geo-referencing information from a GMT NetCDF3 grid file, subsetting if requested.
        :param filename:
          File name of GMT NetCDF3 grid
        :param bounds:
           Tuple of (xmin,xmax,ymin,ymax)
        :returns:
          Tuple of data (2D numpy array of data, possibly subsetted from file) and geodict (see above).
        """
        geodict,xvar,yvar = cls.getNetCDFHeader(filename)
        cdf = netcdf.netcdf_file(filename)
        if bounds is None:
            data = cdf.variables['z'].data.copy()
        else:
            txmin,txmax,tymin,tymax = bounds
            #we're not doing anything fancy with the data here, just cutting out what we need
            xmin = max(geodict['xmin'],txmin)
            xmax = min(geodict['xmax'],txmax)
            ymin = max(geodict['ymin'],tymin)
            ymax = min(geodict['ymax'],tymax)
            gxmin = geodict['xmin']
            gxmax = geodict['xmax']
            gymin = geodict['ymin']
            gymax = geodict['ymax']
            if xmin == gxmin and xmax == gxmax and ymin == gymin and ymax == gymax:
                data = cdf.variables['z'].data.copy()
            else:
                if xmin > xmax:
                    #cut user's request into two regions - one from the minimum to the
                    #meridian, then another from the meridian to the maximum.
                    (region1,region2) = __createSections((xmin,xmax,ymin,ymax),geodict)
                    (iulx1,iuly1,ilrx1,ilry1) = region1
                    (iulx2,iuly2,ilrx2,ilry2) = region2
                    outcols1 = long(ilrx1-iulx1+1)
                    outcols2 = long(ilrx2-iulx2+1)
                    outcols = long(outcols1+outcols2)
                    outrows = long(ilry1-iuly1+1)
                    cdfarray = BinCDFArray(cdf.variables['z'],len(yvar.data),len(xvar.data))
                    section1 = cdfarray[iuly1:ilry1+1,iulx1:ilrx1+1]
                    section2 = cdfarray[iuly2:ilry2+1,iulx2:ilrx2+1]
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
                else:
                    ixmin = np.abs(xvar-xmin).argmin()
                    ixmax = np.abs(xvar-xmax).argmin()
                    iymin = np.abs(yvar-ymin).argmin()
                    iymax = np.abs(yvar-ymax).argmin()
                    geodict['xmin'] = xvar[ixmin].copy()
                    geodict['xmax'] = xvar[ixmax].copy()
                    geodict['ymin'] = yvar[iymin].copy()
                    geodict['ymax'] = yvar[iymax].copy()
                    data = cdf.variables['z'][iymin:iymax+1,ixmin:ixmax+1].copy()
                    geodict['nrows'],geodict['ncols'] = data.shape

        return (data,geodict)

    @classmethod
    def getHDFHeader(cls,hdffile):
        """Get the header information from a GMT NetCDF4 (HDF) file.
        :param fname:
           File name of GMT NetCDF4 grid
        :returns:
          GeoDict specifying spatial extent, resolution, and shape of grid inside NetCDF file.
        """
        geodict = {}
        f = h5py.File(hdffile,'r')
        if 'x' in f.keys():
            xvarname = 'x'
            yvarname = 'y'
        elif 'lon' in cdf.variables.keys():
            xvarname = 'lon'
            yvarname = 'lat'
        if xvarname is not None:
            xvar = f[xvarname][:]
            yvar = f[yvarname][:]
            geodict['nrows'] = len(yvar)
            geodict['ncols'] = len(xvar)
            geodict['xmin'] = xvar[0]
            geodict['xmax'] = xvar[-1]
            geodict['ymin'] = xvar[0]
            geodict['ymax'] = xvar[-1]
            newx = np.linspace(geodict['xmin'],geodict['xmax'],num=geodict['ncols'])
            newy = np.linspace(geodict['ymin'],geodict['ymax'],num=geodict['nrows'])
            geodict['xdim'] = newx[1]-newx[0]
            geodict['ydim'] = newy[1]-newy[0]
        else:
            geodict['xmin'] = f['x_range'][0]
            geodict['xmax'] = f['x_range'][1]
            geodict['ymin'] = f['y_range'][0]
            geodict['ymax'] = f['y_range'][1]
            geodict['ncols'],geodict['nrows'] = (f['dimension'][0],f['dimension'][1])
            xvar = np.linspace(geodict['xmin'],geodict['xmax'],num=ncols)
            yvar = np.linspace(geodict['ymin'],geodict['ymax'],num=nrows)
            geodict['xdim'] = xvar[1] - xvar[0]
            geodict['ydim'] = yvar[1] - yvar[0]
        f.close()
        return (geodict,xvar,yvar)
        
    @classmethod
    def readHDF(cls,hdffile,bounds=None):
        """Read the data and geo-referencing information from a GMT NetCDF4 (HDF) grid file, subsetting if requested.
        :param hdffile:
          File name of GMT NetCDF4 grid
        :param bounds:
           Tuple of (xmin,xmax,ymin,ymax)
        :returns:
          Tuple of data (2D numpy array of data, possibly subsetted from file) and geodict (see above).
        """
        #need a reproducible way of creating netcdf file in HDF format
        geodict,xvar,yvar = cls.getHDFHeader(hdffile)
        f = h5py.File(hdffile,'r')
        if bounds is None:
            data = zvar[:]
        else:
            txmin,txmax,tymin,tymax = bounds
            #we're not doing anything fancy with the data here, just cutting out what we need
            xmin = max(geodict['xmin'],txmin)
            xmax = min(geodict['xmax'],txmax)
            ymin = max(geodict['ymin'],tymin)
            ymax = min(geodict['ymax'],tymax)
            gxmin = geodict['xmin']
            gxmax = geodict['xmax']
            gymin = geodict['ymin']
            gymax = geodict['ymax']
            if xmin == gxmin and xmax == gxmax and ymin == gymin and ymax == gymax:
                data = f['z'][:]
            else:
                if xmin > xmax:
                    (region1,region2) = __createSections((xmin,xmax,ymin,ymax),geodict)
                    (iulx1,iuly1,ilrx1,ilry1) = region1
                    (iulx2,iuly2,ilrx2,ilry2) = region2
                    outcols1 = long(ilrx1-iulx1+1)
                    outcols2 = long(ilrx2-iulx2+1)
                    outcols = long(outcols1+outcols2)
                    outrows = long(ilry1-iuly1+1)
                    section1 = zvar[iuly1:ilry1+1,iulx1:ilrx1+1]
                    section2 = zvar[iuly2:ilry2+1,iulx2:ilrx2+1]
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
                else:
                    ixmin = np.abs(xvar-xmin).argmin()
                    ixmax = np.abs(xvar-xmax).argmin()
                    iymin = np.abs(yvar-ymin).argmin()
                    iymax = np.abs(yvar-ymax).argmin()
                    geodict['xmin'] = xvar[ixmin].copy()
                    geodict['xmax'] = xvar[ixmax].copy()
                    geodict['ymin'] = yvar[iymin].copy()
                    geodict['ymax'] = yvar[iymax].copy()
                    data = f['z'][iymin:iymax+1,ixmin:ixmax+1].copy()
                    geodict['nrows'],geodict['ncols'] = data.shape
        f.close()
        return (data,geodict)

    def save(self,filename,format='netcdf'):
        """Save a GMTGrid object to a file.
        :param filename:
          Name of desired output file.
        :param format:
          One of 'netcdf','hdf' (binary not yet supported).
        :raises DataSetException:
          When format not one of ('netcdf,'hdf')
        """
        if format not in ['netcdf','hdf']:
            raise DataSetException('Only NetCDF3 and HDF (NetCDF4) output supported at this time.')
        if format == 'netcdf':
            f = netcdf.NetCDFFile(filename,'w')
            m,n = self._data.shape
            xdim = f.createDimension('x',n)
            ydim = f.createDimension('y',m)
            x = f.createVariable('x',np.float64,('x'))
            y = f.createVariable('y',np.float64,('y'))
            x[:] = np.linspace(self._geodict['xmin'],self._geodict['xmax'],self._geodict['ncols'])
            y[:] = np.linspace(self._geodict['ymin'],self._geodict['ymax'],self._geodict['nrows'])
            z = f.createVariable('z',self._data.dtype,('y','x'))
            z[:] = self._data
            f.close()
        elif format == 'hdf':
            f = h5py.File(filename,'w')
            xvar = np.linspace(self._geodict['xmin'],self._geodict['xmax'],self._geodict['ncols'])
            yvar = np.linspace(self._geodict['ymin'],self._geodict['ymax'],self._geodict['nrows'])
            x = f.create_dataset('x',data=xvar)
            y = f.create_dataset('y',data=yvar)
            z = f.create_dataset('z',data=self._data)
            f.close()
    
    @classmethod
    def load(cls,gmtfilename,samplegeodict=None,resample=False,method='linear',doPadding=False,padValue=np.nan):
        """Create a GMTGrid object from a (possibly subsetted, resampled, or padded) GMT grid file.
        :param gmtfilename:
          Name of input file.
        :param samplegeodict:
          GeoDict used to specify subset bounds and resolution (if resample is selected)
        :param resample:
          Boolean used to indicate whether grid should be resampled from the file based on samplegeodict.
        :param method:
          If resample=True, resampling method to use ('nearest','linear','cubic','quintic')
        :param doPadding:
          Boolean used to indicate whether, if samplegeodict is outside bounds of grid, to pad values around the edges.
        :param padValue:
          Value to fill in around the edges if doPadding=True.
        :returns:
          GMTgrid instance (possibly subsetted, padded, or resampled)
        :raises DataSetException:
          * When sample bounds are outside (or too close to outside) the bounds of the grid and doPadding=False.
          * When the input file type is not recognized.
        """
        ftype = cls.getFileType(gmtfilename)
        data = None
        geodict = None
        bounds = None
        samplebounds = None
        if samplegeodict is not None:
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
        
        if ftype == 'binary':
            #we're dealing with a binary "native" GMT grid file
            data,geodict = cls.readGMTNative(gmtfilename,samplebounds)
        elif ftype == 'netcdf':
            data,geodict = cls.readNetCDF(gmtfilename,samplebounds)
        elif ftype == 'hdf':
            data,geodict = cls.readHDF(gmtfilename,samplebounds)
        else:
            raise DataSetException('File type "%s" cannot be read.' % ftype)
        
        if doPadding:
            #up to this point, all we've done is either read in the whole file or cut out (along existing
            #boundaries) the section of data we want.  Now we do padding as necessary.
            #_getPadding is a class method inherited from Grid (our grandparent)
            leftpad,rightpad,bottompad,toppad,geodict = super(Grid2D,cls)._getPadding(geodict,samplebounds,padValue)
            data = np.hstack((leftpad,data))
            data = np.hstack((data,rightpad))
            data = np.vstack((toppad,data))
            data = np.vstack((data,bottompad))
        #if the user asks to resample, take the (possibly cut and padded) data set, and resample
        #it using the Grid2D super class
        if resample:
            grid = Grid2D(data,geodict)
            grid.interpolateToGrid(samplegeodict,method=method)
            data = grid.getData()
            geodict = grid.getGeoDict()
        return cls(data,geodict)
            
        
class BinCDFArray(object):
    def __init__(self,array,nrows,ncols):
        self.array = array
        self.nrows = nrows
        self.ncols = ncols

    def __getitem__(self,*args):
        """Allows slicing of CDF data array in the same way as a numpy array."""
        if len(args) == 1 and isinstance(args[0][0],int):
            #user has passed in a tuple of row,col - they only want one value
            row = args[0][0]
            col = args[0][1]
            nrows = self.nrows
            ncols = self.ncols
            if row < 0 or row > nrows-1:
                raise Exception,"Row index out of bounds"
            if col < 0 or col > ncols-1:
                raise Exception,"Row index out of bounds"
            idx = ncols * row + col
            offset = 0
            return self.array[idx]

        if len(args) == 1 and isinstance(args[0][0],slice): #they want a non-scalar subset of the data
            nrows = self.nrows
            ncols = self.ncols
            #calculate offset to first data element
            key1 = args[0][0]
            key2 = args[0][1]
            rowstart = key1.start
            rowend = key1.stop
            rowstep = key1.step
            colstart = key2.start
            colend = key2.stop
            colstep = key2.step
            
            if rowstep is None:
                rowstep = 1
            if colstep is None:
                colstep = 1

            #error checking
            if rowstart < 0 or rowstart > nrows-1:
                raise Exception,"Row index out of bounds"
            if rowend < 0 or rowend > nrows:
                raise Exception,"Row index out of bounds"
            if colstart < 0 or colstart > ncols-1:
                raise Exception,"Col index out of bounds"
            if colend < 0 or colend > ncols:
                raise Exception,"Col index out of bounds"

            colcount = (colend-colstart)
            rowcount = (rowend-rowstart)
            outrows = np.ceil(rowcount/rowstep)
            outcols = np.ceil(colcount/colstep)
            data = np.zeros([outrows,outcols],dtype=self.dtype)
            outrow = 0
            for row in range(int(rowstart),int(rowend),int(rowstep)):
                #just go to the beginning of the row, we're going to read in the whole line
                idx = ncols*row 
                offset = self.dwidth*idx #beginning of row
                line = self.array[idx:idx+ncols]
                data[outrow,:] = line[colstart:colend:colstep]
                outrow = outrow+1
                
        else:
            raise Exception, "Unsupported __getitem__ input %s" % (str(key))
        return(data)

def _save_test():
    try:
        print
        print 'Testing saving and loading from different formats:'
        #make a sample data set
        data = np.arange(0,16).reshape(4,4)
        data = data.astype(np.int32) #arange gives int64 by default, not supported by netcdf3
        m,n = data.shape
        geodict = {'nrows':m,
                   'ncols':n,
                   'xmin':0.5,
                   'xmax':3.5,
                   'ymin':0.5,
                   'ymax':3.5,
                   'xdim':1.0,
                   'ydim':1.0}
        gmtgrid = GMTGrid(data,geodict)
        print gmtgrid

        #save it as netcdf3
        gmtgrid.save('test.grd',format='netcdf')
        gmtgrid2 = GMTGrid.load('test.grd')
        print
        print gmtgrid2

        #save it as HDF
        gmtgrid.save('test.grd',format='hdf')
        gmtgrid2 = GMTGrid.load('test.grd')
        print
        print gmtgrid2
    except:
        pass
    os.remove('test.grd')

def _pad_test():
    try:
        print
        print 'Test padding data with null values:'
        data = np.arange(0,16).reshape(4,4)
        data = data.astype(np.int32) #arange gives int64 by default, not supported by netcdf3
        m,n = data.shape
        geodict = {'nrows':m,
                   'ncols':n,
                   'xmin':0.5,
                   'xmax':3.5,
                   'ymin':0.5,
                   'ymax':3.5,
                   'xdim':1.0,
                   'ydim':1.0}
        gmtgrid = GMTGrid(data,geodict)
        gmtgrid.save('test.grd',format='netcdf')
        print gmtgrid

        newdict = {'xmin':-0.5,'xmax':4.5,'ymin':-0.5,'ymax':4.5,'xdim':1.0,'ydim':1.0}
        gmtgrid2 = GMTGrid.load('test.grd',samplegeodict=newdict,doPadding=True)
        print gmtgrid2
        print gmtgrid2._data
    except:
        pass
    os.remove('test.grd')

def _subset_test():
    try:
        print
        print 'Test subsetting data:'
        gmtgrid = createSampleGrid(4)
        gmtgrid.save('test.grd',format='netcdf')
        print gmtgrid

        newdict = {'xmin':1.5,'xmax':2.5,'ymin':1.5,'ymax':2.5,'xdim':1.0,'ydim':1.0}
        gmtgrid2 = GMTGrid.load('test.grd',samplegeodict=newdict)
        print gmtgrid2
        print gmtgrid2._data

        gmtgrid.save('test.grd',format='hdf')
        print gmtgrid
        gmtgrid3 = GMTGrid.load('test.grd',samplegeodict=newdict)
        print gmtgrid3
        print gmtgrid3._data
    except:
        pass

    os.remove('test.grd')

def _resample_test():
    try:
        print
        print 'Test subsetting data:'
        gmtgrid = createSampleGrid(6)
        gmtgrid.save('test.grd',format='netcdf')
        print gmtgrid
        print gmtgrid._data

        # newdict = {'xmin':2.0,'xmax':4.0,'ymin':2.0,'ymax':4.0,'xdim':1.0,'ydim':1.0}
        # newdict = Grid.fillGeoDict(newdict)
        # gmtgrid3 = GMTGrid.load('test.grd',samplegeodict=newdict,resample=True)
        # print gmtgrid3
        # print gmtgrid3._data

        gmtgrid = createSampleGrid(4)
        gmtgrid.save('test.grd',format='netcdf')
        newdict = {'xmin':-1.0,'xmax':5.0,'ymin':-1.0,'ymax':5.0,'xdim':1.0,'ydim':1.0}
        newdict = Grid.fillGeoDict(newdict)
        gmtgrid3 = GMTGrid.load('test.grd',samplegeodict=newdict,resample=True,doPadding=True)
        print gmtgrid3
        print gmtgrid3._data
    except:
        pass

    os.remove('test.grd')
    
if __name__ == '__main__':
    _save_test()
    _pad_test()
    _subset_test()
    _resample_test()
