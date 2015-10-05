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
#1) Write a test for 180 meridian crossing grids

def createSampleGrid(M,N):
    """Used for internal testing - create an NxN grid with lower left corner at 0.5,0.5, xdim/ydim = 1.0
    :param M:
       Number of rows in output grid
    :param N:
       Number of columns in output grid
    :returns:
       GMTGrid object where data values are an NxN array of values from 0 to N-squared minus 1, and geodict
       lower left corner is at 0.5/0.5 and cell dimensions are 1.0.
    """
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

        #We are going to represent all grids internally as grid-line registered
        #The difference between pixel and gridline-registered grids is depicted well here:
        #http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#grid-registration-the-r-option
        if offset == 1:
            geodict['xmin'] += geodict['xdim']/2.0
            geodict['xmax'] -= geodict['xdim']/2.0
            geodict['ymin'] += geodict['ydim']/2.0
            geodict['ymax'] -= geodict['ydim']/2.0
            
        return geodict

            
    @classmethod
    def readGMTNative(cls,fname,bounds=None,firstColumnDuplicated=False,fmt=None):
        """Read the data and geo-referencing information from a GMT native grid file, subsetting if requested.
        http://gmt.soest.hawaii.edu/doc/5.1.2/GMT_Docs.html#native-binary-grid-files
        :param fname:
          File name of GMT native grid
        :param bounds:
           Tuple of (xmin,xmax,ymin,ymax)
        :param firstColumnDuplicated:
           Boolean - is this a file where the last column of data is the same as the first (for grids that span entire globe).
        :param fmt: 
           Data width, one of: 
             - 'i' (16 bit signed integer)
             - 'l' (32 bit signed integer)
             - 'f' (32 bit float)
             - 'd' (64 bit float)
          Strictly speaking, this is only necessary when the data file is 32 bit float or 32 bit integer, as there is
          no *sure* way to tell from the header or data which data type is contained in the file.  If fmt is None, then
          the code will try to guess as best it can from the data whether it is integer or floating point data. Caveat emptor!
        :returns:
          Tuple of data (2D numpy array of data, possibly subsetted from file) and geodict (see above).
        :raises NotImplementedError:
          For any bounds not None (we'll get to it eventually!)
        """
        #Given that we can't automatically distinguish between 32 bit ints and 32 bit floats, we'll use
        #this value as a cutoff for "detecting" when a value read from disk as a float has an "unreasonably" 
        #high exponent value.  This is NOT guaranteed to work - use the fmt keyword if you want to be sure.
        MAX_FLOAT_EXP = 30
        fsize = os.path.getsize(fname)
        datalen = fsize-892
        f = open(fname,'rb')
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
        #We are going to represent all grids internally as grid-line registered
        #The difference between pixel and gridline-registered grids is depicted well here:
        #http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#grid-registration-the-r-option
        npixels = geodict['nrows']*geodict['ncols']
        lenshort = npixels*2
        lenfloat = npixels*4
        lendouble = npixels*8
        if fmt is None:
            if datalen == lenshort:
                fmt = 'h'
            if datalen == lendouble:
                fmt = 'd'
            if datalen == lenfloat: #let's try to guess whether this is float or int
                fpos = f.tell()
                #read 1 byte, check to see if it's nan or 0 - if it is, then we definitely have a float
                dbytes = struct.unpack('f',f.read(4))[0]
                while dbytes == 0.0:
                    dbytes = struct.unpack('f',f.read(4))[0]
                f.seek(fpos) #go back to where we were
                if np.isnan(dbytes):
                    fmt = 'f'
                elif int(np.abs(np.log10(dbytes))) > MAX_FLOAT_EXP: #does this have a crazy large exponent?
                    fmt = 'i'
                else:
                    fmt = 'f'
        if offset == 1:
            geodict['xmin'] += geodict['xdim']/2.0
            geodict['xmax'] -= geodict['xdim']/2.0
            geodict['ymin'] += geodict['ydim']/2.0
            geodict['ymax'] -= geodict['ydim']/2.0
        if bounds is None:
            sfmt = '%i%s' % (geodict['ncols']*geodict['nrows'],fmt)
            dwidths = {'h':2,'i':4,'f':4,'d':8}
            dwidth = dwidths[fmt]
            dbytes = f.read(geodict['ncols']*geodict['nrows']*dwidth)
            data = np.array(struct.unpack(sfmt,dbytes))
            #data = np.array(data).reshape(geodict['nrows'],-1)
            data.shape = (geodict['nrows'],geodict['ncols'])
            if zscale != 1.0 or offset != 0.0:
                data = (data * zscale) + zoffset
            if firstColumnDuplicated:
                data = data[:,0:-1]
                geodict['xmax'] -= geodict['xdim']
        else:
            #make sure to handle wrap-around data where first column is duplicated
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
        registration = 'gridline'
        if hasattr(cdf,'node_offset') and getattr(cdf,'node_offset') == 1:
            registration = 'pixel'
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

        #We are going to represent all grids internally as grid-line registered
        #The difference between pixel and gridline-registered grids is depicted well here:
        #http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#grid-registration-the-r-option
        if registration == 'pixel':
            geodict['xmin'] += geodict['xdim']/2.0
            geodict['xmax'] -= geodict['xdim']/2.0
            geodict['ymin'] += geodict['ydim']/2.0
            geodict['ymax'] -= geodict['ydim']/2.0
            
        return (geodict,xvar,yvar)
    
    @classmethod
    def readNetCDF(cls,filename,bounds=None,firstColumnDuplicated=False):
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
            if firstColumnDuplicated:
                data = data[:,0:-1]
                geodict['xmax'] -= geodict['xdim']
        else:
            txmin,txmax,tymin,tymax = bounds
            #we're not doing anything fancy with the data here, just cutting out what we need
            xmin = max(geodict['xmin'],txmin)
            xmax = min(geodict['xmax'],txmax)
            ymin = max(geodict['ymin'],tymin)
            ymax = min(geodict['ymax'],tymax)
            #these are the bounds of the whole file
            gxmin = geodict['xmin']
            gxmax = geodict['xmax']
            gymin = geodict['ymin']
            gymax = geodict['ymax']
            xdim = geodict['xdim']
            ydim = geodict['ydim']
            gnrows = geodict['nrows']
            gncols = geodict['ncols']
            if xmin == gxmin and xmax == gxmax and ymin == gymin and ymax == gymax:
                data = cdf.variables['z'].data.copy()
                if firstColumnDuplicated:
                    data = data[:,0:-1]
                    geodict['xmax'] -= geodict['xdim']
            else:
                if xmin > xmax:
                    #cut user's request into two regions - one from the minimum to the
                    #meridian, then another from the meridian to the maximum.
                    (region1,region2) = cls._createSections((xmin,xmax,ymin,ymax),geodict,firstColumnDuplicated)
                    (iulx1,iuly1,ilrx1,ilry1) = region1
                    (iulx2,iuly2,ilrx2,ilry2) = region2
                    outcols1 = long(ilrx1-iulx1)
                    outcols2 = long(ilrx2-iulx2)
                    outcols = long(outcols1+outcols2)
                    outrows = long(ilry1-iuly1)
                    #zvar = cdf.variables['z']
                    section1 = cdf.variables['z'][iuly1:ilry1,iulx1:ilrx1].copy()
                    section2 = cdf.variables['z'][iuly2:ilry2,iulx2:ilrx2].copy()
                    data = np.concatenate((section1,section2),axis=1)
                    outrows,outcols = data.shape
                    xmin = (gxmin + iulx1*xdim)
                    ymax = gymax - iuly1*ydim
                    xmax = gxmin + (ilrx2-1)*xdim
                    ymin = gymin + (gnrows-ilry1)*ydim
                    geodict['xmin'] = xmin
                    geodict['xmax'] = xmax + 360
                    geodict['ymin'] = ymin
                    geodict['ymax'] = ymax
                    geodict['nrows'],geodict['ncols'] = data.shape
                else:
                    ixmin = np.abs(xvar-xmin).argmin()
                    ixmax = np.abs(xvar-xmax).argmin()
                    iymin = np.abs(yvar-ymin).argmin()
                    iymax = np.abs(yvar-ymax).argmin()
                    if firstColumnDuplicated:
                        ixmax -= 1
                    geodict['xmin'] = xvar[ixmin].copy()
                    geodict['xmax'] = xvar[ixmax].copy()
                    geodict['ymin'] = yvar[iymin].copy()
                    geodict['ymax'] = yvar[iymax].copy()
                    data = cdf.variables['z'][iymin:iymax+1,ixmin:ixmax+1].copy()
                    geodict['nrows'],geodict['ncols'] = data.shape
        cdf.close()
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
        registration = 'gridline'
        if f.get('node_offset') is not None and f.attrs['node_offset'][0] == 1:
            registration = 'pixel'
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

        #We are going to represent all grids internally as grid-line registered
        #The difference between pixel and gridline-registered grids is depicted well here:
        #http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#grid-registration-the-r-option
        if registration == 'pixel':
            geodict['xmin'] += geodict['xdim']/2.0
            geodict['xmax'] -= geodict['xdim']/2.0
            geodict['ymin'] += geodict['ydim']/2.0
            geodict['ymax'] -= geodict['ydim']/2.0
        f.close()
        return (geodict,xvar,yvar)
        
    @classmethod
    def readHDF(cls,hdffile,bounds=None,firstColumnDuplicated=False):
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
        zvar = f['z']
        if bounds is None:
            data = zvar[:]
            if firstColumnDuplicated:
                data = data[:,0:-1]
                geodict['xmax'] -= geodict['xdim']
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
                if firstColumnDuplicated:
                    data = data[:,0:-1]
                    geodict['xmax'] -= geodict['xdim']
            else:
                if xmin > xmax:
                    (region1,region2) = __createSections((xmin,xmax,ymin,ymax),geodict,firstColumnDuplicated)
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
                    if firstColumnDuplicated:
                        ixmax -= 1
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
          One of 'netcdf','hdf' or 'binary'.
        :raises DataSetException:
          When format not one of ('netcdf,'hdf','binary')
        """
        if format not in ['netcdf','hdf','binary']:
            raise DataSetException('Only NetCDF3, HDF (NetCDF4), and GMT native output are supported.')
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
        elif format == 'binary':
            f = open(filename,'w')
            f.write(struct.pack('I',self._geodict['ncols']))
            f.write(struct.pack('I',self._geodict['nrows']))
            f.write(struct.pack('I',0)) #gridline registration
            f.write(struct.pack('d',self._geodict['xmin']))
            f.write(struct.pack('d',self._geodict['xmax']))
            f.write(struct.pack('d',self._geodict['ymin']))
            f.write(struct.pack('d',self._geodict['ymax']))
            f.write(struct.pack('d',self._data.min()))
            f.write(struct.pack('d',self._data.max()))
            f.write(struct.pack('d',self._geodict['xdim']))
            f.write(struct.pack('d',self._geodict['ydim']))
            f.write(struct.pack('d',1.0)) #scale factor to multiply data by
            f.write(struct.pack('d',0.0)) #offfset to add to data
            f.write(struct.pack('80s','X units (probably degrees)'))
            f.write(struct.pack('80s','Y units (probably degrees)'))
            f.write(struct.pack('80s','Z units unknown'))
            f.write(struct.pack('80s','')) #title
            f.write(struct.pack('320s','Created with GMTGrid() class, a product of the NEIC.')) #command
            f.write(struct.pack('160s','')) #remark
            if self._data.dtype not in [np.int16,np.int32,np.float32,np.float64]:
                raise DataSetException('Data type of "%s" is not supported by the GMT native format.' % str(self._data.dtype))
            fpos1 = f.tell()
            self._data.tofile(f)
            fpos2 = f.tell()
            bytesout = fpos2 - fpos1
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
        firstColumnDuplicated = False
        if samplegeodict is not None:
            bounds = (samplegeodict['xmin'],samplegeodict['xmax'],samplegeodict['ymin'],samplegeodict['ymax'])
            samplebounds = bounds
            #if the user wants resampling, we can't just read the bounds they asked for, but instead
            #go outside those bounds.  if they asked for padding and the input bounds exceed the bounds
            #of the file, then we can pad.  If they *didn't* ask for padding and input exceeds, raise exception.
            if resample:
                PADFACTOR = 2 #how many cells will we buffer out for resampling?
                filegeodict = cls.getFileGeoDict(gmtfilename)
                xdim = filegeodict['xdim']
                ydim = filegeodict['ydim']
                fbounds = (filegeodict['xmin'],filegeodict['xmax'],filegeodict['ymin'],filegeodict['ymax'])
                hasMeridianWrap = False
                if fbounds[0] == fbounds[1]-360:
                    firstColumnDuplicated = True
                if firstColumnDuplicated or np.abs(fbounds[0]-(fbounds[1]-360)) == xdim:
                    hasMeridianWrap = True
                isOutside = False
                #make a bounding box that is PADFACTOR number of rows/cols greater than what the user asked for
                rbounds = [bounds[0]-xdim*PADFACTOR,bounds[1]+xdim*PADFACTOR,bounds[2]-ydim*PADFACTOR,bounds[3]+ydim*PADFACTOR]
                #compare that bounding box to the file bounding box
                if not hasMeridianWrap:
                    if fbounds[0] > rbounds[0] or fbounds[1] < rbounds[1] or fbounds[2] > rbounds[2] or fbounds[3] < rbounds[3]:
                        isOutside = True
                else:
                    if fbounds[2] > rbounds[2] or fbounds[3] < rbounds[3]:
                        isOutside = True
                if isOutside:
                    if doPadding==False:
                        raise DataSetException('Cannot resample data given input bounds, unless doPadding is set to True.')
                    else:
                        samplebounds = rbounds
                else:
                    samplebounds = rbounds
        
        if ftype == 'binary':
            #we're dealing with a binary "native" GMT grid file
            data,geodict = cls.readGMTNative(gmtfilename,samplebounds,firstColumnDuplicated)
        elif ftype == 'netcdf':
            data,geodict = cls.readNetCDF(gmtfilename,samplebounds,firstColumnDuplicated)
        elif ftype == 'hdf':
            data,geodict = cls.readHDF(gmtfilename,samplebounds,firstColumnDuplicated)
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
            if samplegeodict['xmin'] > samplegeodict['xmax']:
                samplegeodict['xmax'] += 360
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
        print 'Testing saving and loading to/from NetCDF3...'
        #make a sample data set
        gmtgrid = createSampleGrid(4,4)

        #save it as netcdf3
        gmtgrid.save('test.grd',format='netcdf')
        gmtgrid2 = GMTGrid.load('test.grd')
        np.testing.assert_almost_equal(gmtgrid._data,gmtgrid2._data)
        print 'Passed saving and loading to/from NetCDF3.'

        print 'Testing saving and loading to/from NetCDF4 (HDF)...'
        #save it as HDF
        gmtgrid.save('test.grd',format='hdf')
        gmtgrid3 = GMTGrid.load('test.grd')
        np.testing.assert_almost_equal(gmtgrid._data,gmtgrid3._data)
        print 'Passed saving and loading to/from NetCDF4 (HDF)...'

        print 'Testing saving and loading to/from NetCDF4 (HDF)...'
        gmtgrid.save('test.grd',format='binary')
        gmtgrid4 = GMTGrid.load('test.grd')
        np.testing.assert_almost_equal(gmtgrid._data,gmtgrid4._data)
        print 'Passed saving and loading to/from GMT native...'
    except AssertionError,error:
        print 'Failed padding test:\n %s' % error
    os.remove('test.grd')

def _pad_test():
    try:
        print 'Test padding data with null values...'
        gmtgrid = createSampleGrid(4,4)
        gmtgrid.save('test.grd',format='netcdf')

        newdict = {'xmin':-0.5,'xmax':4.5,'ymin':-0.5,'ymax':4.5,'xdim':1.0,'ydim':1.0}
        gmtgrid2 = GMTGrid.load('test.grd',samplegeodict=newdict,doPadding=True)
        output = np.array([[np.nan,np.nan,np.nan,np.nan,np.nan,np.nan],
                           [np.nan,0.0,1.0,2.0,3.0,np.nan],
                           [np.nan,4.0,5.0,6.0,7.0,np.nan],
                           [np.nan,8.0,9.0,10.0,11.0,np.nan],
                           [np.nan,12.0,13.0,14.0,15.0,np.nan],
                           [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]])
        np.testing.assert_almost_equal(gmtgrid2._data,output)
        print 'Passed padding data.'
    except AssertionError,error:
        print 'Failed padding test:\n %s' % error
    os.remove('test.grd')

def _subset_test():
    try:
        print 'Test subsetting data without padding...'
        gmtgrid = createSampleGrid(4,4)
        gmtgrid.save('test.grd',format='netcdf')

        newdict = {'xmin':1.5,'xmax':2.5,'ymin':1.5,'ymax':2.5,'xdim':1.0,'ydim':1.0}
        gmtgrid2 = GMTGrid.load('test.grd',samplegeodict=newdict)
        output = np.array([[5,6],
                           [9,10]])
        np.testing.assert_almost_equal(gmtgrid2._data,output)
        print 'Passed subsetting data without padding.'

        # gmtgrid.save('test.grd',format='hdf')
        # print gmtgrid
        # gmtgrid3 = GMTGrid.load('test.grd',samplegeodict=newdict)
        # print gmtgrid3
        # print gmtgrid3._data
    except AssertionError,error:
        print 'Failed resample test:\n %s' % error

    os.remove('test.grd')

def _resample_test():
    try:
        print 'Test resampling data without padding...'
        gmtgrid = createSampleGrid(7,7)
        gmtgrid.save('test.grd',format='netcdf')

        newdict = {'xmin':3.0,'xmax':4.0,'ymin':3.0,'ymax':4.0,'xdim':1.0,'ydim':1.0}
        newdict = Grid.fillGeoDict(newdict)
        gmtgrid3 = GMTGrid.load('test.grd',samplegeodict=newdict,resample=True)
        output = np.array([[20,21],
                           [27,28]])
        np.testing.assert_almost_equal(gmtgrid3._data,output)
        print 'Passed resampling without padding.'

        print 'Test resampling data with padding...'
        gmtgrid = createSampleGrid(4,4)
        gmtgrid.save('test.grd',format='netcdf')
        newdict = {'xmin':0.0,'xmax':4.0,'ymin':0.0,'ymax':4.0,'xdim':1.0,'ydim':1.0}
        newdict = Grid.fillGeoDict(newdict)
        gmtgrid3 = GMTGrid.load('test.grd',samplegeodict=newdict,resample=True,doPadding=True)
        output = np.array([[np.nan,np.nan,np.nan,np.nan,np.nan],
                           [np.nan,2.5,3.5,4.5,np.nan],
                           [np.nan,6.5,7.5,8.5,np.nan],
                           [np.nan,10.5,11.5,12.5,np.nan],
                           [np.nan,np.nan,np.nan,np.nan,np.nan]])
        np.testing.assert_almost_equal(gmtgrid3._data,output)
        print 'Passed resampling without padding.'
    except AssertionError,error:
        print 'Failed resample test:\n %s' % error

    os.remove('test.grd')

def _meridian_test():
    try:
        print 'Testing resampling of global grid where sample crosses 180/-180 meridian...'
        data = np.arange(0,84).astype(np.int32).reshape(7,12)
        geodict = {'xmin':-180.0,'xmax':150.0,'ymin':-90.0,'ymax':90.0,'xdim':30,'ydim':30,'nrows':7,'ncols':12}
        gmtgrid = GMTGrid(data,geodict)
        gmtgrid.save('test.grd')

        sampledict = {'xmin':105,'xmax':-105,'ymin':-15.0,'ymax':15.0,'xdim':30.0,'ydim':30.0,'nrows':2,'ncols':5}
        gmtgrid5 = GMTGrid.load('test.grd',samplegeodict=sampledict,resample=True,doPadding=True)

        output = np.array([[ 39.5,40.5,35.5,30.5,31.5,32.5],
                           [ 51.5,52.5,47.5,42.5,43.5,44.5]])
        #output = np.random.rand(2,6) #this will fail assertion test
        np.testing.assert_almost_equal(gmtgrid5._data,output)
        print 'Done testing resampling of global grid where sample crosses 180/-180 meridian.'

        print 'Testing resampling of global grid where sample crosses 180/-180 meridian and first column is duplicated by last...'
        data = np.arange(0,84).astype(np.int32).reshape(7,12)
        data = np.hstack((data,data[:,0].reshape(7,1)))
        geodict = {'xmin':-180.0,'xmax':180.0,'ymin':-90.0,'ymax':90.0,'xdim':30,'ydim':30,'nrows':7,'ncols':13}
        gmtgrid = GMTGrid(data,geodict)
        gmtgrid.save('test.grd')

        sampledict = {'xmin':105,'xmax':-105,'ymin':-15.0,'ymax':15.0,'xdim':30.0,'ydim':30.0,'nrows':2,'ncols':5}
        gmtgrid5 = GMTGrid.load('test.grd',samplegeodict=sampledict,resample=True,doPadding=True)

        output = np.array([[ 39.5,40.5,35.5,30.5,31.5,32.5],
                           [ 51.5,52.5,47.5,42.5,43.5,44.5]])
        #output = np.random.rand(2,6) #this will fail assertion test
        np.testing.assert_almost_equal(gmtgrid5._data,output)
        print 'Done testing resampling of global grid where sample crosses 180/-180 meridian and first column is duplicated by last...'
        
    except AssertionError,error:
        print 'Failed meridian test:\n %s' % error
    os.remove('test.grd')
if __name__ == '__main__':
    _save_test()
    _pad_test()
    _subset_test()
    _resample_test()
    _meridian_test()
