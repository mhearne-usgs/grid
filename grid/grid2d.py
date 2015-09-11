#!/usr/bin/env python

from gridbase import Grid
import numpy as np
from scipy import interpolate
import abc

class GridException(Exception):
    """
    Class to represent errors in the Grid class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Grid2D(Grid):
    """
    A partially abstract class to represent 2D lat/lon gridded datasets. Some basic methods
    are implemented here, enough so that all functions of working with data (aside from loading and saving)
    can be used with this class.  Grids are assumed to be pixel-registered - that is, grid coordinates
    represent the value at the *center* of the cells.
    """
    reqfields = set(['xmin','xmax','ymin','ymax','xdim','ydim'])
    def __init__(self,data=None,geodict=None):
        """
        Construct a Grid object.
        
        :param data: 
            A 2D numpy array
        :param geodict: 
            A dictionary containing the following fields:
             - xmin Longitude minimum (decimal degrees) (Center of upper left cell)
             - xmax Longitude maximum (decimal degrees) (Center of upper right cell)
             - ymin Longitude minimum (decimal degrees) (Center of lower left cell)
             - ymax Longitude maximum (decimal degrees) (Center of lower right cell)
             - xdim Cell width (decimal degrees)
             - ydim Cell height (decimal degrees)
             - (optional) nrows Number of rows of input data (will be adjusted if incorrect)
             - (optional) ncols Number of columns of input data (will be adjusted if incorrect)
        :returns:
            A Grid object.  Internal representation of geodict input will have nrows/ncols fields added.
            If these 
        """
        if data is not None and geodict is not None:
            #complain if data is not 2D (squeeze 1d dimensions out)
            dims = data.shape
            if len(dims) != 2:
                raise GridException('Grid data must be 2D.  Input data has shape of %s' % str(data.shape))
            #complain if geodict does not have all required fields
            if not set(geodict.keys()).issuperset(self.reqfields):
                raise GridException('Grid data must be 2D.  Input data has shape of %s' % str(data.shape))
            geodict['nrows'],geodict['ncols'] = data.shape
            self._geodict = geodict.copy()
            self._data = data.copy()
        else:
            self._data = None
            self._geodict = None

    #This should be a @classmethod in subclasses
    @abc.abstractmethod
    def load(filename,bounds=None,resample=False,padValue=None):
        raise NotImplementedError('Load method not implemented in base class')

    #TODO: Figure out how format-specific attributes will be handled (ShakeMap, for example)
    #This should be a @classmethod in subclasses
    @abc.abstractmethod
    def save(self,filename): #would we ever want to save a subset of the data?
        raise NotImplementedError('Save method not implemented in base class')
    
    def getData(self):
        """
        Return a reference to the data inside the Grid
        :returns:
          A reference to a 2D numpy array.
        """
        return self._data #should we return a copy of the data?

    def getGeoDict(self):
        """
        Return a reference to the geodict inside the Grid
        :returns:
          A reference to a dictionary (see constructor).
        """
        return self._geodict #should we return a copy of the geodict?

    def getBounds(self):
        """
        Return the lon/lat range of the data.
        
        :returns:
           Tuple of (lonmin,lonmax,latmin,latmax)
        """
        return (self._geodict['xmin'],self._geodict['xmax'],self._geodict['ymin'],self._geodict['ymax'])

    def trim(self,bounds,resample=False,method='linear'):
        """
        Trim data to a smaller set of bounds, resampling if requested.  If not resampling,
        data will be trimmed to smallest grid boundary possible.
        
        :param bounds:
           Tuple of (lonmin,lonmax,latmin,latmax)
        :param resample:
           Boolean indicating whether the data should be resampled to *exactly* match input bounds.
        :param method:
           If resampling, method used, one of ('linear','nearest','cubic','quintic')
        """
        xmin,xmax,ymin,ymax = bounds
        if not resample:
            uly,ulx = self.getRowCol(ymax,xmin,returnFloat=True)
            lry,lrx = self.getRowCol(ymin,xmax,returnFloat=True)
            uly = int(np.floor(uly))
            ulx = int(np.ceil(ulx))
            lrx = int(np.floor(lrx))
            lry = int(np.ceil(lry))
            self._data = self._data[uly:lry,ulx:lrx]
            newymax,newxmin = self.getLatLon(uly,ulx)
            newymin,newxmax = self.getLatLon(lry,lrx)
            self._geodict['xmin'] = newxmin
            self._geodict['xmax'] = newxmax
            self._geodict['ymin'] = newymin
            self._geodict['ymax'] = newymax
            self._geodict['nrows'],self._geodict['ncols'] = self._data.shape
        else:
            xdim = self._geodict['xdim']
            ydim = self._geodict['ydim']
            indict = {'xmin':xmin,'xmax':xmax,'ymin':ymin,'ymax':ymax,'xdim':xdim,'ydim':ydim}
            ncols = len(np.arange(xmin,xmax+xdim,xdim))
            nrows = len(np.arange(ymin,ymax+ydim,ydim))
            indict['nrows'] = nrows
            indict['ncols'] = ncols
            self.interpolateToGrid(indict,method=method)

    def getValue(self,lat,lon,method='nearest',default=None): #return nearest neighbor value
        """Return numpy array at given latitude and longitude (using nearest neighbor).
        
        :param lat: 
           Latitude (in decimal degrees) of desired data value.
        :param lon: 
           Longitude (in decimal degrees) of desired data value.
        :param method:
           Interpolation method, one of ('nearest','linear','cubic','quintic')
        :param default:
           Default value to return when lat/lon is outside of grid bounds.
        :return: 
           Value at input latitude,longitude position.
        """
        if method == 'nearest':
            row,col = self.getRowCol(lat,lon)
        else:
            row,col = self.getRowCol(lat,lon,returnFloat=True)
        nrows,ncols = self._data.shape
        if (row < 0).any() or (row > nrows-1).any() or (col < 0).any() or (col > ncols-1).any():
            if default is None:
                msg = 'One of more of your lat/lon values is outside Grid boundaries: %s' % (str(self.getRange()))
                raise GridException(msg)
            value = np.ones_like(lat)*default
            return value
        if method == 'nearest':
            return self._data[row,col]
        else:
            raise NotImplementedError('getValue method "%s" not implemented yet' % method)

    def getLatLon(self,row,col):
        """Return geographic coordinates (lat/lon decimal degrees) for given data row and column.
        
        :param row: 
           Row dimension index into internal data array.
        :param col: 
           Column dimension index into internal data array.
        :returns: 
           Tuple of latitude and longitude.
        """
        ulx = self._geodict['xmin']
        uly = self._geodict['ymax']
        xdim = self._geodict['xdim']
        ydim = self._geodict['ydim']
        lon = ulx + col*xdim
        lat = uly - row*ydim
        return (lat,lon)

    def getRowCol(self,lat,lon,returnFloat=False):
        """Return data row and column from given geographic coordinates (lat/lon decimal degrees).
        
        :param lat: 
           Input latitude.
        :param lon: 
           Input longitude.
        :param returnFloat: 
           Boolean indicating whether floating point row/col coordinates should be returned.
        :returns: 
           Tuple of row and column.
        """
        ulx = self._geodict['xmin']
        uly = self._geodict['ymax']
        xdim = self._geodict['xdim']
        ydim = self._geodict['ydim']
        #check to see if we're in a scenario where the grid crosses the meridian
        if self._geodict['xmax'] < ulx and lon < ulx:
            lon += 360
        col = (lon-ulx)/xdim
        row = (uly-lat)/ydim
        if returnFloat:
            return (row,col)
        
        return (np.floor(row).astype(int),np.floor(col).astype(int))
            
    def _getInterpCoords(self,geodict):
        #get the cell coordinates of the grid we want to interpolate to
        dims = self._data.shape
        nrows1 = dims[0]
        ncols1 = dims[1]
        ulx1 = self._geodict['xmin']
        uly1 = self._geodict['ymax']
        xdim1 = self._geodict['xdim']
        ydim1 = self._geodict['ydim']
        
        #extract the geographic information about the grid we're sampling to
        nrows = geodict['nrows']
        ncols = geodict['ncols']
        ulx = geodict['xmin']
        uly = geodict['ymax']
        xdim = geodict['xdim']
        ydim = geodict['ydim']

        #make sure that base grid is completely contained within the grid to be
        #resampled
        lry = geodict['ymin']
        lrx = geodict['xmax']
        lry1 = self._geodict['ymin']
        lrx1 = self._geodict['xmax']

        if (lry < lry1 or lrx > lrx1):
            raise GridException('Error:  Base grid is not completely contained by resampling grid.')
        
        #establish the geographic coordinates of the centers of our pixels...
        #all geostruct grids are what GMT calls "pixel-registered", that is
        #the upper left hand corner position is the geographic position of the 
        #center of the cell, not the upper-left corner of it.
        starty = lry
        endx = lrx
        endy = uly
        startx = ulx

        #we need to handle the meridian crossing here...
        if startx > endx:
            endx += 360
            ulx1 += 360

        gxi = np.arange(startx,endx,xdim,dtype=np.float64)
        gyi = np.arange(endy,starty,-ydim,dtype=np.float64)
        
        #we may wind up with an array that is one shorter than we need...
        #in this case, append the last value.
        if len(gxi) < ncols:
            gxi = np.concatenate((gxi,[endx]))
        if len(gxi) > ncols:
            gxi = gxi[0:-1]
        if len(gyi) < nrows:
            gyi = np.concatenate((gyi,[starty]))
        if len(gyi) > nrows:
            gyi = gyi[0:-1]

        xi = (gxi - ulx1)/xdim1
        yi = (uly1 - gyi)/ydim1

        return (xi,yi)
    
    def interpolateToGrid(self,geodict,method='linear'):
        """
        Given a geodict specifying another grid extent and resolution, resample current grid to match.
        
        :param geodict: 
            geodict dictionary from another grid whose extents are inside the extent of this grid.
        :keyword method: 
            Optional interpolation method - ['linear', 'cubic','quintic','nearest']
        :raises GridException: 
           If the Grid object upon which this function is being called is not completely contained by the grid to which this Grid is being resampled.
        :raises GridException: 
           If the resulting interpolated grid shape does not match input geodict.

        This function modifies the internal griddata and geodict object variables.
        """
        xi,yi = self._getInterpCoords(geodict)

        #now using scipy interpolate functions
        baserows,basecols = self._geodict['nrows'],self._geodict['ncols']
        basex = np.arange(0,basecols) #base grid PIXEL coordinates
        basey = np.arange(0,baserows)
        if method in ['linear','cubic','quintic']:
            if not np.isnan(self._data).any():
                #at the time of this writing, interp2d does not support NaN values at all.
                #switching to griddata, which is slower by ~2 orders of magnitude but supports NaN.
                f = interpolate.interp2d(basex,basey,self._data,kind=method)
                self._data = f(xi,yi)
            else:
                xi,yi = np.meshgrid(xi,yi)
                newrows,newcols = xi.shape
                xi = xi.flatten()
                yi = yi.flatten()
                xnew = np.zeros((len(xi),2))
                xnew[:,0] = xi
                xnew[:,1] = yi
                basex,basey = np.meshgrid(basex,basey)
                basex = basex.flatten()
                basey = basey.flatten()
                xold = np.zeros((len(basex),2))
                xold[:,0] = basex
                xold[:,1] = basey
                self._data = interpolate.griddata(xold,self._data.flatten(),xnew,method=method)
                self._data = self._data.reshape((newrows,newcols))
        else:
            x,y = np.meshgrid(basex,basey)
            f = interpolate.NearestNDInterpolator(zip(x.flatten(),y.flatten()),self._data.flatten())
            newrows = geodict['nrows']
            newcols = geodict['ncols']
            xi = np.tile(xi,(newrows,1))
            yi = np.tile(yi.reshape(newrows,1),(1,newcols))
            self._data = f(zip(xi.flatten(),yi.flatten()))
            self._data = self._data.reshape(xi.shape)
                                                  
            
        nrows,ncols = geodict['nrows'],geodict['ncols']
        dims = self._data.shape
        nrows_new = dims[0]
        ncols_new = dims[1]
        if nrows_new != nrows or ncols_new != ncols:
            msg = "Interpolation failed!  Results (%i,%i) don't match (%i,%i)!" % (nrows_new,ncols_new,nrows,ncols)
            raise GridException(msg)
        #now the extents and resolution of the two grids should be identical...
        self._geodict['nrows'] = geodict['nrows']
        self._geodict['ncols'] = geodict['ncols']
        self._geodict['xmin'] = geodict['xmin']
        self._geodict['xmax'] = geodict['xmax']
        self._geodict['ymin'] = geodict['ymin']
        self._geodict['ymax'] = geodict['ymax']
        self._geodict['xdim'] = geodict['xdim']
        self._geodict['ydim'] = geodict['ydim']

    
def _test():
    xmin = 118.5
    xmax = 120.5
    ymin = 32.0
    ymax = 34.0
    xdim = 0.25
    ydim = 0.25
    ncols = len(np.arange(xmin,xmax+xdim,xdim))
    nrows = len(np.arange(ymin,ymax+ydim,ydim))
    data = np.arange(0,nrows*ncols)
    data.shape = (nrows,ncols)
    geodict = {'xmin':xmin,
               'xmax':xmax,
               'ymin':ymin,
               'ymax':ymax,
               'xdim':xdim,
               'ydim':ydim,
               'nrows':nrows,
               'ncols':ncols}
    grid = Grid2D(data,geodict)
    lat,lon = grid.getLatLon(0,0)
    row,col = grid.getRowCol(lat,lon)
    value = grid.getValue(lat,lon)
    frow,fcol = grid.getRowCol(33.3,119.1,returnFloat=True)
    irow,icol = grid.getRowCol(33.3,119.1,returnFloat=False)
    trimbounds = (xmin+0.3,xmax-0.3,ymin+0.3,ymax-0.3)
    print grid.getData()
    print grid.getData().shape
    grid.trim(trimbounds,resample=True)
    print lat,lon
    print row,col
    print value
    print irow,icol
    print frow,fcol
    print grid.getData()
    print grid.getData().shape
    
if __name__ == '__main__':
    _test()
    
    
        
