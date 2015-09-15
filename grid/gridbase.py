#!/usr/bin/env python

import abc
import numpy as np

#third party imports
from dataset import DataSet

class Grid(DataSet):
    """
    An abstract class to represent lat/lon gridded datasets. Grids are
    assumed to be pixel-registered - that is, grid coordinates
    represent the value at the *center* of the cells.
    """
    @staticmethod
    def getLatLonMesh(geodict):
        lons = np.linspace(geodict['xmin'],geodict['xmax'],num=geodict['ncols'])
        lats = np.linspace(geodict['ymin'],geodict['ymax'],num=geodict['nrows'])
        lon,lat = np.meshgrid(lons,lats)
        return (lat,lon)
    
    @abc.abstractmethod
    def getGeoDict(self):
        """
        Return a reference to the geodict inside the Grid
        
        :returns:
          A reference to a dictionary (see constructor).
        """
        raise NotImplementedError('getGeoDict method not implemented in base class')

    @abc.abstractmethod
    def getLatLon(self,row,col):
        """Return geographic coordinates (lat/lon decimal degrees) for given data row and column.
        
        :param row: 
           Row dimension index into internal data array.
        :param col: 
           Column dimension index into internal data array.
        :returns: 
           Tuple of latitude and longitude.
        """
        raise NotImplementedError('getLatLon method not implemented in base class')

    @abc.abstractmethod
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
        raise NotImplementedError('getRowCol method not implemented in base class')



    
    
    
        
