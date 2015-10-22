#!/usr/bin/env python

#stdlib imports
import sys
import collections
import datetime
import time
import os.path

#third party imports
import h5py
from scipy.io import netcdf 
import numpy as np
from multiple import MultiGrid
from shake import ShakeGrid
from grid2d import Grid2D
from dataset import DataSetException

class MultiHazardGrid(MultiGrid):
    def __init__(self,layers,geodict,origin,header,metadata=None):
        """Construct a ShakeGrid object.
        :param layers:
           OrderedDict containing ShakeMap data layers (keys are 'pga', etc., values are 2D arrays of data).
        :param geodict:
           Dictionary specifying the spatial extent,resolution and shape of the data.
        :param origin:
          Dictionary with elements:
            - id String of event ID (i.e., 'us2015abcd')
            - source String containing originating network ('us')
            - time Float event magnitude
            - lat Float event latitude
            - lon Float event longitude
            - depth Float event depth
            - magnitude Datetime object representing event origin time.
        :param header:
          Dictionary with elements:
            - type Type of multi-layer earthquake induced hazard ('shakemap','gfe')
            - version Integer product version (1)
            - process_time Python datetime indicating when data was created.
            - code_version String version of code that created this file (i.e.,'4.0')
            - originator String representing network that created the hazard grid.
            - product_id String representing the ID of the product (may be different from origin ID)
            - map_status String, one of RELEASED, ??
            - event_type String, one of ['ACTUAL','SCENARIO']
        :param metadata:
          Dictionary of dictionaries containing other metadata users wish to preserve.
        :returns:
           A MultiHazardGrid object.
        """
        self._layers = collections.OrderedDict()
        self._geodict = geodict
        for layerkey,layerdata in layers.iteritems():
            try:
                self._layers[layerkey] = Grid2D(data=layerdata,geodict=geodict)
            except:
                pass
        self._setHeader(header)
        self._setOrigin(origin)
        self._metadata = metadata.copy()

    def _saveDict(self,group,mydict):
        ALLOWED = [str,unicode,int,float,
                   long,list,tuple,np.ndarray,
                   dict,datetime.datetime,
                   collections.OrderedDict]
        for key,value in mydict.iteritems():
            tvalue = type(value)
            if tvalue not in ALLOWED:
                raise DataSetException('Unsupported metadata value type "%s"' % tvalue)
            if not isinstance(value,dict):
                if isinstance(value,datetime.datetime):
                    value = time.mktime(value.timetuple())
                group.attrs[key] = value
            else:
                subgroup = group.create_group(key)
                self._saveDict(subgroup,value)

    @classmethod
    def _loadDict(cls,group):
        tdict = {}
        for key,value in group.attrs.iteritems(): #attrs are NOT subgroups
            if key.find('time') > -1:
                value = value = datetime.datetime.utcfromtimestamp(value)
            tdict[key] = value
        for key,value in group.iteritems(): #these are going to be the subgroups
            tdict[key] = cls._loadDict(value)
        return tdict
                
            
    def save(self,filename,format='hdf'):
        f = h5py.File(filename, "w")
        #create two top-level groups that should always be present
        header = f.create_group('header')
        self._saveDict(header,self._header)

        origin = f.create_group('origin')
        self._saveDict(origin,self._origin)

        #write out any other metadata, creating groups recursively as needed
        metadata = f.create_group('metadata')
        self._saveDict(metadata,self._metadata)

        xvar = np.linspace(self._geodict['xmin'],self._geodict['xmax'],self._geodict['ncols'])
        yvar = np.linspace(self._geodict['ymin'],self._geodict['ymax'],self._geodict['nrows'])
        x = f.create_dataset('x',data=xvar,compression='gzip')
        y = f.create_dataset('y',data=yvar,compression='gzip')
        for layerkey,layer in self._layers.iteritems():
            dset = f.create_dataset(layerkey,data=layer.getData(),compression='gzip')

        f.close()

    @classmethod
    def load(cls,filename):
        f = h5py.File(filename, "r")
        REQUIRED_GROUPS = ['origin','header']
        REQUIRED_DATASETS = ['x','y']
        for group in REQUIRED_GROUPS:
            if group not in f.keys():
                raise DataSetException('Missing required metadata group "%s"' % group)
        for dset in REQUIRED_DATASETS:
            if dset not in f.keys():
                raise DataSetException('Missing required data set "%s"' % dset)

        header = {}
        for key,value in f['header'].attrs.iteritems():
            if key.find('time') > -1:
                value = datetime.datetime.utcfromtimestamp(value)
            header[key] = value

        origin = {}
        for key,value in f['origin'].attrs.iteritems():
            if key.find('time') > -1:
                value = datetime.datetime.utcfromtimestamp(value)
            origin[key] = value

        if 'metadata' in f.keys():
            metadata = cls._loadDict(f['metadata'])

        geodict = {}
        xvar = f['x'][:]
        yvar = f['y'][:]
        geodict['xmin'] = xvar[0]
        geodict['xmax'] = xvar[-1]
        geodict['ymin'] = yvar[0]
        geodict['ymax'] = yvar[-1]
        geodict['nrows'] = len(yvar)
        geodict['ncols'] = len(xvar)
        geodict['xdim'] = xvar[1]-xvar[0]
        geodict['ydim'] = yvar[1]-yvar[0]
        layers = collections.OrderedDict()
        dictDict = {}
        for key in f.keys():
            keytype = str(type(f[key]))
            if keytype.find('Dataset') > -1:
                if key in REQUIRED_DATASETS:
                    continue
                layers[key] = f[key][:]

        f.close()
        cls(layers,geodict,origin,header,metadata=metadata)
        

    def _validateDict(self,tdict):
        ALLOWED = ['str','unicode','int','float','long','list','tuple','numpy.ndarray']
        #input dict can only have strings, numbers, lists, tuples, or numpy arrays as values (no sub-dictionaries)
        for key,value in tdict.iteritems():
            tvalue = type(value)
            if tvalue not in ALLOWED:
                raise DataSetException('Data type "%s" not allowed in MultiHazardGrid extra metadata' % tvalue)

    def _setHeader(self,header):
        self._header = header.copy() #validate later

    def _setOrigin(self,origin):
        self._origin = origin.copy() #validate later

if __name__ == '__main__':
    shakefile = sys.argv[1]
    t1 = datetime.datetime.now()
    sgrid = ShakeGrid.load(shakefile)
    t2 = datetime.datetime.now()
    origin = {}
    origin['id'] = sgrid._eventDict['event_id']
    origin['source'] = sgrid._eventDict['event_network']
    origin['time'] = sgrid._eventDict['event_timestamp']
    origin['lat'] = sgrid._eventDict['lat']
    origin['lon'] = sgrid._eventDict['lon']
    origin['depth'] = sgrid._eventDict['depth']
    origin['magnitude'] = sgrid._eventDict['magnitude']

    header = {}
    header['type'] = 'shakemap'
    header['version'] = sgrid._shakeDict['shakemap_version']
    header['process_time'] = sgrid._shakeDict['process_timestamp']
    header['code_version'] = sgrid._shakeDict['code_version']
    header['originator'] = sgrid._shakeDict['shakemap_originator']
    header['product_id'] = sgrid._shakeDict['shakemap_id']
    header['map_status'] = sgrid._shakeDict['map_status']
    header['event_type'] = sgrid._shakeDict['shakemap_event_type']

    layers = collections.OrderedDict()
    for layername,layerdata in sgrid.getData().iteritems():
        layers[layername] = layerdata.getData()

    tdict = {'name':'fred','family':{'wife':'wilma','daughter':'pebbles'}}
    mgrid = MultiHazardGrid(layers,sgrid.getGeoDict(),origin,header,metadata={'flintstones':tdict})
    mgrid.save('test.hdf')
    t3 = datetime.datetime.now()
    mgrid2 = MultiHazardGrid.load('test.hdf')
    t4 = datetime.datetime.now()
    xmlmb = os.path.getsize(shakefile)/float(1e6)
    hdfmb = os.path.getsize('test.hdf')/float(1e6)
    xmltime = (t2-t1).seconds + (t2-t1).microseconds/float(1e6)
    hdftime = (t4-t3).seconds + (t4-t3).microseconds/float(1e6)
    print 'Input XML file size: %.1f MB (loading time %.1f seconds)' % (xmlmb,xmltime)
    print 'Output HDF file size: %.1f MB (loading time %.1f seconds)' % (hdfmb,hdftime)
    os.remove('test.hdf')    
