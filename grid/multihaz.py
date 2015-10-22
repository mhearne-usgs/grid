#!/usr/bin/env python

#stdlib imports
import sys
from collections import OrderedDict
from datetime import datetime
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
    def __init__(self,layers,geodict,origin,header,dictdict=None):
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
        :param dictdict:
          Dictionary of dictionaries containing other metadata users wish to preserve.
        :returns:
           A MultiHazardGrid object.
        """
        self._layers = OrderedDict()
        self._geodict = geodict
        for layerkey,layerdata in layers.iteritems():
            try:
                self._layers[layerkey] = Grid2D(data=layerdata,geodict=geodict)
            except:
                pass
        self._setHeader(header)
        self._setOrigin(origin)
        self._dictDict = {}
        for tdictname,tdict in dictdict.iteritems():
            self._validateDict(tdict)
            self._dictDict[tdictname] = tdict.copy()

    def save(self,filename,format='hdf'):
        if format == 'hdf':
            f = h5py.File(filename, "w")
            origin = f.create_group('origin')
            for key,value in self._origin.iteritems():
                if isinstance(value,datetime):
                    value = time.mktime(value.timetuple())
                    origin.attrs[key] = value
                    
            header = f.create_group('header')
            for key,value in self._header.iteritems():
                if isinstance(value,datetime):
                    value = time.mktime(value.timetuple())
                header.attrs[key] = value
            for groupname,tdict in self._dictDict.iteritems():
                group = f.create_group(groupname)
                for key,value in tdict.iteritems():
                    group.attrs[key] = value
            xvar = np.linspace(self._geodict['xmin'],self._geodict['xmax'],self._geodict['ncols'])
            yvar = np.linspace(self._geodict['ymin'],self._geodict['ymax'],self._geodict['nrows'])
            x = f.create_dataset('x',data=xvar,compression='gzip')
            y = f.create_dataset('y',data=yvar,compression='gzip')
            for layerkey,layer in self._layers.iteritems():
                dset = f.create_dataset(layerkey,data=layer.getData(),compression='gzip')

            f.close()
        elif format == 'netcdf':
            cdf = netcdf.netcdf_file(filename,'w')
            for key,value in self._origin.iteritems():
                okey = 'origin_%s' % key
                if isinstance(value,datetime):
                    value = time.mktime(value.timetuple())
                setattr(cdf,okey,value)
            for key,value in self._header.iteritems():
                hkey = 'header_%s' % key
                if isinstance(value,datetime):
                    value = time.mktime(value.timetuple())
                setattr(cdf,hkey,value)
            cdf.createDimension('x', self._geodict['ncols'])
            cdf.createDimension('y', self._geodict['nrows'])
            xvar = f.createVariable('x', 'f', ('x',))
            yvar = f.createVariable('y', 'f', ('y',))
            xvar[:] = np.linspace(self._geodict['xmin'],self._geodict['xmax'],self._geodict['ncols'])
            yvar[:] = np.linspace(self._geodict['ymin'],self._geodict['ymax'],self._geodict['nrows'])
            for layerkey,layerdata in self._layers.iteritems():
                var = cdf.createVariable(layerkey,'f',('y','x',))
                var[:] = layerdata
            cdf.close()
        else:
            raise DataSetException('Unsupported file format "%s"' % format)
        

    @classmethod
    def getFileType(cls,grdfile):
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
                ftype = 'unknown'
        f.close()
        return ftype
            
    @classmethod
    def load(cls,filename):
        ftype = cls.getFileType(filename)
        if ftype == 'hdf':
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
                    value = datetime.utcfromtimestamp(value)
                header[key] = value

            origin = {}
            for key,value in f['origin'].attrs.iteritems():
                if key.find('time') > -1:
                    value = datetime.utcfromtimestamp(value)
                origin[key] = value

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
            layers = OrderedDict()
            dictDict = {}
            for key in f.keys():
                keytype = str(type(f[key]))
                if keytype.find('Dataset') > -1:
                    if key in REQUIRED_DATASETS:
                        continue
                    try:
                        layers[key] = f[key][:]
                    except:
                        pass
                else: #we have a group
                    if key in REQUIRED_GROUPS or key in REQUIRED_DATASETS:
                        continue
                    tdict = {}
                    for tkey,tvalue in f[key].attrs.iteritems():
                        tdict[key] = f[key].attrs[tkey]
                    dictDict[key] = tdict.copy()
            f.close()
            
        elif ftype == 'netcdf':
            HEADERKEYS = ['type','version','process_time','code_version','originator','product_id','map_status','event_type']
            ORIGINKEYS = ['id','source','time','lat','lon','depth','magnitude']
            cdf = netcdf.netcdf_file(filename)
            xvar = cdf.variables['x'].data.copy()
            yvar = cdf.variables['y'].data.copy()
            header = {}
            for key in HEADERKEYS:
                hkey = 'header_%s' % key
                if not hasattr(cdf,hkey):
                    raise DataSetException('Missing required attribute %s' % hkey)
                header[key] = getattr(cdf,hkey)
                if key.find('time') > -1:
                    header[key] = datetime.utcfromtimestamp(header[key])
            origin = {}
            for key in ORIGINKEYS:
                okey = 'origin_%s' % key
                if not hasattr(cdf,okey):
                    raise DataSetException('Missing required attribute %s' % okey)
                origin[key] = getattr(cdf,okey)
                if key.find('time') > -1:
                    header[key] = datetime.utcfromtimestamp(header[key])
            layers = OrderedDict()
            for key in cdf.variables.keys():
                if key in ['x','y']:
                    continue
                layers[key] = cdf.variables[key].copy()
            cdf.close()
        else:
            raise DataSetException('Unknown file type for file %s' % filename)
        cls(layers,geodict,origin,header,dictdict=dictDict)
        

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
    sgrid = ShakeGrid.load(shakefile)
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

    layers = OrderedDict()
    for layername,layerdata in sgrid.getData().iteritems():
        layers[layername] = layerdata.getData()
    
    mgrid = MultiHazardGrid(layers,sgrid.getGeoDict(),origin,header,dictdict={'event_specific_uncertainty':sgrid._uncertaintyDict})
    mgrid.save('test.hdf')
    mgrid2 = MultiHazardGrid.load('test.hdf')
    os.remove('test.hdf')    
