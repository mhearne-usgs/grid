#!/usr/bin/env python

#stdlib imports
from xml.dom import minidom
from datetime import datetime
from collections import OrderedDict
import re
import sys
import StringIO

#third party
from gridbase import Grid
from multiple import MultiGrid
from dataset import DataSetException
from grid2d import Grid2D
import numpy as np

GRIDKEYS = {'event_id':'string',
            'shakemap_id':'string',
            'shakemap_version':'int',
            'code_version':'string',
            'process_timestamp':'datetime',
            'shakemap_originator':'string',
            'map_status':'string',
            'shakemap_event_type':'string'}

EVENTKEYS = {'event_id':'string',
             'magnitude':'float',
             'depth':'float',
             'lat':'float',
             'lon':'float',
             'event_timestamp':'datetime',
             'event_network':'string',
             'event_description':'string'}

SPECKEYS = {'lon_min':'float',
            'lon_max':'float',
            'lat_min':'float',
            'lat_max':'float',
            'nominal_lon_spacing':'float',
            'nominal_lat_spacing':'float',
            'nlon':'int',
            'nlat':'int'}

FIELDKEYS = OrderedDict()
FIELDKEYS['pga'] = ('pctg','%.2f')
FIELDKEYS['pgv'] = ('cms','%.2f')
FIELDKEYS['mmi'] = ('intensity','%.2f')
FIELDKEYS['psa03'] = ('pctg','%.2f')
FIELDKEYS['psa10'] = ('pctg','%.2f')
FIELDKEYS['psa30'] = ('pctg','%.2f')
FIELDKEYS['stdpga'] = ('ln(pctg)','%.2f')
FIELDKEYS['urat'] = ('','%.2f')
FIELDKEYS['svel'] = ('ms','%.2f')
          

TIMEFMT = '%Y-%m-%dT%H:%M:%S'

def readElement(element,keys):
    eldict = OrderedDict()
    for key,dtype in keys.iteritems():
        if dtype == 'datetime':
            eldict[key] = datetime.strptime(element.getAttribute(key)[0:19],TIMEFMT)
        elif dtype == 'int':
            eldict[key] = int(element.getAttribute(key))
        elif dtype == 'float':
            eldict[key] = float(element.getAttribute(key))
        else:
            eldict[key] = element.getAttribute(key)
    return eldict

def getXMLText(fileobj):
    tline = fileobj.readline() 
    datamatch = re.compile('grid_data')
    xmltext = ''
    tlineold = ''
    while not datamatch.search(tline) and tline != tlineold:
        tlineold = tline
        xmltext = xmltext+tline
        tline = fileobj.readline()    

    xmltext = xmltext+'</shakemap_grid>'
    return xmltext

def getHeaderData(fileobj):
    xmltext = getXMLText(fileobj)
    root = minidom.parseString(xmltext)
    griddict = OrderedDict()
    gridel = root.getElementsByTagName('shakemap_grid')[0]
    griddict = readElement(gridel,GRIDKEYS)
    eventel = root.getElementsByTagName('event')[0]
    eventdict = readElement(eventel,EVENTKEYS)
    specel = root.getElementsByTagName('grid_specification')[0]
    specdict = readElement(specel,SPECKEYS)
    field_elements = root.getElementsByTagName('grid_field')
    fields = []
    for fieldel in field_elements:
        att = fieldel.getAttribute('name').lower()
        if att in ['lon','lat']:
            continue
        fields.append(att)

    uncertainties = OrderedDict()
    unc_elements = root.getElementsByTagName('event_specific_uncertainty')
    for uncel in unc_elements:
        key = uncel.getAttribute('name')
        value = float(uncel.getAttribute('value'))
        numsta = int(uncel.getAttribute('numsta'))
        uncertainties[key] = (value,numsta)

    return (griddict,eventdict,specdict,fields,uncertainties)

def readShakeFile(fileobj):
    griddict,eventdict,specdict,fields,uncertainties = getHeaderData(fileobj)
    ncols = specdict['nlon']
    nrows = specdict['nlat']
    layers = OrderedDict()

    #use the numpy loadtxt function to read in the actual data
    #we're cheating by telling numpy.loadtxt that the last two lines of the XML file are comments
    data = np.loadtxt(fileobj,comments='<') 
    data = data[:,2:] #throw away lat/lon columns
    for i in range(0,len(fields)):
        field = fields[i]
        layers[field] = data[:,i].reshape(nrows,ncols)

    #create the geodict from the grid_spec element
    geodict = {'xmin':specdict['lon_min'],
               'xmax':specdict['lon_max'],
               'ymin':specdict['lat_min'],
               'ymax':specdict['lat_max'],
               'xdim':specdict['nominal_lon_spacing'],
               'ydim':specdict['nominal_lat_spacing'],
               'nrows':specdict['nlat'],
               'ncols':specdict['nlon']}
    
    return (layers,geodict,eventdict,griddict,uncertainties)

class ShakeGrid(MultiGrid):
    def __init__(self,layers,geodict,eventDict,shakeDict,uncertaintyDict):
        self._layers = OrderedDict()
        self._geodict = geodict
        for layerkey,layerdata in layers.iteritems():
            self._layers[layerkey] = Grid2D(data=layerdata,geodict=geodict)
        self.setEventDict(eventDict)
        self.setShakeDict(shakeDict)
        self.setUncertaintyDict(uncertaintyDict)

    @classmethod
    def getFileGeoDict(cls,shakefile):
        isFileObj = False
        if not hasattr(shakefilename,'read'):
            shakefile = open(shakefilename,'r')
        else:
            isFileObj = True
            shakefile = shakefilename
        griddict,eventdict,specdict,fields,uncertainties = getHeaderData(f)
        if isFileObj:
            shakefile.close()
        geodict = {'xmin':specdict['lon_min'],
                   'xmax':specdict['lon_max'],
                   'ymin':specdict['lat_min'],
                   'ymax':specdict['lat_max'],
                   'xdim':specdict['nominal_lon_spacing'],
                   'ydim':specdict['nominal_lat_spacing'],
                   'nrows':specdict['nlat'],
                   'ncols':specdict['nlon']}
        return geodict

    @classmethod
    def load(cls,shakefilename,samplegeodict=None,resample=False,method='linear',doPadding=False,padValue=np.nan):
        #geodict can have xdim/ydim OR ncols/nrows.  If given both, xdim/ydim will be used to re-calculate nrows/ncols
        isFileObj = False
        if not hasattr(shakefilename,'read'):
            shakefile = open(shakefilename,'r')
        else:
            isFileObj = True
            shakefile = shakefilename

        #fill in nrows/ncols or xdim/ydim, whichever is not specified.  xdim/ydim dictate if both pairs are specified.
        samplegeodict = cls.fillGeoDict(samplegeodict)
        layers,geodict,eventDict,shakeDict,uncertaintyDict = readShakeFile(shakefile)
            
        if not isFileObj:
            shakefile.close()

        if samplegeodict is None:
            pass #everything we need has already been retrieved
        else:
            if doPadding:
                leftpad,rightpad,bottompad,toppad,geodict = _getPadding(geodict,bounds,padValue)
                for layername,layerdata in layers.iteritems():
                    #pad left side
                    layerdata = np.hstack((leftpad,layerdata))
                    #pad right side
                    layerdata = np.hstack((layerdata,rightpad))
                    #pad bottom
                    layerdata = np.vstack((layerdata,bottompad))
                    #pad top
                    layerdata = np.vstack((toppad,layerdata))
                    grid = Grid2D(layerdata,geodict)
                    if resample: #should I just do an interpolateToGrid() here?
                        grid.trim(bounds,resample=resample,method=method)
                    layers[layername] = grid.getData()
            else:
                for layername,layerdata in layers.iteritems():
                    grid = Grid2D(layerdata,geodict)
                    #should I do an interpolateToGrid() here too?
                    grid.trim(bounds,resample=resample,method=method)
                    layers[layername] = grid.getData()
                geodict = grid.getGeoDict().copy()
            
        return cls(layers,geodict,eventDict,shakeDict,uncertaintyDict)
    
    def save(self,filename,version=1):
        isFile = False
        if not hasattr(filename,'read'):
            isFile = True
            f = open(filename,'wt')
        else:
            f = filename    
        SCHEMA1 = 'http://www.w3.org/2001/XMLSchema-instance'
        SCHEMA2 = 'http://earthquake.usgs.gov/eqcenter/shakemap'
        SCHEMA3 = 'http://earthquake.usgs.gov http://earthquake.usgs.gov/eqcenter/shakemap/xml/schemas/shakemap.xsd'

        f.write('<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>')
        fmt = '<shakemap_grid xmlns:xsi="%s" xmlns="%s" xsi:schemaLocation="%s" event_id="%s" shakemap_id="%s" shakemap_version="%i" code_version="%s" process_timestamp="%s" shakemap_originator="%s" map_status="%s" shakemap_event_type="%s">\n'
        tpl = (SCHEMA1,SCHEMA2,SCHEMA3,
               self._shakeDict['event_id'],self._shakeDict['shakemap_id'],self._shakeDict['shakemap_version'],
               self._shakeDict['code_version'],datetime.utcnow().strftime(TIMEFMT),
               self._shakeDict['shakemap_originator'],self._shakeDict['map_status'],self._shakeDict['shakemap_event_type'])
        f.write(fmt % tpl)
        fmt = '<event event_id="%s" magnitude="%.1f" depth="%.1f" lat="%.4f" lon="%.4f" event_timestamp="%s" event_network="%s" event_description="%s"/>\n'
        tpl = (self._eventDict['event_id'],self._eventDict['magnitude'],self._eventDict['depth'],
               self._eventDict['lat'],self._eventDict['lon'],self._eventDict['event_timestamp'].strftime(TIMEFMT),
               self._eventDict['event_network'],self._eventDict['event_description'])
        f.write(fmt % tpl)
        fmt = '<grid_specification lon_min="%.4f" lat_min="%.4f" lon_max="%.4f" lat_max="%.4f" nominal_lon_spacing="%.4f" nominal_lat_spacing="%.4f" nlon="%i" nlat="%i"/>'
        tpl = (self._geodict['xmin'],self._geodict['ymin'],self._geodict['xmax'],self._geodict['ymax'],
               self._geodict['xdim'],self._geodict['ydim'],self._geodict['ncols'],self._geodict['nrows'])
        f.write(fmt % tpl)
        fmt = '<event_specific_uncertainty name="%s" value="%.4f" numsta="%i" />\n'
        for key,unctuple in self._uncertaintyDict.iteritems():
            value,numsta = unctuple
            tpl = (key,value,numsta)
            f.write(fmt % tpl)
        f.write('<grid_field index="1" name="LON" units="dd" />\n')
        f.write('<grid_field index="2" name="LAT" units="dd" />\n')
        idx = 3
        fmt = '<grid_field index="%i" name="%s" units="%s" />\n'
        data_formats = ['%.4f','%.4f']
        for field in self._layers.keys():
            tpl = (idx,field.upper(),FIELDKEYS[field][0])
            data_formats.append(FIELDKEYS[field][1])
            f.write(fmt % tpl)
            idx += 1
        f.write('<grid_data>\n')
        lat,lon = Grid().getLatLonMesh(self._geodict)
        nfields = 2 + len(self._layers)
        data = np.zeros((self._geodict['nrows']*self._geodict['ncols'],nfields))
        data[:,0] = lat.flatten()
        data[:,1] = lon.flatten()
        fidx = 2
        for grid in self._layers.values():
            data[:,fidx] = grid.getData().flatten()
            fidx += 1
        np.savetxt(f,data,delimiter=' ',fmt=data_formats)
        f.write('</grid_data>\n')
        if isFile:
            f.close()

    def _checkType(self,key,dtype):
        if dtype == 'string' and (not isinstance(key,str) and not isinstance(key,unicode)):
            return False
        if dtype == 'int' and not isinstance(key,int):
            return False
        if dtype == 'float' and not isinstance(key,float):
            return False
        if dtype == 'datetime' and not isinstance(key,datetime):
            return False
        return True
    
    def setEventDict(self,eventdict):
        for key,dtype in EVENTKEYS.iteritems():
            if not eventdict.has_key(key):
                raise DataSetException('eventdict is missing key "%s"' % key)
            if not self._checkType(eventdict[key],dtype):
                raise DataSetException('eventdict key value "%s" is the wrong datatype' % str(eventdict[key]))
        self._eventDict = eventdict.copy()

    def setShakeDict(self,shakedict):
        for key,dtype in GRIDKEYS.iteritems():
            if not shakedict.has_key(key):
                raise DataSetException('shakedict is missing key "%s"' % key)
            if not self._checkType(shakedict[key],dtype):
                raise DataSetException('shakedict key value "%s" is the wrong datatype' % str(shakedict[key]))
        self._shakeDict = shakedict.copy()

    def setUncertaintyDict(self,uncertaintyDict):
        self._uncertaintyDict = uncertaintyDict.copy()
        
def _trim_test(shakefile):
    geodict = getShakeDict(shakefile)
    #bring in the shakemap by a half dimension (quarter on each side)
    lonrange = geodict['xmax'] - geodict['xmin']
    latrange = geodict['ymax'] - geodict['ymin']
    newxmin = geodict['xmin'] + lonrange/4.0
    newxmax = geodict['xmax'] - lonrange/4.0
    newymin = geodict['ymin'] + latrange/4.0
    newymax = geodict['ymax'] - latrange/4.0
    newbounds = (newxmin,newxmax,newymin,newymax)
    grid = ShakeGrid.load(shakefile)
    grid.trim(newbounds)

def _save_test():
    pga = np.arange(0,16,dtype=np.float32).reshape(4,4)
    pgv = np.arange(1,17,dtype=np.float32).reshape(4,4)
    mmi = np.arange(2,18,dtype=np.float32).reshape(4,4)
    geodict = {'xmin':0.5,'ymax':3.5,'ymin':0.5,'xmax':3.5,'xdim':1.0,'ydim':1.0}
    layers = OrderedDict()
    layers['pga'] = pga
    layers['pgv'] = pgv
    layers['mmi'] = mmi
    shakeDict = {'event_id':'usabcd1234',
                 'shakemap_id':'usabcd1234',
                 'shakemap_version':1,
                 'code_version':'4.0',
                 'process_timestamp':datetime.utcnow(),
                 'shakemap_originator':'us',
                 'map_status':'RELEASED',
                 'shakemap_event_type':'ACTUAL'}
    eventDict = {'event_id':'usabcd1234',
                 'magnitude':7.6,
                 'depth':1.4,
                 'lat':2.0,
                 'lon':2.0,
                 'event_timestamp':datetime.utcnow(),
                 'event_network':'us',
                 'event_description':'sample event'}
    uncDict = {'pga':(0.0,0),
               'pgv':(0.0,0),
               'mmi':(0.0,0)}
    shake = ShakeGrid(layers,geodict,eventDict,shakeDict,uncDict)
    fobj = StringIO.StringIO()
    shake.save(fobj,version=3)
    filestr = fobj.getvalue()
    fobj.seek(0)
    shake2 = ShakeGrid.load(fobj)
    fobj.seek(0)
    print shake2.getLayer('pga').getData()
    bigbounds = (-0.5,4.5,-0.5,4.5)
    shake3 = ShakeGrid.load(fobj,bounds=bigbounds,resample=False,doPadding=False,padValue=np.nan)
    print pga
    print shake3.getLayer('pga').getData()
    fobj.seek(0)
    shake4 = ShakeGrid.load(fobj,bounds=bigbounds,resample=False,doPadding=True,padValue=np.nan)
    print shake4.getLayer('pga').getData()
    print shake4.getGeoDict()
    fobj.seek(0)
    littlebounds = (1.0,3.0,1.0,3.0)
    shake5 = ShakeGrid.load(fobj,bounds=littlebounds,resample=True,doPadding=False,padValue=np.nan)
    print shake5.getLayer('pga').getData()
    print shake5.getGeoDict()
    fobj.seek(0)
    bigbounds2 = (-1.0,5.0,-1.0,5.0)
    shake6 = ShakeGrid.load(fobj,bounds=bigbounds2,resample=True,doPadding=True,padValue=np.nan)
    print shake6.getLayer('pga').getData()
    print shake6.getGeoDict()
                 
                 
if __name__ == '__main__':
    # shakefile = sys.argv[1]
    # _trim_test(shakefile)
    _save_test()
    
