#!/usr/bin/env python

#stdlib imports
from xml.dom import minidom
from datetime import datetime
from collections import OrderedDict
import re
import sys

#third party
from multiple import MultiGrid
from gridbase import GridException
from grid2d import Grid2D
import numpy as np

GRIDKEYS = {'event_id':'string',
            'shakemap_id':'string',
            'shakemap_version':'string',
            'code_version':'string',
            'process_timestamp':'datetime',
            'shakemap_originator':'string',
            'map_status':'string',
            'shakemap_event_type':'string'}

EVENTKEYS = {'magnitude':'float',
            'depth':'float',
            'lat':'float',
            'lon':'float',
            'event_timestamp':'datetime',
            'shakemap_originator':'string',
            'event_description':'string'}

SPECKEYS = {'lon_min':'float',
            'lon_max':'float',
            'lat_min':'float',
            'lat_max':'float',
            'nominal_lon_spacing':'float',
            'nominal_lat_spacing':'float',
            'nlon':'int',
            'nlat':'int'}

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
        uncertainties[key] = value

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
        layers[field] = np.flipud(data[:,i].reshape(nrows,ncols))

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
        for layerkey,layerdata in layers.iteritems():
            self._layers[layerkey] = Grid2D(data=layerdata,geodict=geodict)
        self.setEventDict(eventDict)
        self.setShakeDict(shakeDict)
        self.setUncertaintyDict(uncertaintyDict)
        
    @classmethod
    def load(cls,shakefilename,bounds=None,resample=False,method='linear',padValue=None,variable='MMI'):
        isFileObj = False
        if not hasattr(shakefilename,'read'):
            shakefile = open(shakefilename,'r')
        else:
            isFileObj = True
            shakefile = shakefilename


        layers,geodict,eventDict,shakeDict,uncertaintyDict = readShakeFile(shakefile)
            
        if not isFileObj:
            shakefile.close()

        #here we need to cut the data down if requested by the user, or expand it with padValue
        if bounds is not None:
            #the requested bounds
            xmin,xmax,ymin,ymax = bounds 
            #the actual bounds
            gxmin,gxmax,gymin,gymax = (geodict['xmin'],geodict['xmax'],geodict['ymin'],geodict['ymax'])
            xdim,ydim = (geodict['xdim'],geodict['ydim'])
            nrows,ncols = (geodict['nrows'],geodict['ncols'])
            if (xmin < gxmin or xmax > gxmax or ymin < gymin or ymax > gymax):
                if padValue is None:
                    raise GridException('Requesting bounds outside ShakeMap grid.')
                padleftcols = int((gxmin - xmin)/xdim)
                padrightcols = int((xmax - gxmax)/xdim)
                padbottomrows = int((gymin - ymin)/ydim)
                padtoprows = int((ymax - gymax)/ydim)
                if padleftcols < 0:
                    padleftcols = 0
                else:
                    geodict['xmin'] = xmin
                if padrightcols < 0:
                    padrightcols = 0
                else:
                    geodict['xmax'] = xmax
                if padbottomrows < 0:
                    padbottomrows = 0
                else:
                    geodict['ymin'] = ymin
                if padtoprows < 0:
                    padtoprows = 0
                else:
                    geodict['ymax'] = ymax
                leftpad = np.ones((nrows,padleftcols))
                rightpad = np.ones((nrows,padrightcols))
                ncols += padrightcols + padleftcols
                bottompad = np.ones((padbottomrows,ncols))
                toppad = np.ones((padtoprows,ncols))
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
                    gxmin,gxmax,gymin,gymax = (geodict['xmin'],geodict['xmax'],geodict['ymin'],geodict['ymax'])
                    if gxmin != xmin or gxmax != xmax or gymin != ymin or gymax != ymax:
                        grid.trim(bounds,resample=resample,method=method)
                    layers[layername] = grid.getData()
            else:
                for layername,layerdata in layers.iteritems():
                    grid = Grid2D(layerdata,geodict)
                    grid.trim(bounds,resample=resample,method=method)
                    layers[layername] = grid.getData()
                geodict = grid.getGeoDict().copy()
            
        return cls(layers,geodict,eventDict,shakeDict,uncertaintyDict)
    
    def save(filename):
        pass

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
                raise GridException('eventdict is missing key "%s"' % key)
            if not self._checkType(eventdict[key],dtype):
                raise GridException('eventdict key value "%s" is the wrong datatype' % str(eventdict[key]))
        self._eventDict = eventdict.copy()

    def setShakeDict(self,shakedict):
        for key,dtype in GRIDKEYS.iteritems():
            if not shakedict.has_key(key):
                raise GridException('shakedict is missing key "%s"' % key)
            if not self._checkType(shakedict[key],dtype):
                raise GridException('shakedict key value "%s" is the wrong datatype' % str(shakedict[key]))
        self._shakeDict = shakedict.copy()

    def setUncertaintyDict(self,uncertaintyDict):
        self._uncertaintyDict = uncertaintyDict.copy()
        
        
if __name__ == '__main__':
    shakefile = sys.argv[1]
    bounds = (-120.0,-118.0,33.0,35.0)
    grid = ShakeGrid.load(shakefile,bounds=bounds)
    grid = ShakeGrid.load(shakefile,bounds=bounds,resample=True)
    bounds = (-122.0,-116.0,32.0,37.0)
    grid = ShakeGrid.load(shakefile,bounds=bounds,resample=True,padValue=np.nan)
    bounds = (-122.0,-116.0,33.0,35.0)
    grid = ShakeGrid.load(shakefile,bounds=bounds,resample=True,padValue=np.nan)
