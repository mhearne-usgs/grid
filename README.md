Introduction
------------

grid is a project designed to provide a library of classes for dealing
with various grid formats, and performing some simple spatial
operations on the data.  The code is documented fairly well, and there
are IPython notebooks in the distribution.  They can be viewed here:

 * GDAL (ESRI format) grids : https://github.com/mhearne-usgs/grid/blob/master/notebooks/GDALGrid.ipynb
 * GMT format grids: https://github.com/mhearne-usgs/grid/blob/master/notebooks/GMTGrid.ipynb
 * ShakeMap format grids: https://github.com/mhearne-usgs/grid/blob/master/notebooks/ShakeMap.ipynb
 * The Grid2D class (superclass of GDAL and GMT grids): https://github.com/mhearne-usgs/grid/blob/master/notebooks/Grid2D.ipynb

Dependencies and Installation
-----------------------------

This library depends on:
 * numpy: <a href="http://www.numpy.org/">http://www.numpy.org/</a>
 * scipy: <a href="http://scipy.org/scipylib/index.html">http://scipy.org/scipylib/index.html</a>
 * h5py: <a href="http://www.h5py.org/">http://www.h5py.org/</a>
 * rasterio: <a href="https://github.com/mapbox/rasterio">https://github.com/mapbox/rasterio</a>
 
These packages are all either installed automatically by the Anaconda scientific Python distribution, or easily installed using the conda command.  The final dependency:

 * openquake: <a href="http://www.globalquakemodel.org/openquake/about/">http://www.globalquakemodel.org/openquake/about/</a>

can be installed by using pip with git:

pip install git+git://github.com/gem/oq-hazardlib.git

To install this package:

pip install git+git://github.com/mhearne-usgs/grid.git

Uninstalling and Updating
-------------------------

To uninstall:

pip uninstall grid

To update:

pip install -U git+git://github.com/mhearne-usgs/grid.git
