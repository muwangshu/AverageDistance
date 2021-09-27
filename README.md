# the avgdist package
Calculate Average Distance between Geometries


## Requirements
~~~
python 3.8+ with numpy and scipy
GDAL
shapely
fiona (for reading GIS data, optional)
~~~

## Tests in the Manuscript 
Test data stored in **TestData\\ResearchArea.shp**: Polygon Shapefile contaning test data. For each polygon, compute the distance between this polygon and the point defined in the using fields *REFX* and *REFY* as *X,Y* coordinates respectively.

Reproduce test result in the Test section using **Test.py**:
&nbsp;&nbsp;&nbsp;&nbsp;**getTestPolygons()**: Get polygons in Figure 4
&nbsp;&nbsp;&nbsp;&nbsp;**EfficiencyTest()**: compare the efficiency of different methods (Figure 6)
&nbsp;&nbsp;&nbsp;&nbsp;**AccuracyTestNMax()**: Test the accuracy of the random points method using n-max points (Table 1)
&nbsp;&nbsp;&nbsp;&nbsp;**GetNHat()**: Calculate the value of n-hat (Table 2)
&nbsp;&nbsp;&nbsp;&nbsp;**AccuracyTestNHat()**: Test the accuracy of the random points method using n-hat points (Table 3)

Reproduce the parallel computing test in the discussion section using **ParallelTest.py**:
&nbsp;&nbsp;&nbsp;&nbsp;**SingleComputeTest()**: Computing using single thread.
&nbsp;&nbsp;&nbsp;&nbsp;**ParallelComputeTest1()**: Parallel computing mode 1
&nbsp;&nbsp;&nbsp;&nbsp;**ParallelComputeTest2()**: Parallel computing mode 2

Reproduce the statistical performance test in the Appendix B using **TestStates.py** and **GenerateGRF**:
&nbsp;&nbsp;&nbsp;&nbsp;**TestCompareGEOS()**: Generate distances in the test
&nbsp;&nbsp;&nbsp;&nbsp;**GenerateGRF** is a R package: 
&nbsp;&nbsp;&nbsp;&nbsp;**SimulationUnif()**: Test when the individual distribution is homogenous 
&nbsp;&nbsp;&nbsp;&nbsp;**SimulationUnif()**: Test when the individual distribution is heterogenous 

## Examples
We provide some examples to show the useage of avgdist, available in **Examples.py**
###Example 1
A simple example calculating distance between a polygon and a point
1. Import the **avgdist** packages
```python
from avgdist import DistCalc
```
2. Import other packages
```python
from shapely.geometry import Point, Polygon, shape
from math import sqrt
```
3. Define a polygon and a point
```python
# Define a polygon (a hexgon centers at (0,0), edge length=1)
a = sqrt(3) / 2
plg = Polygon([(-0.5, -a), (-1, 0), (-0.5, a),
               (0.5, a), (1, 0), (0.5, -a)])
# Define a point
pnt = Point(0, 0)
```
4. Compute the distance between plg and pnt using the closed-form algorithm
```python
dist, time = DistCalc.DistCalcCF(pnt, plg)
print("Distance: %s; Computation time: %ss" % (dist, time))
```
5. Other available methods:
```python
# random points method, nmax random points
dist, time = DistCalc.DistCalcRPT(pnt, plg, method='nmax')
# random points method, nhat random points
dist, time = DistCalc.DistCalcRPT(pnt, plg, method='nhat')
# random points method, 1000 random points
dist, time = DistCalc.DistCalcRPT(pnt, plg, method=1000)
# numerical integration with polygon partitioning
dist, time = DistCalc.DistCalcPART(pnt, plg)
# closed from MABB-based method
dist, time = DistCalc.DistCalcMABBCF(pnt, plg)
```

###Example 2
An example calculating distance between a polygon and a point by reading polygon data from shapefile.

1. Import packages
```python
from avgdist import DistCalc
import fiona
from shapely.geometry import Point, Polygon, shape
from os import path
```

2. Location of the shapefile
```python
shp = path.join('TestData', 'ResearchArea.shp')
```

3. Read shapefile and compute distance
```python
with fiona.open(shp) as source:
    for f in source:
        # polygon
        plg = shape(f['geometry'])
        # coordinate of the refernce point is stored in (REFX,REFY)
        pnt = Point(f['properties']['REFX'], f['properties']['REFY'])
        # print CID, distance, and computation time
        # replace DistCalc.DistCalcCF with other methods
        print((f['properties']['CID'], DistCalc.DistCalcCF(pnt, plg)))
```







