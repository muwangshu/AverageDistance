# Examples
from avgdist import DistCalc
import fiona
from shapely.geometry import Point, Polygon, shape
from math import sqrt
from os import path


def example1():
    """[summary]
    A simple example calculating distance between a polygon and a point
    """
    # Define a polygon (Hexgon centers at (0,0), edge length=1)
    a = sqrt(3) / 2
    plg = Polygon([(-0.5, -a), (-1, 0), (-0.5, a),
                   (0.5, a), (1, 0), (0.5, -a)])
    # Define a point
    pnt = Point(0, 0)

    # Compute the distance between plg and pnt using the closed-form algorithm
    dist, time = DistCalc.DistCalcCF(pnt, plg)
    print("Distance: %s; Computation time: %ss" % (dist, time))

    # Other available methods
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


def example2():
    """[summary]
    An example calculating distance between a polygon and a point
    by reading polygon data from shapefile.
    """

    # data source:
    shp = path.join('TestData', 'ResearchArea.shp')

    # compute distance
    with fiona.open(shp) as source:
        for f in source:
            # polygon
            plg = shape(f['geometry'])
            # coordinate of the refernce point is stored in (REFX,REFY)
            pnt = Point(f['properties']['REFX'], f['properties']['REFY'])
            # print CID, distance, and computation time
            # replace DistCalc.DistCalcCF with other methods
            print((f['properties']['CID'], DistCalc.DistCalcCF(pnt, plg)))


if __name__ == "__main__":
    # Run Example 1
    example1()
    # Run Example 2
    example2()
