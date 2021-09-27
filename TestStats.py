from avgdist import DistCalc

import fiona
from shapely.geometry import Point, Polygon, shape
from numpy import sqrt, std
from math import fsum
import timeit

from os import path


def TestCompareGEOS():
    """[summary]
    Compare Statistic performance with typical distance computation metrics in GEOS
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

            # minimum distance
            distmin = pnt.distance(plg)
            # maximum distance/ hausdorff distance
            distmax = pnt.hausdorff_distance(plg)
            # center distance
            distcent = pnt.distance(plg.centroid)
            # mabb distance
            distmabbcf, t2 = DistCalc.DistCalcMABBCF(
                pnt, plg)
            # average distance
            distcf, t2 = DistCalc.DistCalcCF(pnt, plg)
            # Rndpnt n-max distance
            distnmax, t2 = DistCalc.DistCalcRPT(pnt, plg, method='nmax')
            # Rndpnt n-hat distance
            distnhat, t2 = DistCalc.DistCalcRPT(pnt, plg, method='nhat')
            print((f['properties']['CID'], distcf, distmabbcf, distnmax, distnhat, distmin,
                   distmax, distcent))


if __name__ == "__main__":
    # Run Example 2
    TestCompareGEOS()
