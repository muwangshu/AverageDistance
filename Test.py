from avgdist import DistCalc

import fiona
from shapely.geometry import Point, Polygon, shape
from numpy import sqrt, std
from math import fsum

from os import path


def EfficiencyTest():
    """[summary]
        Test the efficiency of each algorithm.
    """
    shp = path.join('TestData', 'ResearchArea.shp')
    with fiona.open(shp) as source:
        for f in source:
            if f['properties']['CID'] > 0:
                plg = shape(f['geometry'])
                pnt = Point(f['properties']['REFX'], f['properties']['REFY'])
                print((f['properties']['CID'],
                       DistCalc.DistCalcRPT(pnt, plg, method='nmax'), DistCalc.DistCalcRPT(pnt, plg, method='nhat'), DistCalc.DistCalcPART(pnt, plg), DistCalc.DistCalcCF(pnt, plg), DistCalc.DistCalcMABBCF(pnt, plg)))


def getAttributes():
    """[summary]
        Get the number of vertexes and the area of each polygon
    """
    shp = path.join('TestData', 'ResearchArea.shp')
    with fiona.open(shp) as source:
        for f in source:
            plg = shape(f['geometry'])
            print((f['properties']['CID'], len(
                plg.exterior.coords) - 1, plg.area))


def getTestPolygons():
    """[summary]
        Generate polygons in Test 1
    Returns:
        plg: Polygons in Test 1
    """
    plgs = [[] for i in range(4)]
    a = sqrt(3) / 2
    plgs[0] = Polygon([(-a, -0.5), (0, 1), (a, -0.5)])
    b = sqrt(2) / 2
    plgs[1] = Polygon([(-b, -b), (-b, b), (b, b), (b, -b)])
    plgs[2] = Polygon([(-0.5, -a), (-1, 0), (-0.5, a),
                       (0.5, a), (1, 0), (0.5, -a)])
    y1 = sqrt(4 - 0.05**2)
    plgs[3] = Polygon([(-0.05, y1), (0.05, y1), (0.05, y1 - 0.01), (0, y1 - 0.01), (0, 1.01), (0.05, 1.01),
                       (0.05, 1), (-0.05, 1), (-0.05, 1.01), (0, 1.01), (0, y1 - 0.01), (-0.05, y1 - 0.01)])
    return plgs


def AccuracyTestNMax():
    """[summary]
    Test the accuracy using n-max random points
    """
    plgs = getTestPolygons()
    pnt = Point(0, 0)
    # Calculate D using polygon partitioning method
    print("---------D---------------")
    for i in range(len(plgs)):
        print(DistCalc.DistCalcPART(pnt, plgs[i]))
    print("----------------------------")

    # calculate d-bar and MOE of d using nmax random points
    print("---------MOE=0.1 nmax=68------------")
    for i in range(len(plgs)):
        dist = [0 for k in range(100)]
        for k in range(100):
            dist[k], t = DistCalc.DistCalcRPT(pnt, plgs[i], 68)
        print(fsum(dist) / 100, std(dist) * 1.6449)
    print("----------------------------")

    print("---------MOE=0.05 nmax=271------------")
    for i in range(len(plgs)):
        dist = [0 for k in range(100)]
        for k in range(100):
            dist[k], t = DistCalc.DistCalcRPT(pnt, plgs[i], 271)
        print(fsum(dist) / 100, std(dist) * 1.6449)
    print("----------------------------")

    print("---------MOE=0.01 nmax=6764------------")
    for i in range(len(plgs)):
        dist = [0 for k in range(100)]
        for k in range(100):
            dist[k], t = DistCalc.DistCalcRPT(pnt, plgs[i], 6764)
        print(fsum(dist) / 100, std(dist) * 1.6449)
    print("----------------------------")
    return 0


def GetNHat():
    """[summary]
        Calculate the value of n-hat
    Returns:
        no return
    """
    plgs = getTestPolygons()
    pnt = Point(0, 0)
    # Calculate n-hat
    print("---------Calculate n-hat---------------")
    for i in range(len(plgs)):
        sd = DistCalc.getDistSD(pnt, plgs[i], cnt=50)
        nhat1 = 1.6499**2 * sd**2 / 0.1 / 0.1
        nhat2 = 1.6499**2 * sd**2 / 0.05 / 0.05
        nhat3 = 1.6499**2 * sd**2 / 0.01 / 0.01
        print((sd, nhat1, nhat2, nhat3))
    print("----------------------------")

    return 0


def AccuracyTestNHat():
    """
    [summary]
         Calculate d using n-hat random points
    """
    nhat = [[10, 39, 983],
            [11, 43, 1071],
            [13, 51, 1267],
            [66, 264, 6605]]
    result = []
    plgs = getTestPolygons()
    pnt = Point(0, 0)
    print("---------MOE=0.1, 0.05, 0.01 nhat------------")
    for i in range(len(plgs)):  # polygons
        row = []
        for j in range(3):  # MOEs
            dist = [0 for k in range(100)]
            for k in range(100):
                dist[k], t = DistCalc.DistCalcRPT(pnt, plgs[i], nhat[i][j])
            row.append((fsum(dist) / 100, std(dist) * 1.6449))
        result.append(row)
    print(result)
    print("----------------------------")


if __name__ == "__main__":
    # Tests performed in test.py:
    # Get polygons in Figure 4
    getTestPolygons()
    # compare the efficiency of different methods (Figure 6)
    EfficiencyTest()
    # Test the accuracy of the random points method using n-max points (Table 1)
    AccuracyTestNMax()
    # Calculate the value of n-hat (Table 2)
    GetNHat()
    # Test the accuracy of the random points method using n-hat points (Table 3)
    AccuracyTestNHat()
    pass
