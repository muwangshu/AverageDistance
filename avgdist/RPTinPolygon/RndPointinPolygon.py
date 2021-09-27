#Fast algorithm to generate random points in a polygon.
import shapely
from ..EARCut import earcut
from numpy import array, sqrt, abs, where
from shapely.geometry import Polygon, Point
from random import uniform, random


def getRndPointinPolygon(plg, count):
    """[summary]

    Args:
        plg ([type]): [description]
        count ([type]): [description]

    Returns:
        [type]: [description]
    """
    coordlist = list(plg.exterior.coords)[:-1]
    coords = [0 for i in range(2 * len(coordlist))]
    for i in range(len(coordlist)):
        coords[2 * i], coords[2 * i + 1] = coordlist[i]
    triangles = earcut.earcut(coords)

    cumulativeDistribution = array(generateDistribution(triangles, coordlist))
    pnts = []
    for i in range(count):
        rnd = random()
        index = where(cumulativeDistribution > rnd)[0][0]
        t0 = triangles[3 * index]
        t1 = triangles[3 * index + 1]
        t2 = triangles[3 * index + 2]
        p = getRPTinTriangle(coordlist[t0], coordlist[t1], coordlist[t2])
        pnts.append(Point(p[0], p[1]))
    return pnts


def calTrArea(a, b, c):
    """[summary]

    Args:
        a ([type]): [description]
        b ([type]): [description]
        c ([type]): [description]

    Returns:
        [type]: [description]
    """
    return abs(a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])) / 2


def generateDistribution(triangles, coordlist):
    """[summary]

    Args:
        triangles ([type]): [description]
        coordlist ([type]): [description]

    Returns:
        [type]: [description]
    """
    trArea = []
    totalArea = 0
    for i in range(len(triangles) // 3):
        t0 = triangles[3 * i]
        t1 = triangles[3 * i + 1]
        t2 = triangles[3 * i + 2]
        trArea.append(calTrArea(coordlist[t0], coordlist[t1], coordlist[t2]))
        totalArea += trArea[i]

    cumulativeDistribution = []
    lastValue = 0
    for i in range(len(triangles) // 3):
        nextValue = lastValue + trArea[i] / totalArea
        cumulativeDistribution.append(nextValue)
        lastValue = nextValue
    return cumulativeDistribution


def selectRandomTriangle(cumulativeDistribution):
    """[summary]

    Args:
        cumulativeDistribution ([type]): [description]

    Returns:
        [type]: [description]
    """
    rnd = random()
    index = list(filter(lambda i: i > 0.6, cumulativeDistribution))[0]
    return triangles[index]


def getRPTinTriangle(a, b, c):
    """[summary]

    Args:
        a ([type]): [description]
        b ([type]): [description]
        c ([type]): [description]

    Returns:
        [type]: [description]
    """
    # make basis vectors from a->b and a->c
    u = (b[0] - a[0], b[1] - a[1])
    v = (c[0] - a[0], c[1] - a[1])

    # pick a random point in the unit square
    unit = (random(), random())

    # if the point is outside the triangle, remap it inside
    # by mirroring it along the diagonal
    if unit[0] + unit[1] > 1:
        unit = (1 - unit[1], 1 - unit[0])

    # now transform it to fit the basis
    return (unit[0] * u[0] + unit[1] * v[0] + a[0], unit[0] * u[1] + unit[1] * v[1] + a[1])


if __name__ == "__main__":

    y1 = sqrt(4 - 0.05**2)
    plg = Polygon([(-0.05, y1), (0.05, y1), (0.05, y1 - 0.01), (0, y1 - 0.01), (0, 1.01), (0.05, 1.01),
                   (0.05, 1), (-0.05, 1), (-0.05, 1.01), (0, 1.01), (0, y1 - 0.01), (-0.05, y1 - 0.01)])
    coordlist = list(plg.exterior.coords)[:-1]
    coords = [0 for i in range(2 * len(coordlist))]
    for i in range(len(coordlist)):
        coords[2 * i], coords[2 * i + 1] = coordlist[i]
    triangles = earcut.earcut(coords)

    cumulativeDistribution = array(generateDistribution(triangles, coordlist))
    for i in range(500):
        rnd = random()
        index = where(cumulativeDistribution > rnd)[0][0]
        t0 = triangles[3 * index]
        t1 = triangles[3 * index + 1]
        t2 = triangles[3 * index + 2]
        print(getRPTinTriangle(coordlist[t0], coordlist[t1], coordlist[t2]))
