
from shapely.geometry import Polygon, Point
from shapely import affinity

from numpy import sqrt, arctan2, log, abs, arctanh, sign, std
import numpy as np
from math import fsum
from scipy.integrate import dblquad, nquad
from random import uniform, random
from decimal import Decimal
import timeit
import functools

from .RPTinPolygon import RndPointinPolygon


def DistCalcAABB(pnt, plg):
    """[summary]
        Calculate distance using AABB method (numerical integration).
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        distance: distance
        error: error estimation
        time: computation time
    """
    # calculate distance using AABB
    t0 = timeit.default_timer()
    minu = plg.bounds[0]
    minv = plg.bounds[1]
    maxu = plg.bounds[2]
    maxv = plg.bounds[3]
    val = dblquad(EucDistWithBD, minv, maxv, lambda v: minu,
                  lambda v: maxu, args=(pnt.x, pnt.y, plg), epsrel=0.025)
    # val=nquad(EucDistWithBD,[[minu,maxu],[minv,maxv]],args=(pnt.x,pnt.y,plg))
    t1 = timeit.default_timer()

    return (val[0] / plg.area, val[1] / plg.area, t1 - t0)


def DistCalcMABB(pnt, plg):
    """[summary]
        Calculate distance using MABB method (numerical integration).
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        distance: distance
        error: error estimation
        time: computation time, rotation, after rotation, all
    """
    # calculate distance using MABB
    t0 = timeit.default_timer()

    mabb = plg.minimum_rotated_rectangle
    x0 = mabb.exterior.coords[0][0]
    y0 = mabb.exterior.coords[0][1]
    x1 = mabb.exterior.coords[1][0]
    y1 = mabb.exterior.coords[1][1]
    angle = arctan2(y1 - y0, x1 - x0)
    rmabb = affinity.rotate(mabb, -angle, origin=pnt, use_radians=True)
    minu = rmabb.bounds[0] - pnt.x
    minv = rmabb.bounds[1] - pnt.y
    maxu = rmabb.bounds[2] - pnt.x
    maxv = rmabb.bounds[3] - pnt.y
    # val = dblquad(EucDist, minv, maxv, lambda v: minu,
    #              lambda v: maxu, args=(pnt.x, pnt.y))
    t2 = timeit.default_timer()
    val = nquad(EucDistToZero, [[minu, maxu], [minv, maxv]])
    # print(val[1])
    t1 = timeit.default_timer()
    return (val[0] / rmabb.area, val[1], t1 - t0, t1 - t2, t2 - t0)


def DistCalcPART(pnt, plg):
    """[summary]
        Calculate distance using partition method (numerical integration).
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        distance: distance
        time: computation time
    """
    t0 = timeit.default_timer()
    cnt = len(plg.exterior.coords) - 1
    sumdist = 0
    for i in range(cnt):
        pnt1x = plg.exterior.coords[i][0] - pnt.x
        pnt1y = plg.exterior.coords[i][1] - pnt.y
        pnt2x = plg.exterior.coords[i + 1][0] - pnt.x
        pnt2y = plg.exterior.coords[i + 1][1] - pnt.y
        if pnt1x != pnt2x:
            a = (pnt2y - pnt1y) / (pnt2x - pnt1x)
            b = (pnt1y * pnt2x - pnt2y * pnt1x) / (pnt2x - pnt1x)
            sumdist = sumdist + \
                dblquad(EucDistToZero, pnt1x, pnt2x,
                        lambda u: 0, lambda u: a * u + b)[0]
    t1 = timeit.default_timer()
    if sumdist < 0:
        sumdist = -sumdist
    return (sumdist / plg.area, t1 - t0)


def DistCalcPARTNonUnif(pnt, plg, pdf):
    """[summary]
        Calculate distance using partition method with non-uniform pdf (numerical integration).
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from
        pdf (func): probability density function

    Returns:
        distance: distance
        time: computation time
    """
    # important: move the polygon to (0,0) for more accurate result
    t0 = timeit.default_timer()
    cnt = len(plg.exterior.coords) - 1
    sumdist = 0
    dists = []
    pntx = pnt.x - plg.bounds[0]
    pnty = pnt.y - plg.bounds[1]
    distcalc = functools.partial(DistToPntNonUnif, pdf, pntx, pnty)
    for i in range(cnt):
        pnt1x = plg.exterior.coords[i][0] - plg.bounds[0]
        pnt1y = plg.exterior.coords[i][1] - plg.bounds[1]
        pnt2x = plg.exterior.coords[i + 1][0] - plg.bounds[0]
        pnt2y = plg.exterior.coords[i + 1][1] - plg.bounds[1]
        if pnt1x != pnt2x:
            a = (pnt2y - pnt1y) / (pnt2x - pnt1x)
            b = (pnt1y * pnt2x - pnt2y * pnt1x) / (pnt2x - pnt1x)
            #sumdist = sumdist+dblquad(distcalc, pnt1x, pnt2x, lambda u: 0, lambda u: a*u+b)[0]
            dists.append(dblquad(distcalc, pnt1x, pnt2x,
                                 lambda u: 0, lambda u: a * u + b)[0])
    t1 = timeit.default_timer()
    # if sumdist < 0:
    #    sumdist = -sumdist
    sumdist = abs(fsum(dists))
    return (sumdist, t1 - t0)


def CumulativeDensity(plg, pdf):
    """[summary]
        Calculate cumulative probability density within polygon
    Args:
        plg (shapely.geometry.Polygon): Polygon which the distance compute from
        pdf (func): probability density function
    Returns:
        sumcd: cumulative density in polygon
    """
    cnt = len(plg.exterior.coords) - 1
    sumcd = 0
    for i in range(cnt):
        pnt1x = plg.exterior.coords[i][0] - plg.bounds[0]
        pnt1y = plg.exterior.coords[i][1] - plg.bounds[1]
        pnt2x = plg.exterior.coords[i + 1][0] - plg.bounds[0]
        pnt2y = plg.exterior.coords[i + 1][1] - plg.bounds[1]
        if pnt1x != pnt2x:
            a = (pnt2y - pnt1y) / (pnt2x - pnt1x)
            b = (pnt1y * pnt2x - pnt2y * pnt1x) / (pnt2x - pnt1x)
            sumcd = sumcd + dblquad(pdf, pnt1x, pnt2x,
                                    lambda u: 0, lambda u: a * u + b)[0]
            # print(sumcd)
    if sumcd < 0:
        sumcd = -sumcd
    return sumcd


def DistCalcCF(pnt, plg):
    """[summary]
        Calculate distance using closed-form method.
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        distance: distance
        time: computation time
    """
    t0 = timeit.default_timer()
    cnt = len(plg.exterior.coords) - 1
    sumdist = 0
    dist = [0 for i in range(2 * cnt)]
    for i in range(cnt):
        pnt1x = plg.exterior.coords[i][0] - pnt.x
        pnt1y = plg.exterior.coords[i][1] - pnt.y
        pnt2x = plg.exterior.coords[i + 1][0] - pnt.x
        pnt2y = plg.exterior.coords[i + 1][1] - pnt.y

        if pnt1x != pnt2x:
            a = (pnt2y - pnt1y) / (pnt2x - pnt1x)
            b = (pnt1y * pnt2x - pnt2y * pnt1x) / (pnt2x - pnt1x)
            if a != 0 or b != 0:
                dist[2 * i] = EucDistIntegral(a, b, pnt2x)
                dist[2 * i + 1] = -EucDistIntegral(a, b, pnt1x)
    sumdist = fsum(dist)
    if sumdist < 0:
        sumdist = -sumdist
    t1 = timeit.default_timer()
    return (sumdist / plg.area, t1 - t0)


def DistCalcCFParallel(pnt, plg, pool):
    """[summary]
        Calculate distance using closed-form method in parallel mode 1.
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        distance: distance
        time: computation time
    """
    t0 = timeit.default_timer()
    cnt = len(plg.exterior.coords) - 1
    sumdist = 0
    params = []
    for i in range(cnt):
        pnt1x = plg.exterior.coords[i][0] - pnt.x
        pnt1y = plg.exterior.coords[i][1] - pnt.y
        pnt2x = plg.exterior.coords[i + 1][0] - pnt.x
        pnt2y = plg.exterior.coords[i + 1][1] - pnt.y

        if pnt1x != pnt2x:
            a = (pnt2y - pnt1y) / (pnt2x - pnt1x)
            b = (pnt1y * pnt2x - pnt2y * pnt1x) / (pnt2x - pnt1x)
            if a != 0 or b != 0:
                params.append((a, b, pnt2x, 1))
                params.append((a, b, pnt1x, -1))

    dist = pool.map(EucDistIntegral2, params)
    sumdist = fsum(dist)
    if sumdist < 0:
        sumdist = -sumdist
    t1 = timeit.default_timer()
    return (sumdist / plg.area, t1 - t0)


def DistCalcMABBCF(pnt, plg):
    """[summary]
        Calculate distance using MABB closed-form method.
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        distance: distance
        time: computation time
    """
    # calculate distance using MABB CF form integration
    t0 = timeit.default_timer()
    mabb = plg.minimum_rotated_rectangle
    x0 = mabb.exterior.coords[0][0]
    y0 = mabb.exterior.coords[0][1]
    x1 = mabb.exterior.coords[1][0]
    y1 = mabb.exterior.coords[1][1]
    angle = arctan2(y1 - y0, x1 - x0)
    rmabb = affinity.rotate(mabb, -angle, origin=pnt, use_radians=True)
    # print(rmabb)
    minu = rmabb.bounds[0] - pnt.x
    minv = rmabb.bounds[1] - pnt.y
    maxu = rmabb.bounds[2] - pnt.x
    maxv = rmabb.bounds[3] - pnt.y
    val = EucDisRectIntegral(maxu, maxv) - EucDisRectIntegral(maxu, minv) - \
        EucDisRectIntegral(minu, maxv) + EucDisRectIntegral(minu, minv)
    # val=nquad(EucDist,[[minu,maxu],[minv,maxv]],args=(pnt.x,pnt.y))
    # print(val[1])
    t1 = timeit.default_timer()
    return (val / rmabb.area, t1 - t0)


def DistCalcRPT(pnt, plg, method=1000, z=1.6449, moe=10):
    """[summary]
        Calculate distance using random points method.
        
    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from
        method (int or str, optional): method use to calculate distance. Defaults to 1000.
                                "nmax": use nmax random points
                                "nhat": use nhat random points
                                int: use specified random points
        z (float, optional): z-value of specified confidence level. Defaults to 1.6449 (CL=90%).
        moe (int, optional): Margin of error. Defaults to 10.

    Returns:
        distance: distance
        time: computation time
    """
    t0 = timeit.default_timer()
    cnt = 0
    if type(method) is int:
        cnt = method
    elif method == "nmax":
        cnt = getnmax(pnt, plg, z, moe)
    elif method == "nhat":
        cnt = getnhat(pnt, plg, z, moe)
    else:
        return (np.nan, np.nan)
    dist = []
    rpts = RndPointinPolygon.getRndPointinPolygon(plg, cnt)
    for p in rpts:
        dist.append(EucDist(p.x, p.y, pnt.x, pnt.y))
    dist = fsum(dist) / cnt
    t1 = timeit.default_timer()
    return (dist, t1 - t0)


def getDistSD(pnt, plg, cnt=50):
    """[summary]
        calculate standard deviation of distance to random points inside a polygon

    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from
        cnt (int, optional): number of random points. Defaults to 50.

    Returns:
        std: standard deviation
    """
    dist = []
    count = 0
    rpts = RndPointinPolygon.getRndPointinPolygon(plg, cnt)
    for p in rpts:
        dist.append(EucDist(p.x, p.y, pnt.x, pnt.y))
    return std(dist, ddof=1)


def DistCTRID(pnt, plg):
    """[summary]
        Calculate distance to polygon centroid.

    Args:
        pnt ([type]): [description]
        plg ([type]): [description]

    Returns:
        [type]: [description]
    """
    t0 = timeit.default_timer()
    ctrid = plg.centroid
    dist = EucDist(pnt.x, pnt.y, ctrid.x, ctrid.y)
    t1 = timeit.default_timer()
    return (dist, t1 - t0)


def DistNEAR(pnt, plg):
    """[summary]
        Calculate closest distance.

    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        distance: distance
        time: computation time
    """
    t0 = timeit.default_timer()
    dist = plg.distance(pnt)
    t1 = timeit.default_timer()
    return (dist, t1 - t0)


def EucDist(u, v, u0, v0):
    """[summary]
        Calculate Euclidean distance between two points.

    Args:
        u (float): u-coordinate of point
        v (float): v-coordinate of point
        u0 (float): u-coordinate of reference point
        v0 (float): v-coordinate of reference point

    Returns:
        dist: Euclidian distance
    """
    return sqrt((u0 - u)**2 + (v - v0)**2)


def EucDist2(v, u, u0, v0):
    """[summary]
        Calculate Euclidean distance between two points.

    Args:
        v (float): v-coordinate of point
        u (float): u-coordinate of point
        u0 (float): u-coordinate of reference point
        v0 (float): v-coordinate of reference point

    Returns:
        dist: Euclidian distance
    """
    # reverse the u, v parameters
    return sqrt((u0 - u)**2 + (v - v0)**2)


def EucDistToZero(v, u):
    """[summary]
        Calculate Euclidean distance to (0,0).

    Args:
        v (float): v-coordinate of point
        u (float): u-coordinate of point
    Returns:
        dist: Euclidian distance
    """
    return sqrt(u**2 + v**2)


def DistToPntNonUnif(f, pntx, pnty, v, u):
    """[summary]
        Calculate Euclidean distance between two points given pdf.

    Args:
        f (func): probability density function
        pntx (float): x coordinate of point
        pnty (float): y coordinate of point
        v (float): v-coordinate of the reference point
        u (float): u-coordinate of the reference point

    Returns:
        distf: distance*pdf
    """
    return sqrt((u - pntx)**2 + (v - pnty)**2) * f(v, u)


def EucDistIntegral2(t):
    """[summary]
        Calculate Integrated Euclidean distance.

    Args:
        t (tuple): (a,b,x,sign)

    Returns:
        distance
    """
    if t is not None:
        return t[3] * EucDistIntegral(t[0], t[1], t[2])
    else:
        return 0


def EucDistIntegral(a, b, x):
    """[summary]
         Calculate Integrated Euclidean distance.

    Args:
        a (float): a value
        b (float): b value 
        x (float): x value

    Returns:
        val: Integration result
    """
    asq = a * a
    bsq = b * b
    xsq = x * x
    dn = (6 * (1 + asq)**(3 / 2))
    cx = (a * b + x + asq * x) / \
        sqrt((bsq + 2 * a * b * x + (1 + asq) * xsq)) / sqrt((1 + asq))
    if abs(abs(cx) - 1) <= 1E-9 or np.isnan(cx):
        c1 = x * b**2
    else:
        c1 = b**3 * arctanh(np.float(cx))
    c2 = sqrt(bsq + 2 * a * b * x + (1 + asq) * xsq) * \
        (2 * b * x + 2 * asq * b * x + a**3 * xsq + a * (bsq + xsq))
    if x == 0:
        c4 = 0
    else:
        c3 = abs(x) / (b + a * x + sqrt(xsq + (b + a * x)**2))
        if np.isnan(c3) or np.isinf(c3):
            if b == 0:
                c3 = 1 / (sign(x) * a + sqrt(asq + 1))
            else:
                c3 = -2 * b / abs(x)
        c4 = (1 + asq) * x**3 * log(c3)
    return (c1 + sqrt(1 + asq) * (c2 - c4)) / dn


def EucDisRectIntegral(x, y):
    """[summary]
         Calculate Integrated Euclidean distance over rectangle.

    Args:
        x (float): x value
        y (float): x value

    Returns:
        val: integration result
    """
    if x == 0 and y == 0:
        return 0
    c1 = sqrt(x**2 + y**2)
    # return (-x**3+6*x*y*c1+3*y**3*log(x+c1)+3*x**3*log(y+c1))/18
    return (6 * x * y * c1 + 3 * y**3 * log(x + c1) + 3 * x**3 * log(y + c1)) / 18


def EucDistWithBD(u, v, u0, v0, plg):
    """[summary]
         Calculate distance given polygon boundary.

    Args:
        u (float): u-coordinate of point
        v (float): v-coordinate of point
        u0 (float): u-coordinate of reference point
        v0 (float): v-coordinate of reference point
        plg (shapely.geometry.Polygon): Polygon which the distance compute from

    Returns:
        dist: distance if u,v in plg, otherwise 0.
    """
    if plg.contains(Point(u, v)):
        return EucDist(u, v, u0, v0)
    else:
        return 0
    return 0


def getnmax(pnt, plg, z, moe):
    """[summary]
        calculate nmax

    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from
        z (float, optional): z-value of specified confidence level. Defaults to 1.6449 (CL=90%).
        moe (int, optional): Margin of error. Defaults to 10.

    Returns:
        nmax: number of nmax points
    """
    dmin = plg.distance(pnt)
    dmax = plg.hausdorff_distance(pnt)
    nmax = np.ceil(z**2 * (dmax - dmin)**2 / 4 / moe / moe)
    return int(nmax)


def getnhat(pnt, plg, z, moe):
    """[summary]
        calculate nhat

    Args:
        pnt (shapely.geometry.Point): Point which the distance compute from
        plg (shapely.geometry.Polygon): Polygon which the distance compute from
        z (float, optional): z-value of specified confidence level. Defaults to 1.6449 (CL=90%).
        moe (int, optional): Margin of error. Defaults to 10.

    Returns:
        nmax: number of nhat points
    """
    sd = getDistSD(pnt, plg, cnt=50)
    nhat = np.ceil(z**2 * sd**2 / moe / moe)
    return int(nhat)


if __name__ == "__main__":
    pass
