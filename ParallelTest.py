from avgdist import DistCalc

import fiona
from shapely.geometry import Point, Polygon, shape

import timeit
import multiprocessing as mp
from functools import partial

from os import path


def DistCalcCF(a):
    return DistCalc.DistCalcCF(a[0], a[1])


def LoadData():
    """[summary]
    Load data
    Returns:
        data: list of (point,polygon)
    """
    data = []
    shp = path.join('TestData', 'ResearchArea.shp')
    with fiona.open(shp) as source:
        for f in source:
            if f['properties']['CID'] > 0:
                plg = shape(f['geometry'])
                pnt = Point(f['properties']['REFX'], f['properties']['REFY'])
                data.append((pnt, plg))
    return data


def SingleComputeTest(data):
    """[summary]

    Args:
        data ([type]): [description]
    """
    t0 = timeit.default_timer()
    for a in data:
        DistCalcCF(a)
    t1 = timeit.default_timer()
    print("Single Thread: %ss" % (t1 - t0))


def ParallelComputeTest1(data, n):
    """[summary]

    Args:
        data ([type]): [description]
        n ([type]): [description]
    """

    pool = mp.Pool(n)

    t0 = timeit.default_timer()
    for a in data:
        DistCalc.DistCalcCFParallel(a[0], a[1], pool)
    t1 = timeit.default_timer()

    print("Parallel Compute Mode 1 %s Threads: %ss" % (n, t1 - t0))


def ParallelComputeTest2(data, n):
    """[summary]

    Args:
        data ([type]): [description]
        n ([type]): [description]
    """
    pool = mp.Pool(n)

    t0 = timeit.default_timer()
    result2 = pool.map(DistCalcCF, data)
    t1 = timeit.default_timer()
    print("Parallel Compute Mode 2 %s Threads: %ss" % (n, t1 - t0))


if __name__ == "__main__":
    data = LoadData()
    SingleComputeTest(data)


    for i in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]:
        ParallelComputeTest1(data, i)
    for i in [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]:
        ParallelComputeTest2(data, i)
    pass
