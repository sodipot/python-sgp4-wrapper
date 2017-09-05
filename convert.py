from math import sqrt, sin, cos, atan, pi
from sgp4.propagation import _gstime
from sgp4.ext import jday
import numpy as np


def deg2cartesian(lat, lon, h):
    """convert lot & lag to earth centered xyz position"""
    f = 1 / 298.257223563
    a = 6378137
    phi = lat * pi / 180
    lamda = lon * pi / 180
    e = sqrt(2 * f - f * f)
    N = a / sqrt(1 - (e * sin(phi)) ** 2)

    x = (N + h) * cos(phi) * cos(lamda)
    y = (N + h) * cos(phi) * sin(lamda)
    z = (N * (1 - e ** 2) + h) * sin(phi)

    return (x, y, z)


def sat2cartesian(xs, ys, zs, year, month, day, hour, minute, sec):
    """
    Convert sgp4 output(x faces vernal equinox) to earth centerd
    cartesian position(x faces prime meridian)
    """

    Tg = _gstime(jday(year, month, day, hour, minute, sec))
