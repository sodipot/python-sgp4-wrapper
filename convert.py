from math import sqrt, sin, cos, atan, pi
from sgp4.propagation import _gstime
from sgp4.ext import jday
from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv

import numpy as np


def deg2cartesian(lat, lon, h):
    """convert lot & lag to earth centered xyz position"""
    f = 1 / 298.257223563
    a = 6378137
    phi = np.deg2rad(lat) 
    lamda = np.deg2rad(lon) 
    e = sqrt(2 * f - f * f)
    N = a / sqrt(1 - (e * sin(phi)) ** 2)

    x = (N + h) * cos(phi) * cos(lamda)
    y = (N + h) * cos(phi) * sin(lamda)
    z = (N * (1 - e ** 2) + h) * sin(phi)

    return (x, y, z)

def deg2geometry(position,a=6378137.0,f=(1/298.257222101)):
    """this function convert decalt coordinates to Geodetic coordinates"""
    """http://www2.nc-toyama.ac.jp/~mkawai/htmlsample/coordinate/transformcg/index.html"""
    L = np.rad2deg(np.arctan(position[1]/position[0]))
    if position[0] < 0 :
        L += 180

    r = sqrt(position[0]**2+position[1]**2)

    if position[2] >= 0:
        b = a*(1-f)
    else :
        b = a*(f-1)
    E = ((position[2]+b)*b/a-a)/r
    F = ((position[2]-b)*b/a+a)/r
    P = 4*(E*F+1)/3
    Q = 2*(E**2-F**2)
    D = P**3+Q**2

    if D >= 0:
        S = sqrt(D) + Q
        S = (S/np.absolute(S))*np.exp(np.log(np.absolute(S))/3)
        V = (-1)*(2*Q+(P/S-S)**3)/(3*P)
    else :
        V = 2*sqrt((-1)*P)*np.cos(np.arccos(Q/sqrt((-1)*P**3))/3)

    G = (E + sqrt(E**2+V))/2
    T = sqrt(G**2 + (F-V*G)/(2*G-E))-G
    B = np.rad2deg(np.arctan((1-T**2)*a/(2*b*T)))
    H = (r-a*T)*np.cos(np.deg2rad(B))+(position[2]-b)*np.sin(np.deg2rad(B))
    coordinate = np.array([B,L,H])

    return coordinate

def sat2cartesian(xs, ys, zs, year, month, day, hour, minute, sec):
    """
    Convert sgp4 output(x faces vernal equinox) to earth centerd
    cartesian position(x faces prime meridian)
    """

    Tg = _gstime(jday(year, month, day, hour, minute, sec))

    conv = np.array([[np.cos(Tg), np.sin(Tg), 0],[-np.sin(Tg), np.cos(Tg), 0], [0, 0, 1]])
    source = np.array([[xs], [ys], [zs]])
    result = np.zeros((3,1))

    np.dot(conv, source, out = result)

    result = np.ndarray.tolist(result)

    return result




