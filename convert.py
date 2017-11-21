from math import sqrt, sin
from sgp4.propagation import _gstime
from sgp4.ext import jday

import numpy as np


def deg2geometry(position,a=6378137.0,f=(1/298.257223563)):
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

def calcTg(date):

    jdut1 = jday(date.year, date.month, date.day, date.hour, date.minute, date.second)
    Tg = _gstime(jdut1)

    return Tg


def deg2ecef(observer, a=6378137.0, f=(1 / 298.257223563)):
    """
    This function returns observer's position in cartesian.

    input:
        observer(array): [latitude, longtitude, height]

    output: [x, y, z]
    """
    phi = np.deg2rad(observer[0])
    lamda = np.deg2rad(observer[1])
    e = np.sqrt(2 * f - f * f)
    N = a / np.sqrt(1 - (e * sin(phi)) ** 2)

    ground_x = (N + observer[2]) * np.cos(phi) * np.cos(lamda) / 1000
    ground_y = (N + observer[2]) * np.cos(phi) * np.sin(lamda) / 1000
    ground_z = (N * (1 - e ** 2) + observer[2]) * sin(phi) / 1000

    return [ground_x, ground_y, ground_z]

def sat2ecef(satellite, Tg):
    """
    Calculate satellite's position in ecef(Earth Centered Earth Fixed)

    input:
        satellite(array): [x, y, z]

    output:
        position(array): [x, y, z]
    """

    conv1 = np.matrix([[np.cos(Tg), np.sin(Tg), 0],\
                      [-np.sin(Tg), np.cos(Tg), 0],\
                      [0, 0, 1]])
    source = np.matrix([satellite[0], satellite[1], satellite[2]]).transpose()

    coordinates = conv1 * source
    return [coordinates[0, 0], coordinates[1, 0], coordinates[2, 0]]


def sat2direction(observer, satposition):
    """
    Return Satellites direction.

    input
        observer(touple): (latitude, longtitude, height)
        satposition(touple): (x, y, z)

    output
        direction(touple): (azimuth, elevation)
    """

    phi = np.deg2rad(observer[0])
    lamda = np.deg2rad(observer[1])
    height = observer[2]
    observer_xyz = deg2ecef(observer)

    rvec = np.asmatrix(satposition) - np.asmatrix(observer_xyz)
    Er = np.linalg.norm(rvec)
    rvec = rvec.transpose()

    #calcultate the relative position of satellite
    conv2 = np.matrix([[np.cos(lamda), np.sin(lamda), 0],\
                      [-np.sin(lamda), np.cos(lamda), 0],\
                      [0, 0, 1]])
    conv3 = np.matrix([[np.sin(phi), 0, -np.cos(phi)],\
                      [0, 1, 0],\
                      [np.cos(phi), 0, np.sin(phi)]])

    rdirection = conv3 * conv2 * rvec

    azimuth = 180.0 - np.rad2deg(np.arctan2(rdirection[1, 0], rdirection[0, 0]))
    if azimuth >= 360.0:
        azimuth -= 360.0

    E1 = rdirection[2, 0] / Er
    elevation = np.arctan(E1 / np.sqrt(1 - E1 * E1))
    elevation = np.rad2deg(elevation)

    return (azimuth, elevation)