from math import sqrt, sin, cos, atan, pi
from sgp4.propagation import _gstime
from sgp4.ext import jday
from sgp4.io import twoline2rv

from datetime import datetime

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

def deg2ecef(observer, a = 6378137.0, f = (1 / 298.257223563)):
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

    ground_x = (N + observer[2]) * np.cos(phi) * np.cos(lamda)
    ground_y = (N + observer[2]) * np.cos(phi) * np.sin(lamda)
    ground_z = (N * (1 - e ** 2) + observer[2]) * sin(phi)

    return [ground_x, ground_y, ground_z]

def sat2ecef(satellite, date):
    """
    Calculate satellite's position in ecef(Earth Centered Earth Fixed)

    input:
        satellite(array): [x, y, z]
        date(datetime)

    output:
        position(numpy array): [x, y, z]
    """
    jdut1 = jday(date.year, date.month, date.day, date.hour, date.minute, date.second)
    Tg = _gstime(jdut1)

    conv1 = np.array([[np.cos(Tg), np.sin(Tg), 0],\
                      [-np.sin(Tg), np.cos(Tg), 0],\
                      [0, 0, 1]])
    source = np.array([[satellite[0]], [satellite[1]], satellite[2]])

    coordinates = np.dot(conv1, source)

    return coordinates


def sat2direction(viewposition, satposition, datetime, a = 6378137.0, f = (1 / 298.257223563)):
    """
    Return Satellites direction.

    input
        viewposition(touple): (latitude, longtitude, height)
        satposition(touple): (x, y, z)
        datetime:(touple): (year, month, date, hour, minute, second)

    output
        direction(touple): (azimuth, elevation)
    """

    #calculate the relative position vector and the distance from surface
    rvec = coordinate - np.array([ground_x, ground_y, ground_z]).reshape(3,1)
    Er = np.linalg.norm(rvec)

    #calcultate the relative position of satellite
    conv2 = np.array([[np.cos(lamda), np.sin(lamda), 0],\
                      [-np.sin(lamda), np.cos(lamda), 0],\
                      [0, 0, 1]])
    conv3 = np.array([[np.sin(phi), 0, -np.cos(phi)],\
                      [0, 1, 0],\
                      [np.cos(phi), 0, np.sin(phi)]])

    rdirection = np.dot(conv3, np.dot(conv2, rvec))
    if(rdirection[0] >= 0):
        azimuth = float(np.rad2deg(np.arctan(-rdirection[1] / rdirection[0])) - 180)
    else:
        azimuth = float(np.rad2deg(np.arctan(-rdirection[1] / rdirection[0])) + 180)

    E1 = rdirection[2] / Er
    elevation = float(np.arctan(E1 / np.sqrt(1 - np.power(E1, 2))))

    return (azimuth, elevation)

    