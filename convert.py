from math import sqrt, sin, cos, atan, pi
from sgp4.propagation import _gstime
from sgp4.ext import jday
from sgp4.io import twoline2rv

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

    #calculate the viewer's position in xyz
    phi = np.deg2rad(viewposition[0]) 
    lamda = np.deg2rad(viewposition[1]) 
    e = np.sqrt(2 * f - f * f)
    N = a / np.sqrt(1 - (e * sin(phi)) ** 2)

    ground_x = (N + viewposition[2]) * np.cos(phi) * np.cos(lamda)
    ground_y = (N + viewposition[2]) * np.cos(phi) * np.sin(lamda)
    ground_z = (N * (1 - e ** 2) + viewposition[2]) * sin(phi)
    print("観測点の地心直交座標は" + str((ground_x, ground_y, ground_z)))


    #calculate the satellite's position in xyz(x faces prime meridian)
    Tg = np.deg2rad(_gstime(jday(datetime[0], datetime[1], datetime[2], datetime[3], datetime[4], datetime[5])))

    conv1 = np.array([[np.cos(Tg), np.sin(Tg), 0],[-np.sin(Tg), np.cos(Tg), 0], [0, 0, 1]])
    source = np.asarray(satposition).reshape(3,1)

    coordinate = np.dot(conv1, source)
    print("衛星の地心直交座標は" + str(coordinate))


    #calculate the relative position vector and the distance from surface
    rvec = coordinate - np.array([ground_x, ground_y, ground_z]).reshape(3,1)
    Er = np.linalg.norm(rvec)
    print("衛星までの距離は" + str(Er) + "m")

    #calcultate the relative position of satellite
    conv2 = np.array([[np.cos(lamda), np.sin(lamda), 0], [-np.sin(lamda), np.cos(lamda), 0], [0, 0, 1]])
    conv3 = np.array([[np.sin(phi), 0, -np.cos(phi)], [0, 1, 0], [np.cos(phi), 0, np.sin(phi)]])

    rdirection = np.dot(conv3, np.dot(conv2, rvec))
    if(rdirection[0] >= 0):
        print("x>=0であった。")
        print("未補正方位角は" + str(float(np.rad2deg(np.arctan(-rdirection[1] / rdirection[0])))))
        azimuth = float(np.rad2deg(np.atan(-rdirection[1] / rdirection[0])) - 180)
    else:
        print("x<0であった。")
        print("未補正方位角は" + str(float(np.rad2deg(np.arctan(-rdirection[1] / rdirection[0])))))
        azimuth = float(np.rad2deg(np.arctan(-rdirection[1] / rdirection[0])) + 180)

    E1 = rdirection[2] / Er
    elevation = float(np.arctan(E1 / np.sqrt(1 - np.power(E1, 2))))

    return (azimuth, elevation)


if __name__ == "__main__":
    from sgp4.earth_gravity import wgs84
    from sgp4.io import twoline2rv

    viewposition = (35.885908, 139.844885, 40.883)
    datetime = (2017, 10, 4, 15, 13, 50)
    
    line1 = '1 41932U 98067KU  17272.13766631  .00039990  00000-0  44145-3 0  9992'
    line2 = '2 41932 051.6389 231.8839 0003424 312.5715 047.4990 15.62655094039776'

    satellite = twoline2rv(line1, line2, wgs84)
    position, velocity = satellite.propagate(2017, 10, 4, 15, 13, 50)
    if(satellite.error == 0):
        azimuth, elevation = sat2direction(viewposition, position, datetime)
        print("衛星名称: TO-89")
        print("方位角: " + str(azimuth) + "度")
        print("仰角: " + str(elevation)+ "度")
    else:
        print("error at propagate")
    