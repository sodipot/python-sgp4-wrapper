from math import sqrt, sin, cos, pi


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

    return (x,y,z)


