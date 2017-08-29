from math import sqrt, sin, cos


def deg2cartesian(phi, lamda, h):
    """convert lot & lag to earth centered xyz position"""
    f = 1 / 298.257223563
    e = sqrt(f * (2 - f))
    N = 6378.137 / (2 * (1 - e ** 2 * (sin(phi) ** 2)))

    x = (N + h) * cos(phi) * cos(lamda)
    y = (N + h) * cos(phi) * sin(lamda)
    z = (N * (1 - e ** 2) + h) * sin(phi)

    return (x,y,z)


