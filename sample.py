from convert import sat2ecef, deg2ecef, calcTg, sat2direction
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs84
from datetime import datetime

line1 = '1 41932U 98067KU  17306.58554006  .00023005  00000-0  23835-3 0  9993'
line2 = '2 41932 051.6397 057.6982 0004387 055.3512 304.7896 15.64577202045164'

satellite = twoline2rv(line1, line2, wgs84)

position, velocity = satellite.propagate(2017, 11, 6, 5, 22, 49)

Tg = calcTg(datetime(2017,11,6,5,22,49)) 
satECEF = sat2ecef(position, Tg)
print("衛星の地心直交座標")
print(satECEF)
print("観測点の地心直交座標")
print(deg2ecef([36,140,30]))

print(sat2direction([36,140,30], satECEF))
