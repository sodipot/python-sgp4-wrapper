from convert import sat2ecef, deg2ecef, calcTg, sat2direction
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs84
from datetime import datetime, timezone
from time import sleep

line1 = '1 41932U 98067KU  17317.56946614  .00017631  00000-0  17996-3 0  9992'
line2 = '2 41932 051.6402 002.0561 0004752 087.9238 272.2302 15.65153909046886'

satellite = twoline2rv(line1, line2, wgs84)
for i in range(0,60):
    observeTime = datetime.now(timezone.utc)
    #observeTime = datetime(2017,11,15,7,24,i)
    position, velocity = satellite.propagate(2017, 11, 15, observeTime.hour, observeTime.minute, observeTime.second)

    Tg = calcTg(observeTime) 
    satECEF = sat2ecef(position, Tg)
    result = sat2direction([36,140,30], satECEF)
    print(observeTime.strftime("%c")+" Azimuth:"+str(result[0]) + " Elevation:"+str(result[1]))
    sleep(1.0)
