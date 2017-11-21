ITF_SatTracker

---

# Overview
Read TLE(Tow Line Elements) and caluculate a position of a satellite.

## Description
This software caluculates a position(latitude & longitude) of a satellite. 
In addition, this software show us how we see the satellite from ground.
Data source is NASA Tow Line Elements(TLE). This format includes satellite's 
orbital elements. 

## Requirements
* Python 3.5 or 3.6
* sgp4 package
* numpy package

You can install sgp4 and numpy package with pip.
```
#use pip or pip3(according to your environment)
pip3 install sgp4
pip3 install numpy
```
## Usage
```python
from convert import sat2ecef, deg2ecef, calcTg, sat2direction
from sgp4.io import twoline2rv
from sgp4.earth_gravity import wgs84
from datetime import datetime, timezone
from time import sleep

line1 = '1 41932U 98067KU  17317.56946614  .00017631  00000-0  17996-3 0  9992'
line2 = '2 41932 051.6402 002.0561 0004752 087.9238 272.2302 15.65153909046886'

satellite = twoline2rv(line1, line2, wgs84)

observeTime = datetime.now(timezone.utc)
position, velocity = satellite.propagate(observeTime.year, observeTime.month, observeTime.date, observeTime.hour, observeTime.minute, observeTime.second)

Tg = calcTg(observeTime) 
satECEF = sat2ecef(position, Tg)
result = sat2direction([36,140,30], satECEF)
print(observeTime.strftime("%c")+" Azimuth:"+str(result[0]) + " Elevation:"+str(result[1]))
sleep(1.0)
```
## Output example
```
Wed Nov 15 09:19:47 2017 Azimuth:333.319747494 Elevation:-68.9659189319
```