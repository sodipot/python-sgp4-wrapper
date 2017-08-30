ITF_SatTracker

---

# Overview
TLE(2行軌道要素)を読み込んで、任意の時間の衛星位置を計算するソフトウェア。

Read TLE(Tow Line Elements) and caluculate a position of a satellite.

## Description
NASA TLEから人工衛星の軌道要素を読み込んで、任意の時間における衛星の緯度経度、および地上の任意の地点からの
相対的な見え方(仰角、方位)や速度を計算する。

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
pip3 install sgp4
pip3 install numpy
```
