import convert
import sgp4
import unittest

class TestConvertMethods(unittest.TestCase):
    from sgp4.earth_gravity import wgs84

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_jday(self):
        datetime = [2071, 10, 14, 15, 13, 50]
        reslut = sgp4.ext.jday(datetime[0], datetime[1], datetime[2], datetime[3], datetime[4], datetime[5]) - 2400000.5
        self.assertEqual(reslut, 58030.64346)

    def test_sat2direction(self):
        viewposition = (35.885908, 139.844885, 40.883)
        datetime = [2071, 10, 14, 15, 13, 50]
        line1 = '1 41932U 98067KU  17272.13766631  .00039990  00000-0  44145-3 0  9992'
        line2 = '2 41932 051.6389 231.8839 0003424 312.5715 047.4990 15.62655094039776'

        satellite = sgp4.io.twoline2rv(line1, line2, sgp4.earth_gravity.wgs84)
        position, velocity = satellite.propagate(2017, 10, 4, 15, 13, 50)
        position = list(position)
        position[0] *= 1000
        position[1] *= 1000
        position[2] *= 1000
    
        if(satellite.error == 0):
            azimuth, elevation = convert.sat2direction(viewposition, position, datetime)


if __name__ == '__main__':
    unittest.main()