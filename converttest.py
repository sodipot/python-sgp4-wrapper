import convert
import sgp4
import unittest

class TestConvertMethods(unittest.TestCase):
    #TLE using:
    #TO-89
    #1 41932U 98067KU  17306.58554006  .00023005  00000-0  23835-3 0  9993
    #2 41932 051.6397 057.6982 0004387 055.3512 304.7896 15.64577202045164
    from sgp4.earth_gravity import wgs84

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_deg2cartesian(self):
        observer = [36, 140, 30]
        result = convert.deg2cartesian(observer)
        self.assertAlmostEqual(result[0], -3957000, delta = 1000)
        self.assertAlmostEqual(result[1], 3320000, delta = 1000)
        self.assertAlmostEqual(result[2], 3728000, delta = 1000)


if __name__ == '__main__':
    unittest.main()