#! /usr/bin/env ipython

import scipy.io as sio
from scipy import interpolate
import numpy as np

def load(tStart, tEnd, Ts, datFolder='data/smc_20140430_232938_dmhe_testing'):
    #
    # Load data
    #

    names = ['cable_length', # 50 Hz, approx half period delayed
             'control_surfaces', # 500 Hz
             'encoder', # 1 kHz
             'imu', # 500 Hz
#            'led_tracker', # 12.5 Hz, 1 sample delayed
             'line_angle_sensor'] # 1 kHz

    rawData = {}
    for name in names:
        rawData[ name ] = sio.loadmat(datFolder + '/' + name + '.mat')
        _data = rawData[ name ]
        _indices = np.where((_data["time"] > tStart * 0.95) & (_data["time"] < tEnd * 1.05))[ 0 ]
        for k, v in _data.items():
                if k.startswith( "__" ):
                        continue
                rawData[ name ][ k ] = v[ _indices ]

    #
    # Resample
    #
    _interval = np.arange(tStart, tEnd, Ts)
    data = {}
    for name in names:
        data[ name ] = {}
    #   print name
        _time = rawData[ name ][ "time" ]
        for k in rawData[ name ].keys():
                if k == "time" or k.startswith( "__" ):
                        continue
    #           print k
                _x, _y = _time.flatten(), rawData[ name ][ k ].flatten()
                _f = interpolate.interp1d(_x, _y)
                data[ name ][ k ] = _f( _interval )
    return (data, _interval, rawData)


if __name__=='__main__':
    #
    # Testing
    #
    import matplotlib.pyplot as plt

    tStart = 0.5 # [sec]
    tEnd   = 30.0 # [sec]
    Ts = 1. / 100. #

    (data, interval, rawData) = load(tStart, tEnd, Ts)

    plt.figure()
    plt.subplot(2, 1, 1)
    plt.step(rawData['imu']['time'], rawData['imu']['gyro_x'], 'r')
    plt.step(interval, data['imu']['gyro_x'], 'b')
    plt.ylabel('gyro_x')
    plt.subplot(2, 1, 2)
    plt.step(interval, data['encoder']['speed_rpm'], 'b')
    plt.ylabel('speed_rpm')

    plt.show()
