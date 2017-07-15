# -*- coding: utf-8 -*-

def xy2bl(L0, x, y):
    import math
    x = x - 500000
    K0 = 0.157046064172e-6
    K1 = 0.005051773759
    K2 = 0.000029837303
    K3 = 0.000000238189
    PAI = 3.1415926535898
    a = 6378245.0
    e2 = 0.00669342162297
    el2 = 0.00673852541468
    p0 = 180.0 / PAI
    p2 = 3600.0 * 180.0 / PAI
    P0 = PAI / 180.0

    B0f = K0 * y * p0
    sinB0f = math.sin(B0f * P0)
    sinB0f2 = math.sin(B0f * P0) * math.sin(B0f * P0)
    Temp1 = p0 * (math.cos(B0f * P0)) * sinB0f * (K1 - K2 * sinB0f2 + K3 * sinB0f2 * sinB0f2)
    Bf = B0f + Temp1
    tf = math.tan(Bf * P0)
    Ita2 = el2 * math.cos(Bf * P0) * math.cos(Bf * P0)
    Wf = 1 - e2 * math.sin(Bf * P0) * math.sin(Bf * P0)
    Mf = a * (1 - e2) / (Wf * Wf * Wf)
    Nf = a / math.sqrt(1 - e2 * math.sin(Bf * P0) * math.sin(Bf * P0))
    lon = L0 + (p2 * x / (math.cos(Bf * P0) * Nf) + p2 * (-1. / 6.) * (1. + 2. * tf * tf + Ita2) * x * x * x / (
    Nf * Nf * Nf * math.cos(Bf * P0)) + p2 * (1. / 120.) * (
                5. + 28. * tf * tf + 24. * tf * tf * tf * tf + 6. * Ita2 + 8. * Ita2 * tf * tf) * x * x * x * x * x / (
                Nf * Nf * Nf * Nf * Nf * math.cos(Bf * P0))) / 3600.
    lat = Bf + (p2 * (-1. / 2.) * tf * x * x / (Mf * Nf) + p2 * (1. / 24.) * tf * (
    5. + 3. * tf * tf + Ita2 - 9. * Ita2 * tf * tf) * x * x * x * x / (Mf * Nf * Nf * Nf) + p2 * (-1. / 720.) * tf * (
                61. + 90. * tf * tf + 45. * tf * tf * tf * tf) * x * x * x * x * x * x / (
                Mf * Nf * Nf * Nf * Nf * Nf)) / 3600.
    lonlat = [lon, lat]
    return lonlat


def bl2xy(L0, lon, lat):
    # [x,y]=bl2xy(L0,lon,lat)
    # convert the  coor of longitude and latitude to 54 coor(xy)
    # L0 the center coor of longitude and latitude
    # lon,lat the longitude coor and latitude coor
    import math
    PAI = 3.1415926535898
    a = 6378245
    e2 = 0.00669342162297
    el2 = 0.00673852541468
    P0 = PAI / 180
    C0 = 6367558.49686
    C1 = 32005.79642
    C2 = 133.86115
    C3 = 0.7031

    B = lat * P0
    L = lon * P0
    Ita2 = el2 * math.cos(B) * math.cos(B)
    ll = (L - L0 * P0)
    t = math.tan(B)
    m0 = math.cos(B) * ll
    N = a / math.sqrt(1 - e2 * math.sin(B) * math.sin(B))
    X1 = C0 * B - math.cos(B) * math.sin(B) * (
    C1 + C2 * math.sin(B) * math.sin(B) + C3 * math.sin(B) * math.sin(B) * math.sin(B) * math.sin(B))
    y = X1 + 1. / 2. * N * t * m0 * m0 + 1. / 24. * N * t * (
    5. - t * t + 9 * Ita2 + 4 * Ita2 * Ita2) * m0 * m0 * m0 * m0 + 1 / 720 * N * t * (
    61 - 58 * t * t + t * t * t * t + 270 * Ita2 - 330 * Ita2 * t * t) * m0 * m0 * m0 * m0 * m0 * m0
    x = 500000 + N * m0 + 1. / 6. * N * (1 - t * t + Ita2) * m0 * m0 * m0 + 1. / 120. * N * (
    5 - 18 * t * t + t * t * t * t + 14 * Ita2 - 58 * Ita2 * t * t) * m0 * m0 * m0 * m0 * m0
    xy = [x, y]
    return xy


if __name__ == '__main__':
    import numpy as np

    print '===== test coordinate converting in __main__======'
    L0 = 123
    data = [[389240.35, 3465459.19],
            [428301.30, 3450155.96],
            [349394.96, 3484091.40],
            [417319.21, 3458011.34],
            [376462.15, 3451563.01],
            [376462.16, 3451563.01],
            [376462.16, 3451563.01],
            [376462.16, 3451563.01],
            [376462.23, 3451563.03],
            [376462.24, 3451563.04],
            [376526.28, 3451612.51],
            [381745.52, 3455998.57]]
    a = np.array(data)
    (m, n) = np.shape(a)
    lonlat = np.zeros([m, n])
    b = np.zeros_like(a)
    c = np.zeros_like(a)

    for i in range(m):
        lonlat = xy2bl(L0, a[i, 0], a[i, 1])
        b[i, 0] = lonlat[0]
        b[i, 1] = lonlat[1]

    for i in range(m):
        xy = bl2xy(L0, b[i, 0], b[i, 1])
        c[i, 0] = xy[0]
        c[i, 1] = xy[1]

