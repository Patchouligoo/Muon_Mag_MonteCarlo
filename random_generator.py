#!/usr/bin/env python

import numpy as np
from math import *
import random
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import hyp2f1
#interval = 0.3 # deg

def getinterval(max, min):
    interval = (max - min)/4
    return interval


def integrand(theta, E_max, E_min):

    temp1 = 0.14 * -0.588235 * (E_max **-1.7) * hyp2f1(-1.7, 1, -0.7, -1.11 * cos(theta) * E_max/(115 * pow(10, 9)))

    temp2 = 0.14 * 0.054 * -0.588235 * (E_max **-1.7) * hyp2f1(-1.7, 1, -0.7, -1.11 * cos(theta) * E_max/(850 * pow(10, 9)))

    temp3 = 0.14 * -0.588235 * (E_min **-1.7) * hyp2f1(-1.7, 1, -0.7, -1.11 * cos(theta) * E_min/(115 * pow(10, 9)))

    temp4 = 0.14 * 0.054 * -0.588235 * (E_min **-1.7) * hyp2f1(-1.7, 1, -0.7, -1.11 * cos(theta) * E_min/(850 * pow(10, 9)))

    return temp1 + temp2 - temp3 - temp4

def get_P(zenith_min, zenith_max, azi_half, E_min, E_max):

    integral = quad(integrand, zenith_min, zenith_max, args=(E_max, E_min))[0] * 2 * azi_half

    return integral


def check_in_range(E, theta, rand):
    val = (0.14 * pow(E, -2.7)) * ( 1/ (1 + 1.1 * E * np.cos(theta) / (115 * pow(10, 9)))
                              + 0.054 / ( 1 + 1.1 * E * np.cos(theta) / (850 * pow(10, 9))) )
    if rand < val:
        return True
    else:
        return False


def get_rand_ini(theta, E_min, E_max):

    len1 = 6369500
    len2 = 6369000
    R = 6371000
    A = pi - theta

    b = len2 * cos(A) + sqrt(pow(len2, 2) * pow(cos(A), 2) - (pow(len2, 2) - pow(R, 2)))

    zenith_interval = asin(500 / b) * 2

    B = acos((R**2 + len2**2 - b**2)/(2 * R * len2))

    base_ang = (theta - zenith_interval/2) - B

    zenith_min = base_ang
    zenith_max = base_ang + zenith_interval

    factor = 0.14 * pow(E_min, -2.7) * (1 + 0.054)

    checker = False

    zenith_dir = 0
    energy = 0

    while not checker:
        interval = degrees(getinterval(zenith_max, zenith_min))
        energy = np.random.random() * (E_max - E_min) + E_min
        zenith_dir = np.random.random() * (zenith_max - zenith_min + 2*radians(interval)) + zenith_min - radians(interval)
        curr = np.random.random() * factor

        checker = check_in_range(energy, zenith_dir, curr)

    """---------azimuth part---------"""
    azi_half = radians(degrees(asin(500/b)) + interval)
    azimuth_dir = np.random.random() * azi_half * 2 + (pi - azi_half)


    solid_ang = (cos(zenith_min - radians(interval)) - cos(zenith_max + radians(interval))) * (2 * (azi_half + radians(interval)))

    integral = get_P(zenith_min - radians(interval), zenith_max + radians(interval), azi_half + radians(interval), E_min, E_max)
    base = get_P(0, pi/2, pi, 1e-10, E_max * 10)

    

    return energy, zenith_dir, azimuth_dir, solid_ang, integral, base


"""-------------------------------- get radndom pos ------------------------------------"""
def convert_theta(icecube_theta):
    len2 = 6369000
    R = 6371000
    A = pi - icecube_theta

    b = len2 * cos(A) + sqrt(pow(len2, 2) * pow(cos(A), 2) - (pow(len2, 2) - pow(R, 2)))

    cos_theta = (R**2 + len2**2 - b**2) / (2 * R * len2)

    theta = acos(cos_theta)

    return theta


def get_rand_spherical(theta_min, theta_max):

    theta_min = convert_theta(theta_min)
    theta_max = convert_theta(theta_max)

    phi = 2*pi*random.random()

    u_min = (cos(theta_min) + 1)/2

    u_max = (cos(theta_max) + 1)/2

    u = random.random() * (u_max - u_min) + u_min

    theta = acos(2*u - 1)

    return theta, phi


def convert_final_dir(dir):

    dir = -1 * dir

    x = dir[0]
    y = dir[1]

    z_axis = np.array([0, 0, 1])

    theta = acos(np.dot(dir, z_axis)/np.linalg.norm(dir))

    rho = np.array([dir[0], dir[1], 0])

    rho = rho/np.linalg.norm(rho)

    x_axis = np.array([1, 0, 0])

    phi = acos(np.dot(rho, x_axis))

    if x < 0 and y < 0:
        phi = 2*pi - phi

    elif x > 0 and y < 0:
        phi = 2*pi - phi

    return theta, phi


"""
list = []

for i in range(1000):
    print(i)
    list.append(get_rand_ini(pi/2, pow(10, 15), pow(10, 13))[0])

plt.hist(list)

plt.show()
"""

"""-----------------------------------------------------------------------------------------"""

"""
def getnewE_inc(E, dis):
    b = 0.000363

    a = 0.259

    density = 0.9167
    # todo: check this: the unit for a and b are /mwe, which is meter water equivalence, as the result
    # todo: density is used here to convert into ice.
    # todo: however on the original paper a and b are plotted for ice, so I don't know if density is needed
    # todo: paper is at  FIG.21 from arXiv:hep-ph/0407075v3

    # convert E in GeV to use this eq
    A = (E / pow(10, 9) * b * density) + a * density
    # convert E_final back in eV
    E = (A * np.power(e, b * 1 * dis * density) - a * density) / (b * density) * pow(10, 9)

    return E

def getDis(theta):

    len2 = 6369000
    R = 6371000
    A = pi - theta

    b = len2 * cos(A) + sqrt(pow(len2, 2) * pow(cos(A), 2) - (pow(len2, 2) - pow(R, 2)))

    return b

cos_list = np.linspace(1, 0, 1000)
theta_list = np.arccos(cos_list)

E_list = []
for i in range(len(cos_list)):
    E_list.append(getnewE_inc(pow(10, 9), getDis(theta_list[i])))

f = plt.figure()
plt.plot(cos_list, E_list)
plt.yscale('log')
plt.xlabel('cos_theta')
plt.ylabel('E log scale')
plt.xlim(1, 0)
plt.grid()
plt.show()

"""
