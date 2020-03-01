#!/usr/bin/env python

import numpy as np
from math import *
import random
import matplotlib.pyplot as plt
from scipy.integrate import simps

interval = 0.3 # deg

def get_rand_ini(theta, E):

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

    zenith_list = np.linspace(radians(degrees(zenith_min) - interval), radians(degrees(zenith_max) + interval), 5000)

    factor = 1 / (1 + 1.1 * E * np.cos(zenith_list) / (115 * pow(10, 9))) + 0.054 / (
            1 + 1.1 * E * np.cos(zenith_list) / (850 * pow(10, 9)))

    factor = np.array(factor) / np.sum(np.array(factor))

    zenith_dir = np.random.choice(zenith_list, p=factor)

    """---------azimuth part---------"""
    azi_half = radians(degrees(asin(500/b)) + interval)
    azimuth_dir = np.random.random() * azi_half * 2 + (pi - azi_half)


    solid_ang = (cos(zenith_min) - cos(zenith_max)) * (2 * azi_half)

    return zenith_dir, azimuth_dir, solid_ang


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


#theta = radians(84.84)
#print(getDis(theta))
#print(getnewE(pow(10, 9), getDis(theta) - 1500))





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













"""--------------------------------------------- get random E--------------------------------------------"""

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

print(getnewE_inc(pow(10, 9), getDis(radians(72.54))))
"""

