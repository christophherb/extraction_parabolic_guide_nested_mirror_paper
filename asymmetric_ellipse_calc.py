import numpy as np
import matplotlib.pyplot as plt

def calc_ks(pz, pr2, LStart, LEnd):
    """return k1 k2 and k3 such that r = sqrt(k1+z*k2+z*z*k3) from a point pz, pr**2 on the mirror surface
    and the position of both focal points

    Args:
        pz (float): z-coordinate of the point on the mirror
        pr2 (float): x**2+y**2, square of the distance of the point on the mirror to the optical axis
        LStart (float): first focal point
        LEnd (float): second focal point

    Returns:
        tuple: (k1, k2, k3) tuple such that r = sqrt(k1+z*k2+z*z*k3)
"""
    c = (LEnd-LStart)/2
    u = (pz+c-LEnd)
    a = (u*u+c*c+pr2+((u*u+c*c+pr2)**2-4*c*c*u*u)**0.5)**0.5/2**0.5
    k3 = c*c/(a*a)-1
    k2 = 2*k3*(c-LEnd)
    k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a
    return k1, k2, k3

def return_r(z, k1, k2, k3):
    return (k1+ k2*z+k3*z**2)**0.5

nummirrors = 10

r_at_z0 = []
r_at_zend = []
#r_at_lStart = []
pz0 = 0.0
r0 = 0.218/2

z0 = 0

pr2 = r0**2
LStart = -5000
LEnd = 6

lStart = -0
lEnd = 0.6

for mirror_ind in range(nummirrors):
    k1, k2, k3 = calc_ks(pz0, pr2, LStart=LStart, LEnd=LEnd)
    r_at_z0.append((k1+z0*k2+z0**2*k3)**0.5)
    r2lEnd = (k1 + lEnd*k2 + lEnd*lEnd*k3)
    rlEnd = r2lEnd**0.5
    r_at_zend.append(rlEnd)
    rlStart = rlEnd*(lStart-LStart)/(lEnd-LStart)
    pr2 = rlStart**2
    pz0 = lStart

#------------------------ Testzone -----------#
def return_all_b0s(nummirrors, z0, r0, LStart, LEnd, lStart, lEnd):
    r_at_z0 = []
    r_at_zend = []
    pz0 = z0
    pr2 = r0**2
    for mirror_ind in range(nummirrors):
        k1, k2, k3 = calc_ks(pz0, pr2, LStart=LStart, LEnd=LEnd)
        r_at_z0.append((k1+z0*k2+z0**2*k3)**0.5)
        r2lEnd = (k1 + lEnd*k2 + lEnd*lEnd*k3)
        rlEnd = r2lEnd**0.5
        r_at_zend.append(rlEnd)
        rlStart = rlEnd*(lStart-LStart)/(lEnd-LStart)
        pr2 = rlStart**2
        pz0 = lStart
    return r_at_z0, r_at_zend

# calculation of all the employed mirrors in the ess extraction mirror code
focal_length = 6
incoming_length = 5000
z0 = 0
lStart = 0
lEnd = 0.6

#mirror 1
mirror_ind = 1
print('this is mirror {}'.format(mirror_ind))
mirror_ind+=1
r0 = 0.10815692497298013
LStart = -(focal_length-0.6)
LEnd =incoming_length+0.6
print(return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd))

#mirror 2
print('this is mirror {}'.format(mirror_ind))
mirror_ind+=1
r0 = 0.1087016643197862
LStart = -focal_length
LEnd = incoming_length
print(return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd))
#
#mirror3
print('this is mirror {}'.format(mirror_ind))
mirror_ind+=1
r0 = 0.218/2+0.005
LStart = -incoming_length+0.6
LEnd = focal_length+0.6
print(return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd))
#
#
#mirror 4
print('this is mirror {}'.format(mirror_ind))
mirror_ind+=1
r0 = 0.218/2+0.005
LStart = -incoming_length
LEnd = focal_length
print(return_all_b0s(10, z0, r0, LStart, LEnd, lStart, lEnd))