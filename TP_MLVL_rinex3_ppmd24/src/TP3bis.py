# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""
import os
import numpy as np
import gpsdatetime as gpst




""" Import du module orbits """
# version locale
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import orbits

# version installee
import gnsstoolbox.orbits as orbits

def tp3_pos_sat_brdc(orb, const, prn, mjd):
    """ Calcul de la postion du satellite "const/prn"
    à un instant donné mjd """
    X_ECEF, Y_ECEF, Z_ECEF, dte = 0,0,0,0

    """ A COMPLETER (TP3) """
    Eph = orb.getEphemeris(const,prn,mjd)
    print("TOE : ", Eph.mjd)
    
    t = gpst.gpsdatetime(mjd=mjd)
    delta_t = t - gpst.gpsdatetime(mjd=Eph.mjd)
    print("t : ", t)
    print("delta_t : ", delta_t)
    
    print("demi-grand axe a=%.3fm" % (Eph.sqrt_a**2))
    print(Eph.mjd, mjd)

    n0 = np.sqrt(3.986005e14/(Eph.sqrt_a**2)**3)
    delta_n = Eph.delta_n
    n = n0 + delta_n
    print("n : ", n)

    M0 = Eph.M0
    M = M0 + n*delta_t
    print("M : ", M)

    E0 = M
    E1 = M + Eph.e*np.sin(E0)
    cpt = 0
    alpha = 10e-3 / 30.10e6
    while np.abs(E1 - E0) > alpha :
        E0 = E1
        E1 = M + Eph.e*np.sin(E0)
        cpt += 1 
    print("E_k+1 : ", E1)
    print("cpt : ", cpt)

    e = Eph.e
    v = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E1/2))
    print("v : ", v)

    a = Eph.sqrt_a**2
    r = (a)*(1 - e*np.cos(E1))
    print("r : ", r)

    omega = Eph.omega
    print("omega : ", omega)
    phi = omega + v
    print("phi : ", phi)

    delta_phi = Eph.cus*np.sin(2*phi) + Eph.cuc*np.cos(2*phi)
    print('delta_phi : ', delta_phi)  

    delta_r = Eph.crs*np.sin(2*phi) + Eph.crc*np.cos(2*phi)
    print('delta_r : ', delta_r) 

    x = (r + delta_r) * np.cos(phi + delta_phi)
    y = (r + delta_r) * np.sin(phi + delta_phi) 
    print("x : ", x)
    print("y : ", y)

    i_0 = Eph.i0
    i_dot = Eph.IDOT
    delta_i = Eph.cis*np.sin(2*phi) + Eph.cic*np.cos(2*phi)
    print("delta_i : ", delta_i)
    i = i_0 + i_dot*delta_t + delta_i
    print("i : ", i)

    Omega = Eph.OMEGA + Eph.OMEGA_DOT*delta_t
    print("Omega : ", Omega)
    rotZ = np.zeros((3,3))
    rotZ[0][0] = np.cos(Omega)
    rotZ[0][1] = -np.sin(Omega)
    rotZ[1][0] = np.sin(Omega)
    rotZ[1][1] = np.cos(Omega)
    rotZ[2][2] = 1

    rotX = np.zeros((3,3))
    rotX[0][0] = 1
    rotX[1][1] = np.cos(i)
    rotX[1][2] = -np.sin(i)
    rotX[2][1] = np.sin(i)
    rotX[2][2] = np.cos(i)

    vect = np.zeros((3,1))
    vect[0][0] = x
    vect[1][0] = y
    vect[2][0] = 0

    coord_ECI = rotZ@rotX@vect
    print("coord_ECI : ", coord_ECI)

    t = t.wsec

    ROTZ_ecef = np.zeros((3,3))
    ROTZ_ecef[0][0] = np.cos(-7.2921151467e-5*t)
    ROTZ_ecef[0][1] = -np.sin(-7.2921151467e-5*t)
    ROTZ_ecef[1][0] = np.sin(-7.2921151467e-5*t)
    ROTZ_ecef[1][1] = np.cos(-7.2921151467e-5*t)
    ROTZ_ecef[2][2] = 1

    coord_ECEF = ROTZ_ecef@coord_ECI
    print("coord_ecef : ", coord_ECEF)

    X_ECEF = coord_ECEF[0][0]
    Y_ECEF = coord_ECEF[1][0]
    Z_ECEF = coord_ECEF[2][0]


    return X_ECEF, Y_ECEF, Z_ECEF, dte, E1


if __name__ == "__main__":

    tic = gpst.gpsdatetime()

    Orb = orbits.orbit()

    """ lecture du fichier BRDC """
    dir_orb = '../data'
    Orb.loadRinexN(os.path.join(dir_orb, 'MLVL00FRA_R_20212400000_01D_GN.rnx'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpsdatetime(yyyy=2021,doy=240,dsec=5435)
    print(t)

    """ SATELLITE """
    constellation = 'G'
    prn = 14

    try:
        Eph = Orb.getEphemeris(constellation,prn,t.mjd)
        print(Eph)
        print("TOC : ",Eph.tgps.st_iso_epoch())
             
    except:
        print("Unable to find satellite !!! ")

    print("\nCalcul avec la fonction integee a la toolbox")
    X,Y,Z,dte = Orb.calcSatCoord(constellation,prn,t.mjd) 
    print("Solution pygnssToolbox\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\ndte = %.9f s" % (X,Y,Z,dte))
    
    """ impression des elements de debug """
#    print(Orb.debug.__dict__)

    print("\nCalcul avec la fonction developpee lors du TP")
    X,Y,Z,dte, E1 = tp3_pos_sat_brdc(Orb, constellation,prn,t.mjd)
    print("Solution TP\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\ndte = %.9f s" % (X,Y,Z,dte))

    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))



