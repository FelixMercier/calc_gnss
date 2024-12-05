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
    
    
    
    Eph = orb.getEphemeris(const,prn,mjd)
    print("TOE : ", Eph.mjd)
    
    t_eph = gpst.gpsdatetime(mjd=Eph.mjd)
    t = gpst.gpsdatetime(mjd=mjd)
    
    delta_t = t - t_eph
    print()
    
    n = np.sqrt(mu / (Eph.sqrt_a**6) ) + Eph.delta_n
    print("n = ", n)

    M = Eph.M0 + n*delta_t
    print("M = ", M)
    
    E0=Eph.M0
    E = E0+1
    R = 30e6
    
    while np.abs(R*E - R*E0) > 1e-3:
        E0 = E
        E = M + Eph.e * np.sin(E0)
    print("E = ", E)
    
    v = 2*np.arctan(np.sqrt((1 + Eph.e)/(1 - Eph.e)) * np.tan(E/2))
    print("v = ", v)
    
    r = (Eph.sqrt_a**2)*(1 - Eph.e * np.cos(E))
    print("r = ", r)
    
    phi = Eph.omega + v
    dphi = Eph.cus*np.sin(2*phi) + Eph.cuc*np.cos(2*phi)
    print("dphi = " , dphi)
    
    dr = Eph.crs*np.sin(2*phi) + Eph.crc*np.cos(2*phi)
    print("dr = ", dr)
    
    x = (r + dr) * np.cos(phi + dphi)
    y = (r + dr) * np.sin(phi + dphi)
    print("x = ", x)
    print("y = ", y)
    
    di = Eph.cis*np.sin(2*phi) + Eph.cic*np.cos(2*phi)
    
    i = Eph.i0 + Eph.IDOT * delta_t + di
    Om = Eph.OMEGA + Eph.OMEGA_DOT * delta_t
    
    matRotOm = np.array([[np.cos(Om), - np.sin(Om), 0],
                       [np.sin(Om) , np.cos(Om), 0],
                       [0, 0, 1]])
    matRoti = np.array([[1, 0, 0],
                         [0, np.cos(i), -np.sin(i)],
                         [0, np.sin(i), np.cos(i)]])
    POS = matRotOm @ matRoti @ np.array([[x], [y], [0]])
    X, Y, Z = POS[0, 0], POS[1, 0], POS[2, 0]
    print("X = ", X)
    print("Y = ", Y)
    print("Z = ", Z)
    
    Om_dot_e = -7.2921151467e-5
    wk_sec = t.wsec
    
    if t.wk != Eph.tgps.wk:
        if t.wk > Eph.tgps.wk:
            wk_sec += 86400
        else:
            wk_sec -= 86400
    
    arg = Om_dot_e * wk_sec
    
    matRotOme = np.array([[np.cos(arg), - np.sin(arg), 0],
                       [np.sin(arg) , np.cos(arg), 0],
                       [0, 0, 1]])
    POS_ECEF = matRotOme @ POS
    
    X_ECEF, Y_ECEF, Z_ECEF = POS_ECEF[0, 0], POS_ECEF[1, 0], POS_ECEF[2, 0]
    
    print("demi-grand axe a=%.3fm" % (Eph.sqrt_a**2))
    print(Eph.mjd, mjd)

    return X_ECEF, Y_ECEF, Z_ECEF, dte


if __name__ == "__main__":

    mu = 3.986005e14
    
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
    X,Y,Z,dte = tp3_pos_sat_brdc(Orb, constellation,prn,t.mjd)
    print("Solution TP\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\ndte = %.9f s" % (X,Y,Z,dte))

    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))

