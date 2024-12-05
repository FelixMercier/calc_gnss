# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""
import re
import numpy as np
import gpsdatetime as gpst
import math
import pip
import os

""" Import du module orbits """
# version locale
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import orbits

# version installee
import gnsstoolbox.orbits as orbits

def pos_sat_rk4(orb, const, prn, mjd):
    """ Calcul de la postion du satellite "const/prn"
    à un instant donné mjd """
    X_ECEF, Y_ECEF, ZECEF, dte = 0,0,0,0

    """ A COMPLETER (TP4) """

    return X_ECEF, Y_ECEF, ZECEF, dte


if __name__ == "__main__":

    tic = gpst.gpsdatetime()

    Orb = orbits.orbit()
    
    """ lecture du fichier BRDC """
    dir_orb = '../data'
    Orb.loadRinexN(os.path.join(dir_orb, 'MLVL00FRA_R_20212400000_01D_RN.rnx'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpsdatetime(yyyy=2021,doy=240,dsec=5400)

    """ Ecart temps Glonass/GPS """ 
    t -= 18
    print(t)

    """ SATELLITE """
    constellation = 'R'
    prn = 1

    try:
        Eph = Orb.getEphemeris(constellation,prn,t.mjd)
        print(Eph)
        print("TOC : ",Eph.tgps.st_iso_epoch())
    except:
        print("Unable to find satellite !!! ")

    XSP3 =[  15877666.256,  -3409654.358, -19677679.894,     84.494967]

    """ Calcul avec la fonction integree a la toolbox """
    X,Y,Z, VX,VY,VZ,dte, deb = Orb.calcSatCoordGlonassNav(constellation,prn,t.mjd)
    print("\nSolution pygnssToolbox\nX = %13.3f \nY = %13.3f \nZ = %13.3f \ndte = %.9f us" % (X,Y,Z,dte))
    
#    print(X - XSP3[0], Y-XSP3[1], Z-XSP3[2])

    """ Calcul avec la fonction developpee lors du TP """
    X,Y,Z,dte = pos_sat_rk4(Orb, constellation, prn, t.mjd)
    print("\nSolution TP\nX = %13.3f \nY = %13.3f \nZ = %13.3f \ndte = %.9f us" % (X,Y,Z,dte))


    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))

