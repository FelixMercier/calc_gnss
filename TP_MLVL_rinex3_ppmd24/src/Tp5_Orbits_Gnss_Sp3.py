# -*- coding: utf-8 -*-
"""
@author: BEILIN Jacques - IGN/ENSG
@date  : %(date)s
"""

import re
import numpy as np
import math
import os

import gpsdatetime as gpst


""" Import du module orbits """
# version locale
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import orbits

# version installee
import gnsstoolbox.orbits as orbits


def pos_sat_sp3(sp3, const, prn, mjd, ordre):
    """ Calcul de la postion du satellite "const/prn"
    à un instant donné mjd """
    X_ECEF, Y_ECEF, ZECEF, dte = 0,0,0,0

    """ A COMPLETER (TP5) """
    (orb,nl) = sp3.getSp3(const,prn)
    print(orb[0:3,:])
    print(orb.shape)

    return X_ECEF, Y_ECEF, ZECEF, dte

if __name__ == "__main__":

    tic = gpst.gpsdatetime()

    mysp3 = orbits.orbit()

    dir_orb = '../data'
    mysp3.loadSp3(os.path.join(dir_orb,'GFZ0OPSULT_20212400600_02D_05M_ORB.SP3'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpsdatetime(yyyy=2021,doy=240,dsec=5435.000)
#    print(t)

    """ SATELLITE """
    constellation = 'G'
    prn = 14

    """ Ordre du polynome """
    ordre = 9

    """ Calcul avec la fonction integree a la toolbox """
    (X_sp3,Y_sp3,Z_sp3,dte_sp3)	 = mysp3.calcSatCoordSp3(constellation,prn,t.mjd,ordre)
    print("\nSolution pygnssToolbox\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\ndte = %.9f s" % (X_sp3,Y_sp3,Z_sp3,dte_sp3))

    (X1,Y1,Z1,dte1)	= mysp3.calcSatCoordSp3(constellation,prn,t.mjd-0.001/86400,ordre)
    (X2,Y2,Z2,dte2)	= mysp3.calcSatCoordSp3(constellation,prn,t.mjd+0.001/86400,ordre)
    V = np.array([[X2-X1],[Y2-Y1],[Z2-Z1]]) / 0.002
    print("V = %.3f m/s" % np.linalg.norm(V))

    """ Fonction developpee lors du TP """
    X,Y,Z,dte = pos_sat_sp3(mysp3, constellation, prn, t.mjd, ordre)
    print("\nFonction developpee lors du TP\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\ndte = %.9f s" % (X,Y,Z,dte))

    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))

