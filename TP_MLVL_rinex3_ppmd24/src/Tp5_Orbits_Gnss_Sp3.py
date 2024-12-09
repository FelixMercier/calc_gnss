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
    
    (orb,nl) = sp3.getSp3(const,prn)
    orb[:,1:4] *=1000
    orb[:,4] /= 1e6
    
    #trouver ti
    id=0
    for i in range(nl):
        if orb[i, 0] >= mjd: 
            id = i-1
            break

    #extrait des m dates
    mp1s2 = int((ordre+1)/2)
    

    #effet de bords 
    if id + mp1s2 -1 > len(orb):
        pos = orb[-10:, 1:]
        t_list = orb[-10:, 0]
    elif id - mp1s2 < 0:
        pos = orb[:10, 1:]
        t_list = orb[:10, 0]
    else:
        t_list = orb[id - mp1s2 : id + mp1s2, 0]
        pos = orb[id - mp1s2 : id + mp1s2, 1:]
    
    
    #interpolation 
    Pos = Pol_Lagrange(ordre, mjd, t_list, pos)

    X_ECEF, Y_ECEF, Z_ECEF, dte = Pos[0], Pos[1], Pos[2], Pos[3]

    return X_ECEF, Y_ECEF, Z_ECEF, dte

def Lagrange(j, m, t_list, t):
    Lj=1
    for k in range(m):
        if  k!= j:
            Lj *= ((t-t_list[k]) / (t_list[j] - t_list[k]))
    return Lj

def Pol_Lagrange(m, t, t_list, Y):
    y=0
    for i in range(m):
        y += Lagrange(i, m, t_list, t) * Y[i, :]
    return y



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
    
    (X1,Y1,Z1,dte1)	= pos_sat_sp3(mysp3, constellation,prn,t.mjd-0.001/86400,ordre)
    (X2,Y2,Z2,dte2)	= pos_sat_sp3(mysp3, constellation,prn,t.mjd+0.001/86400,ordre)
    V = np.array([[X2-X1],[Y2-Y1],[Z2-Z1]]) / 0.002
    print("V = %.3f m/s" % np.linalg.norm(V))

    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))

