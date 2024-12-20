# -*- coding: utf-8 -*-
"""
@author: BEILIN Jacques - IGN/ENSG
@date  : %(date)s
"""

import numpy as np
import math
import os
import json

import gpsdatetime as gpst


""" Import du module """
# version locale
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import gnss_process as proc

# version installee
import gnsstoolbox.gnss_process as proc

def loadJson(jsonFile):
    try:
        with open(jsonFile, "r") as read_file:
            data = json.load(read_file)
    except Exception as err:
        print(err)
        return
        
    sats = data["Sat"]
    PosSat = []
    Dobs = []
    ElevSat = []
    SatIndex = []
    for sat in sats:
        PosSat.append([sat["Xs"], sat["Ys"], sat["Zs"]])
        Dobs.append(sat["Dist"])
        ElevSat.append(sat["ElevSat"])
        SatIndex.append(sat["SatIndex"])
        
    X0 = np.array(data["X0"])
    return np.array(PosSat), np.array(Dobs), X0, np.array(ElevSat), np.array(SatIndex)
        

def estimation(PosSat, Dobs, ElevSat, SatIndex, X0):
    """ Calcul d'une trilatération simple avec un offset
    correspondant à l'erreur d'horloge du récepteur """
    
    print(PosSat)
    print(Dobs)
    print(X0)

    SigmaX = []
    V = []
    X = 0
    Y = 0
    Z = 0
    cdtr = 0
    sigma0_2 = 0
    return X,Y,Z,cdtr,sigma0_2,V,SigmaX


if __name__ == "__main__":

    tic = gpst.gpsdatetime()
 
    dir_data = '../data'
    PosSat, Dobs, X0, ElevSat, SatIndex = loadJson(os.path.join(dir_data, "data_tp9.json"))


    X,Y,Z,cdtr,cGGTP, cGPGL, sigma0_2,V,SigmaX = proc.trilatGnssPonderationElev(PosSat,Dobs,X0,satIndex,ElevSat)
    print("\nSolution pygnssToolbox\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\nc*dtr = %.9f m" % (X,Y,Z,cdtr))
    print("s0 = ", sigma0_2,"\nV = ", V)
    print("SigmaX", SigmaX)

    X,Y,Z,cdtr,sigma0_2,V,SigmaX = estimation(PosSat, Dobs, ElevSat, SatIndex, X0)
    print("\nSolution TP\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\nc*dtr = %.9f m" % (X,Y,Z,cdtr))
    print("sigma0_2 = %.3f" % (sigma0_2))
    print("V = ",V,"\nSigmaX = ",SigmaX)

    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))
