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

def Alpha(X0, PosSat):
    x0, y0, z0 = X0[0], X0[1], X0[2]
    xs, ys, zs = PosSat[:, 0], PosSat[:, 1], PosSat[:, 2]
    return np.sqrt((x0-xs)**2 + (y0-ys)**2 + (z0-zs)**2 )

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
    for sat in sats:
        PosSat.append([sat["Xs"], sat["Ys"], sat["Zs"]])
        Dobs.append(sat["Dist"])
        
    PosSat = np.array(PosSat)
    Dobs = np.array(Dobs)
    X0 = np.array(data["X0"])
    return PosSat, Dobs, X0
        

def estimation(PosSat,Dobs,X0):
    """ Calcul d'une trilatération simple avec un offset
    correspondant à l'erreur d'horloge du récepteur """
    sigma0_2 = 1
    Xchap = np.zeros((4, 1))
    
    nitermax = 20
    iter = 0
    while sigma0_2 > 1e-6:
        iter += 1
        
        X0 += Xchap
        alpha = Alpha(X0, PosSat).reshape(-1, 1)
        B = Dobs - alpha - X0[-1]
        
        jx = X0[0] - PosSat[:, 0].reshape(-1, 1)
        jx = jx * 1/alpha
        jy = X0[1] - PosSat[:, 1].reshape(-1, 1)
        jy = jy * 1/alpha
        jz = X0[2] - PosSat[:, 2].reshape(-1, 1)
        jz = jz * 1/alpha
        jdt = np.ones((len(Dobs), 1))
        A = np.hstack((jx, jy, jz, jdt))

        sig = 2
        sigmaB = sig * np.eye(len(Dobs))
        P = np.linalg.inv(sigmaB)
        
        N = A.T@P@A
        K = A.T@P@B
        Xchap = np.linalg.inv(N) @ K

        
        V = B - A@Xchap
        
        sigma0_2 = (V.T @ P @ V) / (len(Dobs) - len(Xchap))
        
        
        if iter > nitermax: break

    X0 += Xchap
    X, Y, Z, cdtr = X0[0], X0[1], X0[2], X0[3]

    
    return X,Y,Z,cdtr,sigma0_2,V


if __name__ == "__main__":

    tic = gpst.gpsdatetime()
 
    dir_data = '../data'
    PosSat, Dobs, X0 = loadJson(os.path.join(dir_data, "data_tp7.json"))
    Dobs = Dobs.reshape(-1, 1)
    X0 = X0.reshape(-1, 1)
    X0 = X0.astype('float64')
#    PosSat = np.genfromtxt("possat.txt",skip_header=1,delimiter=",")
#    Dobs = np.genfromtxt("dobs.txt",skip_header=1,delimiter=",")
#    
#    L = []
#    for i in range(len(Dobs)):
#        S = PosSat[i,:]
#        sat = {"Xs" : S[0], "Ys" : S[1], "Zs" : S[2], "Dist": Dobs[i]}
#        L += [sat]
#    
#    with open("data_tp7.json", "wt") as f:
#        s = json.dumps({"Sat": L}, indent=4)
#        f.write(s)

#
    X,Y,Z,cdtr,sigma0_2,V,SigmaX = proc.trilatGps(PosSat,Dobs,X0)
    print("\nSolution pygnssToolbox\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\nc*dtr = %.9f m" % (X,Y,Z,cdtr))
    print("sigma0_2 = ", sigma0_2,"\nV = ", V)
    print("SigmaX", SigmaX)

    X,Y,Z,cdtr,sigma0_2,V = estimation(PosSat,Dobs,X0)
    print("\nSolution TP\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\nc*dtr = %.9f m" % (X,Y,Z,cdtr))
    print("sigma0_2 = %.3f" % (sigma0_2))
    # print("V = ",V,"\nSigmaX = ",SigmaX)

    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))
