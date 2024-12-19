# -*- coding: utf-8 -*-
"""
@author: BEILIN Jacques - IGN/ENSG
@date  : %(date)s
"""

import numpy as np
import math
import os
import json
import matplotlib.pyplot as plt
from scipy.stats import norm

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
    for sat in sats:
        PosSat.append([sat["Xs"], sat["Ys"], sat["Zs"]])
        Dobs.append(sat["Dist"])
        
    PosSat = np.array(PosSat)
    Dobs = np.array(Dobs)
    X0 = np.array(data["X0"])
    return PosSat, Dobs, X0

def create_matB(Dobs, PosSat, X0):
    # if X0.shape == (4,):
    #     X0 = X0.reshape(-1,1)
    # print("x0 : ", X0[1][0])
    Pcorr = Dobs.reshape(-1,1)
    Xsat = PosSat[:, 0].reshape(-1,1)
    Ysat = PosSat[:, 1].reshape(-1,1)
    Zsat = PosSat[:, 2].reshape(-1,1)
    # print("Xsat : ", Xsat)
    Pcorr0 = np.sqrt( (X0[0][0] - Xsat)**2 + (X0[1][0] - Ysat)**2 + (X0[2][0] - Zsat)**2) + X0[3][0]
    B = Pcorr - Pcorr0
    return B

def create_matA(Dobs, PosSat, X0):
    # if X0.shape == (4,):
    #     X0 = X0.reshape(-1,1)
    Xsat = PosSat[:, 0]
    Ysat = PosSat[:, 1]
    Zsat = PosSat[:, 2]
    A = np.ones((Dobs.shape[0], 4))
    A[:, 0] = (X0[0][0] - Xsat) / np.sqrt((X0[0][0] - Xsat)**2 + (X0[1][0] - Ysat)**2 + (X0[2][0] - Zsat)**2)
    A[:, 1] = (X0[1][0] - Ysat) / np.sqrt((X0[0][0] - Xsat)**2 + (X0[1][0] - Ysat)**2 + (X0[2][0] - Zsat)**2)
    A[:, 2] = (X0[2][0] - Zsat) / np.sqrt((X0[0][0] - Xsat)**2 + (X0[1][0] - Ysat)**2 + (X0[2][0] - Zsat)**2)
    return A

def create_matP(n, err):
    P = np.eye(n)*err
    return P




def Gauss_Newton(Dobs, PosSat, X0, err):
    
    X0 = X0.reshape(-1,1)
    # print("X0 or boucle : ", X0.shape)
    matA = create_matA(Dobs, PosSat, X0)
    # print('matA : ', matA)
    matB = create_matB(Dobs, PosSat, X0)
    # print('matB : ', matB)
    P = create_matP(Dobs.shape[0],err)
    # print('matP : ', P)
    dX = np.linalg.inv(matA.T@P@matA)@matA.T@P@matB 
    # print('dX : ', dX)
    X_chap = X0 + np.linalg.inv(matA.T@P@matA)@matA.T@P@matB
    # print('X_chap or boucle : ', X_chap)
    Vchap =  matB - matA@dX
    # print("Vchap or boucle : ", Vchap)
    n=matA.shape[0] 
    p=matA.shape[1] 
    sigma02_old = 0
    sigma02=(Vchap.T@P@Vchap)[0,0]/float(n-p)
    # print('sigma02 : ', sigma02)
    varVchap = sigma02*(np.linalg.inv(P)-matA@np.linalg.inv(matA.T@P@matA)@matA.T)
    Vnor=np.linalg.inv(np.sqrt(np.diag(np.diag(varVchap))))@Vchap
    cpt = 0
    print('itération : ', cpt)
    # print('X_chap : ', X_chap)
    while np.abs(sigma02-sigma02_old)>10e-6:
        X0 = X_chap
        # print("X0 dans boucle: ", X0.shape)
        sigma02_old = sigma02
        matB = create_matB(Dobs, PosSat, X0)
        matA = create_matA(Dobs, PosSat, X0)
        n=matA.shape[0] 
        p=matA.shape[1] 
        X_chap = X0 + np.linalg.inv(matA.T@P@matA)@matA.T@P@matB
        # print("N : ", matA.T@P@matA)
        # print("K : ", matA.T@P@matB)
        # print('X_chap : ', X_chap.shape)
        dX = np.linalg.inv(matA.T@P@matA)@matA.T@P@matB
        print('dX : ', dX)
        Vchap =  matB - matA@dX
        # print("Vchap dans boucle : ", Vchap)
        sigma02=(Vchap.T@P@Vchap)[0,0]/float(n-p)
        print('sigma02 : ', sigma02)
        varVchap = sigma02*(np.linalg.inv(P)-matA@np.linalg.inv(matA.T@P@matA)@matA.T)
        # print('varVchap : ', varVchap)
        Vnor=np.linalg.inv(np.sqrt(np.diag(np.diag(varVchap))))@Vchap
        # print('Vnor : ', Vnor)
        cpt += 1
        print('itération : ', cpt)
    SigmaX = sigma02 * np.linalg.inv(matA.T@P@matA)
    return X_chap, Vchap, sigma02, Vnor, SigmaX


        




def estimation(PosSat,Dobs,X0,err):
    """ Calcul d'une trilatération simple avec un offset
    correspondant à l'erreur d'horloge du récepteur """
    
    # print("PosSat : ", PosSat.shape)
    # print("Dobs : ", Dobs.shape)
    # print("X0 : ", X0.shape)

    SigmaX = []
    V = []
    X = 0
    Y = 0
    Z = 0
    cdtr = 0
    sigma0_2 = 0

    X_chap, Vchap, sigma0_2, Vnor, SigmaX = Gauss_Newton(Dobs, PosSat, X0, err)
    X = X_chap[0][0]
    Y = X_chap[1][0]
    Z = X_chap[2][0]
    cdtr = X_chap[3][0]
    V = Vchap

    return X,Y,Z,cdtr,sigma0_2,V,SigmaX


if __name__ == "__main__":

    tic = gpst.gpsdatetime()
 
    dir_data = '../data'
    PosSat, Dobs, X0 = loadJson(os.path.join(dir_data, "data_tp7.json"))
    print('Dobs  : ', Dobs.shape[0])
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

    # B = create_matB(Dobs, PosSat, X0)

    X,Y,Z,cdtr,sigma0_2,V, SigmaX = estimation(PosSat,Dobs,X0,err=0.25)
    print("\nSolution TP\nX = %13.3f m\nY = %13.3f m\nZ = %13.3f m\nc*dtr = %.9f m" % (X,Y,Z,cdtr))
    print("sigma0_2 = %.3f" % (sigma0_2))
    print("V = ",V,"\nSigmaX = ",SigmaX)

    X_chap, Vchap, sigma02, Vnor, SigmaX = Gauss_Newton(Dobs, PosSat, X0, err=0.25)
    print('X_chap : ', X_chap)
    # print('Vchap : ', Vchap)
    print('sigma02 : ', sigma02)

    # Histogramme pour Vnor
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.hist(Vnor, bins=30, density=True, alpha=0.6, color='b', label='Vnor Histogram')
    # Calcul de la courbe gaussienne
    mu, std = norm.fit(Vnor)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2, label='Gaussian Fit')
    plt.title('Histogram of Vnor with Gaussian Fit : residual normalized')
    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.legend()

    # Graphique pour Vchap
    plt.subplot(1, 2, 2)
    plt.scatter(range(len(Vchap)), Vchap, color='r', s=5, label='Vchap Points : residual')  # s=5 pour la taille des points
    plt.title('Vchap Values : residual')
    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.legend()

    plt.tight_layout()
    plt.show()

    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))





