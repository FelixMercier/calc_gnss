
#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

# GNSS observation class
# Beilin Jacques
# 2014-06-02


import re
import math
import numpy as np
import copy
import os
import gpsdatetime as gpst
import Tp5_Orbits_Gnss_Sp3 as tp5
import Tp3_Orbits_Gps_Gal_Nav as tp3
import TP7clem as tp7
# import Tp7_Trilat as tp7

#sys.path.append('../TD07_trilat_GPS')
#import Tp7_Trilat as TP7

""" Import du module orbits """
# version locale
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import orbits
#import rinex_o as rx
#import gnss_process as proc

# version installee
import gnsstoolbox.orbits as orbits
import gnsstoolbox.rinex_o as rx
import gnsstoolbox.gnss_const as const
from gnsstoolbox.gnss_const import d2r
import gnsstoolbox.gnsstools as tool
import gnsstoolbox.gnss_process as proc


class gnss_process_TP():
    """GNSS proccess class"""


    def __init__(self):

        """Options par defaut"""
        self.process="spp" #"DGNSS","phase"
        self.freq="C1"
        self.cut_off=3*const.d2r
        self.X0 = np.zeros((3,1))
        self.iono = 'none'
        self.nav = 'brdc'
        self.const= 'GRE'
        self.type="o"
        self.constraint=0 # contraintes sur la position, utilisé pour solution DGNSS

        self.nb_sat = 0
        self.nb_GPS = 0
        self.nb_GLO = 0
        self.nb_GAL = 0

        self.nb_GPS_calc = 0
        self.nb_GLO_calc = 0
        self.nb_GAL_calc = 0

    def spp(self,epoch,brdc):
        """Process epoch using spp"""

        Dobs = []
        PosSat=[]

        self.nb_sat=len(epoch.satellites)

        c = 299792458.0
        mu = 3.986005e14
        F = -4.442807633e-10
        T = 23.934 * 3600  #jour sidéral en secondes

#        print(epoch.satellites[0].__dict__.keys())

        mjd = epoch.tgps.mjd
        for S in epoch.satellites:

            """ Test si le satellite S est bien dans une constellation qu'on souhaite utiliser """
            if not(re.search(S.const,self.const)):
                continue

            observable = 'C1C'
#            print(S.obs)
            S.PR = S.obs.get(observable)

            """ Test si l'observation est bien cohérente """
            if S.PR < 15e6:
                continue

            """ Dispose-t-on d'un message de navigation pour ce satellite/epoch ? """
            Eph = brdc.getEphemeris(S.const,S.PRN,mjd)
            if not(hasattr(Eph,'mjd')):
                print("No orbit for %1s%2s" % (S.const,S.PRN))
                continue

#            print(S.const,S.PRN)
#            print("mjd = ",epoch.tgps.mjd,observable, "=", S.obs.get(observable))

            """ A ce stade on est certain de disposer d'une observation et d'un
            message de navigation utilisable pour le satellite S à l'époque qui nous interesse"""

            """ Calcul de la date de réception (mjd) """
            S.tr = epoch.tgps.mjd

            """ Calcul du temps de vol et de la date d'émission """
            S.tv = S.PR / c            #en sec
            S.tv = (S.tv / 86400)      #en jour
            S.te = S.tr - S.tv         #en jour 
            

            if self.nav == 'brdc':
                """ Correction de la derive d'horloge du satellite """
                S.dte = Eph.alpha0 + Eph.alpha1 * ((S.te - Eph.TOC)*86400) + Eph.alpha2*((S.te-Eph.TOC)*86400)**2

                S.te -= S.dte/86400
                S.PR += c*S.dte
                """ Calcul de l'effet relativiste"""
                t_eph = gpst.gpsdatetime(mjd=S.te)
                t = gpst.gpsdatetime(mjd=mjd)
                
                delta_t = t - t_eph
                
                n = np.sqrt(mu / (Eph.sqrt_a**6) ) + Eph.delta_n

                M = Eph.M0 + n*delta_t
                
                E0=M
                E = E0+1
                R = 30e6
                
                while np.abs(R*E - R*E0) > 1e-3:
                    E0 = E
                    E = M + Eph.e * np.sin(E0)
                S.dtrelat = F * Eph.e * Eph.sqrt_a * np.sin(E)
                
                S.te -= S.dtrelat/86400
                S.PR += c*S.dtrelat

            if self.nav == 'sp3':
                """ Correction de la derive d'horloge du satellite """
                S.X, S.Y, S.Z, S.dte = tp5.pos_sat_sp3(brdc, S.const, S.PRN, S.te, 9)
                
                S.te -= S.dte
                S.PR += c*S.dte
                
                """ Calcul de l'effet relativiste """                
                (X1,Y1,Z1,dte1)	= tp5.pos_sat_sp3(brdc, S.const,S.PRN,S.te.mjd-0.001/86400, 9)
                (X2,Y2,Z2,dte2)	= tp5.pos_sat_sp3(brdc, S.const,S.PRN,S.te.mjd+0.001/86400, 9)
                V = np.array([[X2-X1],[Y2-Y1],[Z2-Z1]]) / 0.002
                
                S.dtrelat = -2 * (S.X * V[0,0] + S.Y * V[1,0] + S.Z*V[2,0]) / (c**2)
                
                
                S.te -= S.dtrelat/86400
                S.PR += c*S.dtrelat
                
                
            """ Calcul de la position des satellites a te """
            S.X, S.Y, S.Z, dte = tp3.tp3_pos_sat_brdc(brdc, Eph.const, Eph.PRN, S.te)
          

            
            """ Et pourtant elle tourne """
            S.tv = (S.tr - S.te) * 86400
            
            alpha = (2*np.pi * S.tv/T)
            
            S.X, S.Y, S.Z = tool.toolRotZ(S.X, S.Y, S.Z, -alpha)
            
            print("PosSat", S.X, S.Y, S.Z)
            
            """ Sauvegarde des donnees """
            Prcorr = S.PR
            Xs = S.X
            Ys = S.Y
            Zs = S.Z
            Dobs.append(Prcorr)
            PosSat.append([Xs,Ys,Zs])

        Dobs = np.array(Dobs).reshape(-1, 1)
        PosSat = np.array(PosSat)
        # print("PosSat = ", PosSat)
        print("Dobs = ", Dobs)
        X0 = np.zeros((4, 1))
        
        X,Y,Z,cdtr,sigma0_2,V,SigmaX = tp7.estimation(PosSat, Dobs, X0, err=1)
        # X,Y,Z,cdtr,sigma0_2,V, iter = tp7.estimation(PosSat, Dobs, X0)
        print(X, Y, Z, cdtr)
        print("écarts : ", X - 4201603.471, Y - 189865.522, Z - 4779084.965, cdtr -28.453 )
        
        
        return X, Y, Z


def main():

    """ lecture du fichier BRDC """
    Orb = orbits.orbit()
    dir_orb = '../data'
    Orb.loadRinexN(os.path.join(dir_orb, 'BRDC00IGN_R_20212400000_01D_MN.rnx'))

    """ lecture du fichier SP3 """
    sp3 = orbits.orbit()
    sp3.loadSp3(os.path.join(dir_orb,'GFZ0OPSULT_20212400600_02D_05M_ORB.SP3'))

    """ Definition de l'instant pour lequel on cherche une position """
    t = gpst.gpsdatetime(yyyy=2021,doy=240,dsec=5400)
    print(t)

    rnx = rx.rinex_o()
    filename = os.path.join(dir_orb,'MLVL00FRA_R_20212400000_01D_30S_MO.21o')
    ret = rnx.loadRinexO(filename)

    if ret<0:
        print(ret)
        return

    Ep = rnx.getEpochByMjd(t.mjd)
    
    """ Calcul avec pygnsstoolbox """
    spp2 = proc.gnss_process()
    spp2.const='G'

    print("Calcul sur les orbites brdc")
    Ep2 = spp2.spp(Ep,Orb)
    print("X = %.2fm Y = %.2fm Z = %.2fm " % (Ep2.X, Ep2.Y, Ep2.Z))
#
    print("PRN           %-16s %-16s %-16s %-16s" % ('Xs', 'Ys', 'Zs', 'PR'))
    for s in Ep2.satellites:
        print("%1s%02d %16.3f %16.3f %16.3f %16.3f" % (s.const,s.PRN, s.Xs, s.Ys, s.Zs, s.PR))
    print("V",Ep2.V)

    print("Calcul sur les orbites sp3")
    Ep3 = spp2.spp(Ep,sp3)
    print("X = %.2fm Y = %.2fm Z = %.2fm " % (Ep3.X, Ep3.Y, Ep3.Z))

    print("PRN           %-16s %-16s %-16s %-16s" % ('Xs', 'Ys', 'Zs', 'PR'))
    for s in Ep3.satellites:
        print("%1s%02d %16.3f %16.3f %16.3f %16.3f" % (s.const,s.PRN, s.Xs, s.Ys, s.Zs, s.PR))
    print("V",Ep3.V)
#    print(Ep2.__dict__)

#    [E,N,U]= tool.tool_cartloc_GRS80(float(spp2.X0[0]),float(spp2.X0[1]),float(spp2.X0[2]),E2.X,E2.Y,E2.Z)
#    print("dE = %.2fm dN = %.2fm dU = %.2fm " % (E,N,U))


    """ Calcul développé au cours du TP """

    spp1 = gnss_process_TP()
    spp1.const='G'
    spp1.constraint=0
    spp1.cut_off = 10 * d2r
    spp1.X0[0]=rnx.headers[0].X
    spp1.X0[1]=rnx.headers[0].Y
    spp1.X0[2]=rnx.headers[0].Z
    spp1.spp(Ep,Orb)
    print()
    

if __name__ == '__main__':


    tic = gpst.gpsdatetime()
    main()
    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))