
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

""" Import du module orbits """
# version locale
#import sys
#sys.path.append("/home/beilin/progs/Prog_scientifique/Packages/yagnss_toolbox/python/src/pygnsstoolbox/gnsstoolbox")
#import orbits
#import rinex_o as rx
#import gnss_process as proc
#import gnsstools as tool

# version installee
import gnsstoolbox.orbits as orbits
import gnsstoolbox.rinex_o as rx
import gnsstoolbox.gnss_const as const
from gnsstoolbox.gnss_const import d2r
import gnsstoolbox.gnsstools as tool
import gnsstoolbox.gnss_process as proc

class gnss_process():
    """GNSS proccess class"""


    def __init__(self):

        """Options par defaut"""
        self.process="spp" #"DGNSS","phase"
        self.freq="C1"
        self.cut_off=3*const.d2r
        self.X0 = np.zeros((3,1))
        self.iono = 'none'
        self.nav = 'sp3'
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

        self.nb_sat=len(epoch.satellites)

        c = 299792458.0

        print(epoch.satellites[0].__dict__.keys())

        mjd = epoch.tgps.mjd
        for S in epoch.satellites:

            """ Test si le satellite S est bien dans une constellation qu'on souhaite utiliser """
            if not(re.search(S.const,self.const)):
                continue

            observable = 'C1'
            S.PR = S.getObs(observable)

            """ Test si l'observation est bien cohérente """
            if S.PR<15e6:
                continue

            """ Dispose-t-on d'un message de navigation pour ce satellite/epoch ? """
            Eph = brdc.getEphemeris(S.const,S.PRN,mjd)
            if not(hasattr(Eph,'mjd')):
                continue

            print(S.const,S.PRN)
            print("mjd = ",epoch.tgps.mjd,observable, "=", S.getObs(observable))

            """ A ce stage on est certain de disposer d'une observation et d'un
            message de navigation utilisable pour le satellite S à l'époque qui nous interesse"""

            """ Calcul de la date de réception (seconde dans la semaine) """
            S.tr = epoch.tgps.mjd

            """ Calcul du temps de vol et de la date d'émission """

            """ Correction de la derive d'horloge du satellite """

            """ Calcul de l'effet relativiste """

            """ Calcul de la position des satellites a te """

            """ Et pourtant elle tourne """

            """ Sauvegarde des donnees """


        return


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

    spp1 = gnss_process()
    spp1.const='G'
    spp1.constraint=0
    spp1.cut_off = 10 * const.d2r
    spp1.X0[0]=rnx.headers[0].X
    spp1.X0[1]=rnx.headers[0].Y
    spp1.X0[2]=rnx.headers[0].Z
    spp1.spp(Ep,Orb)

    spp2 = proc.gnss_process()
    spp2.const='G'
    print(spp2.X0, rnx.headers[0].X)
    spp2.X0[0][0] = rnx.headers[0].X
    spp2.X0[1][0] = rnx.headers[0].Y
    spp2.X0[2][0] = rnx.headers[0].Z
    Ep2 = spp2.spp(Ep,sp3)
    print("X = %.2fm Y = %.2fm Z = %.2fm " % (Ep2.X, Ep2.Y, Ep2.Z))
    print("V=",Ep2.V)


if __name__ == '__main__':


    tic = gpst.gpsdatetime()
    main()
    toc = gpst.gpsdatetime()
    print ('%.3f sec elapsed ' % (toc-tic))