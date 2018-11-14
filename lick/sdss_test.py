# -*- coding: utf-8 -*-
""" 

Created on 13/11/18

Author : Carlos Eduardo Barbosa

Testing Lick indices using SDSS data

"""
from __future__ import print_function, division

import os

import numpy as np
from astropy.table import Table, hstack
import matplotlib.pyplot as plt

from lick import Lick
from convolve import broad2lick

if __name__ == "__main__":
    work_dir = "/home/kadu/Dropbox/ESPECTROS/ORIGINALS"
    os.chdir(work_dir)
    spectra = [_ for _ in os.listdir(work_dir) if _.endswith(".fits")]
    # Reading ascii file
    table_vel = Table.read("lick-fot-SN-cmodel.dat",
                           format="ascii")
    vels = Table([table_vel["col1"], table_vel["col68"]], names=["filename",
                                                                  "vel"])
    bandsfile = "bands.txt"
    bands = np.loadtxt(bandsfile, usecols=(2, 3, 4, 5, 6, 7,))
    bandsnames = np.loadtxt(bandsfile, usecols=(0,), dtype=str )
    results = []
    for spec in spectra:
        print(spec)
        table = Table.read(spec, hdu=1)
        wave = np.power(10, table["loglam"]).data
        flux = table["flux"].data
        idx = np.where(vels["filename"] == spec)[0][0]
        v = vels["vel"][idx]
        flux_lick = broad2lick(wave, flux, 2.75, v)
        ll = Lick(wave, flux_lick, bands, vel=v)
        ll.classic_integration()
        results.append(ll.Ia)
    results = np.array(results)
    results_table = Table(results, names=bandsnames)
    newtable = hstack([Table([spectra], names=["spec"]), results_table])
    newtable.write("results.txt", format="ascii", overwrite=True)


