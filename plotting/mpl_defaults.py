#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 12:32:32 2018

@author: alexisklimasewski
"""

import matplotlib.pyplot as plt

import matplotlib as mpl

def plot_defaults():
    plt.style.use("classic")
    mpl.rcParams['font.size'] = 20
    mpl.rcParams['xtick.labelsize'] = 20
    mpl.rcParams['ytick.labelsize'] = 20
    mpl.rcParams['xtick.major.size'] = 7
    mpl.rcParams['ytick.major.size'] = 7
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['ytick.major.width'] = 2
    mpl.rcParams['figure.subplot.hspace'] = 0.3
    mpl.rcParams['figure.subplot.left'] = 0.08
    mpl.rcParams['figure.subplot.right'] = 0.99
    mpl.rcParams['figure.subplot.top'] = 0.95
    mpl.rcParams['figure.subplot.bottom'] = 0.1