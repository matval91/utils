#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 11:02:10 2019

@author: vallar
"""
import numpy as np

def cyl2euc(Arho, Aphi, phi):
    Ax = np.cos(phi)*Arho+np.sin(phi)*Aphi
    Ay = -1.*np.sin(phi)*Arho+np.cos(phi)*Aphi
    return Ax, Ay