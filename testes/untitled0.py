# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 14:21:20 2015

@author: Danielbraga
"""
from __future__ import print_function
import os
for file in os.listdir("../res/blast_without_note"):
    if file.endswith(".xml"):
        print(file)
