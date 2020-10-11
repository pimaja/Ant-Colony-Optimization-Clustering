# -*- coding: utf-8 -*-
"""
Created on Sun Mar 15 14:53:23 2020

@author: Lucija
"""

# Importing the libraries
import numpy as np
import pandas as pd

# Importing the dataset
dataset = pd.read_csv('processed.cleveland.data', header=None, delimiter=',')
dataset[11] = pd.to_numeric(dataset[11], errors='coerce') 
dataset[12] = pd.to_numeric(dataset[12], errors='coerce') 
X = dataset.iloc[:, 0:13].values
y = dataset.iloc[:, 13].values