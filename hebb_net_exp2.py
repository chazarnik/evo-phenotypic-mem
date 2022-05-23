"""This script replicates experiments with the Hebb learning rule from 
The Evolution of Phenotypic Correlations and 'Developmental Memory' - 
Watson et al, Figures 1C and 2C"""

import numpy as np
import matplotlib.pyplot as plt
import os, sys

from generate_plots import plot_B_matrix_instance

INDIV_SIZE = 8

TARGET_1 = np.array([1,1,-1,-1,-1, 1,-1,1]) # Exp1
TARGET_2 = np.array([1,1,-1,-1,-1,1,-1,1]) # Exp2
TARGET_3 = np.array([1,-1,1,-1,1,-1,-1,-1]) # Exp3

def hebb_net(target, iterations=40, l_rate=0.05):
    """Implement a network working with Hebb's learning rule.
    This function returns the weights of the trained network."""
    w_matrix = np.zeros((target.shape[0], target.shape[0]))
    w_temp = np.zeros((w_matrix.shape))
    for i in range(iterations):
        for j in range(target.shape[0]):
            w_matrix[j,:] = w_matrix[j,:] + l_rate*(target[j]*target)

    return w_matrix

w_target1 = hebb_net(TARGET_1)
plot_B_matrix_instance(w_target1, "Hebb Network Weights Matrix \n Single Environment")

w_target2 = hebb_net(TARGET_2)
w_target3 = hebb_net(TARGET_3)
w_exp2 = w_target2+w_target3
plot_B_matrix_instance(w_exp2, "Hebb Network Weights Matrix \n Two-Target Environment")