#!/usr/bin/env python
# coding: utf-8

#==================================================================================
# Striped and Spotted Binary Image Sweeping Filtration
# asociated with the manuscript titled:
#   "Quantifying Topological Features and Irregularities in Zebrafish 
#   Patterns Using the Sweeping-plane Filtration"
#
#   Version 1.0
#   (C) 2025/09/30
#    Nour Khoudari, John T. Nardini, Alexandria Volkening
#
#   A function that loads binary images of zebrafish skin pattern generated from data
#   simulations of an agent-based model: [A. Volkening and B. Sandstede, Iridophores
#   as a source of robustness in zebrafish stripes and variability in Danio patterns,
#   Nat. Commun., 2018] to compute persistence homology of sweeping filtrations in
#   all directions: top-to-bottom(TB), bottom-to-top(BT), left-to-right(LR), and
#   right-to-left(RL), with periodicity in x-direction in TB and BT only.
#
#   We use two data set of 1000 images each, the user can choose between
#   dataset = 1 for wild-type zebrafish skin patterns
#   dataset = 2 for pfeffer mutant zebrafish skin patterns
#   The binary images are generated for pixel size stepX = 80 \mum using
#   generate_binary_image.m and the simulation data used can be found on Figshare:
#   M. McGuirl, A. Volkening, and B. Sandstede, Zebrafish simulation data, 2020.
#   
#   Outputs:
#   .txt files containing persistence homology of sweeping in all four
#   directions across binary images in the data set
#==================================================================================



import numpy as np
import imageio.v2
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas as pd
import glob, pdb, os, sys
from ripser import ripser
from persim import plot_diagrams
import gudhi as gd
from methodsNardini_periodic import param_sweep,Persist_im, betti_curve


def main():
    descriptor = "applies sweeping and flooding to zebrafish images binned like an on-lattice model"

    # dataset represents the data set to be loaded, dataset = 1 referes to wild-type, dataset = 2 referes to pfeffer
    dataset = 1 # user should choose 1 or 2
        
    # Establishing folder to draw from and standard file name of an image
    if dataset == 1: 
        folderName = "./WT_binnedImages/" 
        dataForAnalysis = "centered_wt_day45_sim" 
        
    elif dataset == 2:
        folderName = "./WT_spottedImages/"
        dataForAnalysis = "pfef_day45_sim"
        
    else:
        raise ValueError("Unknown dataset number")
    
    
    # numFiles represents the number of files to loop over in the data set 
    numFiles = 1001 #files+1

    # Looping through our data
    for dataj in range(numFiles-1):

        for datai in range(80,85,5):#range(30,155,5):
        
            # Gathering our png image file 
            inFile = folderName + dataForAnalysis + str(dataj+1) + '_step' + str(datai) +'.png'
            binaryIm = imageio.v2.imread(inFile)
        
            # Note that the image gets flipped on its side in this process (so xNum is the image height and yNum is its width)
            xNum,yNum = binaryIm.shape         # size of image in pixels
            xSteps = np.linspace(0,1,xNum)     # vector of steps from 0 to 1 (inclusive of both endpoints, with step width 0/xNum)
            ySteps = np.linspace(0,1,yNum)     # similar vector of steps from 0 to 1 a
            Y,X = np.meshgrid(ySteps,xSteps)   # used for slicing through the image
        
            # Orients refer to the direction in which we are filtering FROM: this means 'left' describes
            # filtering FROM the left to the right, and 'top' means filtering down from the TOP, for example.
            # This represents a change from the style of Nardini et al's paper.
            orients = ['left','right','top','bottom']
        
            # For-loop processing the data by filtering in each direction in turn
            for orient in orients:
                if orient == "top":
                    plane_dir = 'less'
                    indep_var = X
                    periodic_y = True
                    r1b = xNum
                elif orient == "bottom":
                    plane_dir = 'greater'
                    indep_var = X
                    periodic_y = True
                    r1b = xNum
                elif orient == "right":
                    plane_dir = 'greater'
                    indep_var = Y
                    periodic_y = False
                    r1b = yNum
                elif orient == "left":
                    plane_dir = 'less'
                    indep_var = Y
                    periodic_y = False
                    r1b = yNum

                # Computing persistence diagram for the given data and filter direction
                diag = param_sweep(binaryIm,indep_var,iter_num=r1b, plane_dir=plane_dir,filename=folderName + dataForAnalysis + '_plane_' + orient + '_bd_data_' + str(dataj+1),periodic_y = periodic_y) 

                # Writing and saving the barcodes in dimension 0
                with open(folderName + dataForAnalysis + '_plane_' + orient + '_dim0_' + str(dataj+1) + '_step' + str(datai) + '.txt', 'w') as f:
                    for feature in diag:
                        if feature[0] == 0:
                            line = f"{feature[1][0]},{feature[1][1]}"
                            f.write(line)
                            f.write('\n')
            
                # Repeating the same thing for the barcode in dimension 1
                with open(folderName + dataForAnalysis + '_plane_' + orient + '_dim1_' + str(dataj+1) + '_step' + str(datai) + '.txt', 'w') as f:
                    for feature in diag:
                        if feature[0] == 1:
                            line = f"{feature[1][0]},{feature[1][1]}"
                            f.write(line)
                            f.write('\n')

                #write_PD_text(diag,folderName + dataForAnalysis + '_plane_' + orient + '_dim0_' + str(dataj+1)
            
                # Computing and saving Betti curves
                b0,b1,r = betti_curve(diag, r0 = 0, r1 = r1b, rN = r1b+1,filename_save=folderName + dataForAnalysis + '_plane_' + orient + '_Betti_' + str(dataj+1)+ '_step' + str(datai))


        

if __name__=="__main__":
    main()


