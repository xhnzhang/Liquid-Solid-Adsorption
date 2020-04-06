#!/usr/bin/env python
# Calculate diffusion/sticking duration
# Xiaohong Zhang
# Getman Research Group
# Jun 26, 2018

import sys, os, argparse
import re, time
import numpy as np
import csv, ast
from collections import OrderedDict
import matplotlib.pyplot as plt
from matplotlib import rc_params

# Usage: put 2_analyze_data.py inside O_results/ folder
# python 2_analyze_data.py 
# output file: resTimeSepDist.csv, resTimeDist.csv, hopEventDict.csv, rateConstDict.csv
##################### manu input starts #######################
DUMP_FREQ = 10  # dumping frequency of coords 
TIME_STEP = 1.0 # [fs]
TIME_UNIT = DUMP_FREQ*TIME_STEP/1000.0  # convert [fs] time unit to [ps]
##################### manu input ends #########################



def runLengthEncoding(input_layer_array):
    """ encode duration time of molecule in each integer layer """
    length = len(input_layer_array)                    # input_layer_array is np array
    ia = np.array(input_layer_array[:,1],dtype='int')  # ignor time, only extract layer number
    
    if (length == 0):
        return dict(run_length = None,
                    positions = None,
                    layer_num = None)

    else:
        y = np.array(ia[1:] != ia[:-1])      # pairwise unequal (string safe)
        i = np.append(np.where(y), length-1)      # must include last element posi
        z = np.diff(np.append(-1, i))        # run lengths
        p = np.cumsum(np.append(0, z))[:-1]  # positions
        return dict(run_length = np.asarray(z),   # duration time of every crossing
                    positions = np.asarray(p),
                    layer_num = np.asarray(ia[i]))   


def calcResTimeDistribution(RLE_dict):
    """ calculate distribution of residence time in each integer layer, unit [10fs] """ 
    """ record leave each layer, rather than enter each layer """
    resTime_dict = OrderedDict()  # distribution of residence time in each layer
    resTime_separate_dict = OrderedDict() # distribution of residence time for forward and reverse, separately
    layer_num = RLE_dict['layer_num']   # 1D numpy array of the residence layer along the simulation    
    
    for i in range(len(layer_num)-1):
        layer_idx = layer_num[i]
        residence_time = RLE_dict['run_length'][i] 
        f_idx = 'B'+str(layer_idx-1)+'F'   # forward
        r_idx = 'B'+str(layer_idx)+'R' # reverse 
             
        """ separate res time distribution with forward and reverse directions """
        if layer_idx not in resTime_separate_dict:  # initialize distribution of a layer into sub-dict
            resTime_separate_dict[layer_idx] = {f_idx:[], r_idx:[]}  # i.e., L2, f_idx is B1F, r_idx is B2R
 
        if (layer_idx > layer_num[i+1]):    # forward hopping, to the surface
            resTime_separate_dict[layer_idx][f_idx].append(residence_time)
        else:   # reverse hopping, away from the surface
            resTime_separate_dict[layer_idx][r_idx].append(residence_time)
    """ assume final layer is reverse hopping to future layer """
    final_layer = layer_num[-1]
    final_idx = 'B'+str(final_layer)+'R' # reverse
    final_resTime = RLE_dict['run_length'][-1]
    resTime_separate_dict[final_layer][final_idx].append(final_resTime)
    

    """ lumped total res time """
    for key in resTime_separate_dict:   # loop over integer layer as key
        #print(key)
        if key not in resTime_dict:
            resTime_dict[key]=[]    # initialize distribution list for each layer key
        for subkey in resTime_separate_dict[key]:   # loop over forword or reverse as subkey
            resTime_dict[key].append(resTime_separate_dict[key][subkey])

    for key in resTime_dict:   # flat out the list of list
        resTime_dict[key] = [item for sublist in resTime_dict[key] for item in sublist]

    
#    print(resTime_separate_dict)
    resTime_separate_dict = OrderedDict(sorted(resTime_separate_dict.items())) 
    resTime_dict = OrderedDict(sorted(resTime_dict.items()))   # dict was converted to list, then back to dict
#    print(resTime_separate_dict)
    
#    for k,v in residence_time_distribution_dict.items():
#        print (k, v)

    # save to disk
    with open('resTimeSepDist.csv',"w") as f:
        w = csv.writer(f)
        for k,v in resTime_separate_dict.items():
            w.writerow([k,v])
            
    with open('resTimeDist.csv',"w") as f:
        w = csv.writer(f)
        for k,v in resTime_dict.items():
            w.writerow([k,v])    

    return resTime_separate_dict,resTime_dict




def calcRateConstant(RLE_dict, int_layer_count):
    """ calculate apparent rate constant for leaving each integer layer, e.g., k_B1F = hop_up_B1F / residence time in L2 """

    hopEvent_dict = OrderedDict() # store hop information into a dictionary
    rateConst_dict = OrderedDict()   # divided by total res time in each layer
    layer_num = RLE_dict['layer_num']   # 1D numpy array of the residence layer along the simulation    
    
    """ hopping event with forward and reverse directions """
    for i in range(len(layer_num)-1):
        layer_idx = layer_num[i]
#        print(layer_idx)

        f_idx = 'B'+str(layer_idx-1)+'F'   # forward
        r_idx = 'B'+str(layer_idx)+'R' # reverse 
        
        if layer_idx not in hopEvent_dict:  # initialize hopEvent_dict into sub-dict
            hopEvent_dict[layer_idx] = {f_idx:0, r_idx:0}  # i.e., L2, f_idx is B1F, r_idx is B2R
            
        if (layer_idx > layer_num[i+1]):    # forward hopping, to the surface
            hopEvent_dict[layer_idx][f_idx] += 1
        else:   # reverse hopping, away from the surface
            hopEvent_dict[layer_idx][r_idx] += 1
    """ assume final layer is reverse hopping to future layer """
    final_layer = layer_num[-1]
    final_idx = 'B'+str(final_layer)+'R' # reverse
    hopEvent_dict[final_layer][final_idx] += 1
            
    """ apparent rate constant, using combined total residence time """    
    for i in range(len(int_layer_count)):     #  loop over integer layer
        layer_idx = int_layer_count[i][0]   
#        print(layer_idx)
        f_idx = 'B'+str(layer_idx-1)+'F'   # forward
        r_idx = 'B'+str(layer_idx)+'R' # reverse 
        
                
        hop_forward_count = hopEvent_dict[layer_idx][f_idx]    # corresponding layer number
        hop_reverse_count = hopEvent_dict[layer_idx][r_idx]
        
        if layer_idx not in rateConst_dict:
            rateConst_dict[layer_idx] = {f_idx:0, r_idx:0}   # initialize each layer into forward/reverse hop
        #### convert time from [10 fs] to [ps], time*10fs/1000
        rateConst_dict[layer_idx][f_idx] = np.round(hop_forward_count/(int_layer_count[i][1]*TIME_UNIT), 2)
        rateConst_dict[layer_idx][r_idx] = np.round(hop_reverse_count/(int_layer_count[i][1]*TIME_UNIT), 2)
    

    hopEvent_dict = OrderedDict(sorted(hopEvent_dict.items()))
    rateConst_dict = OrderedDict(sorted(rateConst_dict.items()))
    
    
    with open('hopEventDict.csv',"w") as f:
        w = csv.writer(f)
        for k,v in hopEvent_dict.items():
            w.writerow([k,v])
            
    with open('rateConstDict.csv',"w") as f:
        w = csv.writer(f)
        for k,v in rateConst_dict.items():
            w.writerow([k,v])            
    
    return hopEvent_dict, rateConst_dict



def main():
    
    input_layer_array = np.loadtxt('./intLayer.dat',skiprows=1,dtype='int')
    int_layer_count = np.loadtxt('./intCount.dat',skiprows=1,dtype='int')


    RLE_dict = runLengthEncoding(input_layer_array)
    
    resTime_separate_dict,resTime_dict = calcResTimeDistribution(RLE_dict)
    
    hopEvent_dict, rateConst_dict = calcRateConstant(RLE_dict, int_layer_count)

    
    
if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- Runtime: %.2f seconds ---" % (time.time() - start_time))