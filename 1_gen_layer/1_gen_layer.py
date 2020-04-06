#!/usr/bin/env python
# Xiaohong Zhang
# Getman Research Group
# Mar 27, 2018

import sys, os
import re
import numpy as np
from collections import OrderedDict
import argparse
import time


# Usage: put 1_gen_layer.py and O_COM.out in the same folder or type full path of O_COM.out file folder
# python 1_gen_layer.py O_COM.out #C_COM.out COM_COM.out ...
# output files: decLayer.dat, decCount.dat, intLayer.dat, intCount.dat

##################### manu input starts ##################################
DUMP_FREQ = 10                                        # dumping frequency of coords 
TIME_STEP = 1.0                                       # [fs]
TIME_UNIT = DUMP_FREQ*TIME_STEP/1000.0                # convert [fs] time unit to [ps]
MOL_DIAMETER = 3.7                                    # this is MeOH diameter, used as single layer width
DEC_LAYER = 0.1                                       # decimal layer
INT_LAYER = 1.0                                       # integer layer
DEC_LAYER_WIDTH = MOL_DIAMETER * DEC_LAYER 
INT_LAYER_WIDTH = MOL_DIAMETER * INT_LAYER 
DEC_START_Z = 5.958035                                # there are Pt surfaces below
INT_START_Z = 6.698035                                # after count on Pt radius and the first 0.2 layer of MeOH, start from 0.3 layer
##################### manu input ends ##################################


def make_folder(dir_name):
    try:
        os.makedirs(dir_name)
    except OSError:
        if not os.path.isdir(dir_name):
            raise
            
                
def readFile(file_name, output_folder):
    """Read from LAMMPS center-of-mass output"""
    data = np.genfromtxt(file_name, skip_header=3, dtype='unicode')
#    print(data)
    data_length = len(data)
    data_array = np.empty((int(data_length/4), 4), dtype='object')               

        
    for i in range(data_length):                       # Convert data structure every 4 lines into 1 line
        row = int(i / 4)
        col = int(i % 4)

        if (col == 0):
            data_array[row][col] = data[i][0]
#            print(type(data_array[row][col]))
        else:
            data_array[row][col] = data[i][1]

    np.savetxt(output_folder+'/dataArr.dat', data_array, header="TimeStep\tx\ty\tz", fmt="%s\t%15s\t%15s\t%15s")
#    print(data_array)
    return data_array




def calcLayerNumber(data_array, output_folder):
    """ cal both decimal and integer layer number based on z height """
        
    dec_layer_arr = np.empty((len(data_array),2), dtype='object')      # save time step + decimal layer number 
    int_layer_arr = np.empty((len(data_array),2), dtype='object')      # save time step + integer layer number 
    
    dec_layer_count_dict = OrderedDict()               # count total number in each decimal layer
    int_layer_count_dict = OrderedDict()               # count total number in each integer layer
        
    
    for i in range(len(data_array)):                   # loop over time steps
        z = float(data_array[i][3])                    # extract z value

                
        """ decimal layer calculation """
        dec_layer =  str('{:.1f}'.format(np.ceil((z-DEC_START_Z)/DEC_LAYER_WIDTH) * DEC_LAYER))    # get float(layer_number)
#        print(type(dec_layer))
#        print(dec_layer)
        dec_layer_arr[i][0] = data_array[i][0]     # get time step
        dec_layer_arr[i][1] = dec_layer # corresponding layer number
        
        if dec_layer in dec_layer_count_dict:
            dec_layer_count_dict[dec_layer]+=1          # update each value (layer count) to a key (layer number, int type)
        else:
            dec_layer_count_dict[dec_layer]=1           # assign initial value 1 (layer count) to a key (layer number, int type)

            
        """ integer layer calculation """
        int_layer = str('{:.0f}'.format(np.ceil((z-INT_START_Z)/INT_LAYER_WIDTH))) # get integer layer number
#        print(type(int_layer))
#        print(int_layer)
        if (int_layer == '0' or int_layer == '-0'):                         # count 0.2 layer into interger layer 1
            int_layer = '1'        
        
        int_layer_arr[i][0] = data_array[i][0]         # get time step
        int_layer_arr[i][1] = int_layer                # corresponding layer number
        
        if int_layer in int_layer_count_dict:
            int_layer_count_dict[int_layer]+=1         # update each value (layer count) to a key (layer number, int type)
        else:
            int_layer_count_dict[int_layer]=1          # assign initial value 1 (layer count) to a key (layer number, int type)

    dec_layer_count_dict = OrderedDict(sorted(dec_layer_count_dict.items()))    
    int_layer_count_dict = OrderedDict(sorted(int_layer_count_dict.items())) 
            
            
#    print(dec_layer_arr)        
#    print(int_layer_arr)          
    """ save decimal layer file """
    np.savetxt(output_folder + '/decLayer.dat', dec_layer_arr, header="TimeStep\tDecLayer", fmt='%s\t%s') 
    with open(output_folder + '/decCount.dat',"w") as f:
        f.writelines("#DecLayer\Count\n")
        f.writelines('{} {}\n'.format(k,v) for k,v in dec_layer_count_dict.items())    
                         
                         
    """ save integer layer file """
    np.savetxt(output_folder + '/intLayer.dat', int_layer_arr, header="TimeStep\tIntLayer", fmt='%s\t%s') 
    with open(output_folder + '/intCount.dat',"w") as f:
        f.writelines("IntLayer\Count\n")
        f.writelines('{} {}\n'.format(k,v) for k,v in int_layer_count_dict.items())                           

    return dec_layer_arr,int_layer_arr

def main():   
    
    for i in range(len(sys.argv) - 1):   # loop over all files
    
        file_name = sys.argv[i+1]
        
        output_folder = re.split(r'_',file_name)[0]  + '_results'
        make_folder(output_folder)
        
        data_array = readFile(file_name,output_folder)    # convert COM 4 line format into 1 line
                    
        dec_layer_arr,int_layer_arr = calcLayerNumber(data_array, output_folder)         # write both decimal and integer layer to separate files
    
    
if __name__ == "__main__":
    start_time = time.time()
    main()            
    print("--- Runtime: %.2f seconds ---" % (time.time() - start_time))    
    
