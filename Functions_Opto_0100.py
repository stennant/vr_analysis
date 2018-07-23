# -*- coding: utf-8 -*-
"""
    
    @author: Sarah Tennant
    
    Includes all core functions which are imported into main code for analysis.
    
    """

# Import packages
import matplotlib.pyplot as plt
import numpy as np
import h5py
from matplotlib import rcParams
from Functions_Params_0100 import STOP_THRESHOLD, DIST, HDF_LENGTH, BINNR, SHUFFLE_N
from scipy.stats import uniform
from scipy import stats
import random


# -------------------------------------------------------------------------------------------------------------- #


def round_down(num, divisor):
	return num - (num%divisor)



# -------------------------------------------------------------------------------------------------------------- #



# CALCULATE FIRST STOPPING LOCATION PER TRIAL
def FirstStops_hist( trarray,stops ):
    data = []
    for row in trarray: # for each trial
        tarray = stops[stops[:,2] ==row,:] # get data only for each trial
        for row in tarray: # go through row of data for specified trial
            if float(row[0]) > 3.2 and float(row[0]) <= 11.2: # if stop is between black boxes
                data.append([float(row[0]), int(row[1]), int(row[2])]) # append data
                break # break so only get first stop then goes onto next trial
    return np.array(data)



# CALCULATE FIRST STOPPING LOCATION PER TRIAL, SPLIT OPTO AND NON OPTO
def FirstStops_Opto( trarray,stops ):
    data_l = []
    data_nl = []
    indicies = get_opto_trial_indicies(stops)
    #print(indicies.shape, 'indicies')
    for row in trarray: # for each trial
        tarray = stops[stops[:,2] ==row,:] # get data only for each trial
        opto = indicies[indicies[:,0] == row, 1]
        #print(opto,'opto')
        if opto == 1:
            for row in tarray: # go through row of data for specified trial
                if float(row[0]) > 3.2 and float(row[0]) <= 11.5: # if stop is between black boxes
                    data_l.append([float(row[0]), int(row[1]), int(row[2])]) # append data
                    break # break so only get first stop then goes onto next trial
        if opto == 0:
            for row in tarray: # go through row of data for specified trial
                if float(row[0]) > 3.2 and float(row[0]) <= 11.2: # if stop is between black boxes
                    data_nl.append([float(row[0]), int(row[1]), int(row[2])]) # append data
                    break # break so only get first stop then goes onto next trial

    return np.array(data_l),np.array(data_nl)



def get_opto_trial_indicies(stops):
    trials = np.unique(stops[:,2]) # array of unique trial numbers
    rounded_total_trials = round_down(np.amax(trials), 10)
    #print(np.amax(trials), rounded_total_trials)
    iterations = rounded_total_trials/10
    iterations_array = np.arange(0,iterations,1)
    trial_indicies = []
    for x in iterations_array:
        if x == 0:
            trial_indicies = np.append(trial_indicies, np.repeat(0,20))
        elif x % 2 != 0:
            trial_indicies = np.append(trial_indicies, np.repeat(1,10))
        else:
             trial_indicies = np.append(trial_indicies, np.repeat(0,10))
    trial_indicies = np.asarray(trial_indicies)
    trials = np.arange(1,len(trial_indicies)+1, 1)
    #print(trial_indicies.shape,trials.shape)
    data = np.vstack((trials,trial_indicies))
    data = np.transpose(data)
    #print(data.shape, 'iterations')
    #print(trial_indicies)           
    return data


# input: 0001111111100000111111 
def first_indicies(indicies):
        first_stim_indicies = np.zeros((len(indicies)))
        last_stim_indicies = np.zeros((len(indicies)))
        first_nonstim_indicies = np.zeros((len(indicies)))
        last_nonstim_indicies = np.zeros((len(indicies)))
        for rowcount,row in enumerate(indicies):
                difference = indicies[rowcount,1] + indicies[rowcount-1,1]
                #print(difference, indicies[rowcount,1])
                if difference == 1 and indicies[rowcount,1] ==0:
                        first_stim_indicies[rowcount+1]=1
                        last_nonstim_indicies[rowcount]=1
                elif difference== 1 and indicies[rowcount,1] ==1:
                         last_stim_indicies[rowcount]=1
                         first_nonstim_indicies[rowcount+1]=1

        trials = np.arange(1,len(indicies)+1, 1)
        data = np.vstack((trials,first_stim_indicies, last_stim_indicies, first_nonstim_indicies, last_nonstim_indicies))
        data = np.transpose(data)

        
        return data
              
                        
                




# CALCULATE FIRST STOPPING LOCATION PER TRIAL, SPLIT OPTO AND NON OPTO
def FirstStops_Opto_split( trarray,stops ):
    data_f_s = []
    data_f_ns = []
    data_l_s = []
    data_l_ns = []
    indicies = get_opto_trial_indicies(stops)
    data = first_indicies(indicies)
    #print(indicies.shape, 'indicies')
    for row in trarray: # for each trial
        tarray = stops[stops[:,2] ==row,:] # get data only for each trial
        #load indicies to mark stimulation events
        first_stim = data[data[:,0] == row, 1]
        last_stim = data[data[:,0] == row, 2]
        first_nonstim = data[data[:,0] == row, 3]
        last_nonstim = data[data[:,0] == row, 4]
       #print(opto,'opto')
        if first_stim == 1:
            for row in tarray: # go through row of data for specified trial
                if float(row[0]) > 3.2 and float(row[0]) <= 17.5: # if stop is between black boxes
                    data_f_s.append([float(row[0]), int(row[1]), int(row[2])]) # append data
                    break # break so only get first stop then goes onto next trial
        elif last_stim == 1:
            for row in tarray: # go through row of data for specified trial
                if float(row[0]) > 3.2 and float(row[0]) <= 17.2: # if stop is between black boxes
                    data_l_s.append([float(row[0]), int(row[1]), int(row[2])]) # append data
                    break # break so only get first stop then goes onto next trial

        elif first_nonstim == 1:
            for row in tarray: # go through row of data for specified trial
                if float(row[0]) > 3.2 and float(row[0]) <= 17.2: # if stop is between black boxes
                    data_f_ns.append([float(row[0]), int(row[1]), int(row[2])]) # append data
                    break # break so only get first stop then goes onto next trial

        elif last_nonstim == 1:
            for row in tarray: # go through row of data for specified trial
                if float(row[0]) > 3.2 and float(row[0]) <= 17.2: # if stop is between black boxes
                    data_l_ns.append([float(row[0]), int(row[1]), int(row[2])]) # append data
                    break # break so only get first stop then goes onto next trial

    return np.array(data_f_s),np.array(data_f_ns), np.array(data_l_s),np.array(data_l_ns)



#-------------------------------------------------------------------------------------------------------------- #





# Input: array[:,4] (columns: location, time, trialno, reward, empty)
# Output: array[trialnumbers, locationbins]
# Function: Creates histogram of stops in bins
# BIN STOPS INTO 20, 10 CM LOCATION BINS
def create_srdata( stops, trialids ):
    if stops.size == 0:
        return np.zeros((BINNR,))
    
    # create histogram
    posrange = np.linspace(0, HDF_LENGTH, num=BINNR+1) # 0 VU to 20 VU split into 20
    trialrange = trialids
    trialrange = np.append(trialrange, trialrange[-1]+1)  # Add end of range
    values = np.array([[trialrange[0], trialrange[-1]],[posrange[0], posrange[-1]]])

    H, bins, ranges = np.histogram2d(stops[:,2], stops[:,0], bins=(trialrange, posrange), range=values)
    H[np.where(H[::]>1)] = 1

    return H



# SHUFFLE STOPS
def shuffle_stops2( stops,n ):
    shuffled_stops = np.copy(stops) # this is required as otherwise the original dataset would be altered
    # create an array that contains the amount by which every stop will be shuffled
    rand_rotation = uniform.rvs(loc=0, scale=HDF_LENGTH, size=stops.shape[0])
    # add random value
    shuffled_stops[:,0] = rand_rotation
    shuffled_stops[:,2] = n
    
    return shuffled_stops



# Input: array[:,4] (columns: location, time, trialno, reward, zeros), array[unique trialnumbers]
# Output: array[20], array[20], array[20], array[20]
# Function: creates shuffled stops datasets from real dataset
# CREATE SHUFFLED DATASETS
def shuffle_analysis_pertrial3(stopsdata, trialids):
    if stopsdata.size == 0:
        return np.zeros((BINNR, )), np.zeros((BINNR, )), np.zeros((BINNR, )), np.zeros((BINNR, ))
    SHUFFLE1 = 100
    # Calculate stop rate for each section of the track
    srbin = create_srdata( stopsdata, trialids )                        # Array(BINNR, trialnum)
    srbin_mean = np.mean(srbin, axis=0)                                 # Array(BINNR)
    srbin_std = stats.sem(srbin, axis=0)                                 # Array(BINNR)
    # Shuffling random 100 trials 1000 times
    shuffled_srbin_mean = np.zeros((SHUFFLE_N, BINNR))
    for i in range(SHUFFLE_N): # for i in 1000
        shuffledtrials = np.zeros((SHUFFLE1, 5))
        shuffleddata =np.zeros((SHUFFLE1, BINNR))
        for n in range(SHUFFLE1): # Create sample data with 100 trials
            trial = random.choice(trialids) # select random trial from real dataset
            data = stopsdata[stopsdata[:,2] ==trial,:] # get data only for each tria
            shuffledtrial = shuffle_stops2(data,n) # shuffle the locations of stops in the trial
            shuffledtrials = np.vstack((shuffledtrials,shuffledtrial)) # stack shuffled stops
        trialids2 = np.unique(shuffledtrials[:, 2]) # find unique trials in the data
        shuffled_srbin = create_srdata( shuffledtrials, trialids2 ) #
        shuffled_srbin_mean[i] = np.mean(shuffled_srbin, axis=0)        # Array(BINNR)
    # Mean of the mean stops in the shuffled data for each bin    
        
    shuffled_mean = np.mean(shuffled_srbin_mean, axis=0)                # Array(BINNR)
    shuffled_std = np.std(shuffled_srbin_mean, axis=0)                  # Array(BINNR)
    return srbin_mean, srbin_std, shuffled_mean, shuffled_std




# Input: array[20,trialids], array[20,trialids]

def split_light_nolight(stops, trarray):
    data_l = np.zeros((0,5))
    data_nl = np.zeros((0,5))
    indicies = get_opto_trial_indicies(stops)

    l_indicies = np.delete(indicies, np.where(indicies[:,1] == 1),0) 
    nl_indicies = np.delete(indicies, np.where(indicies[:,1] == 0),0)
    l_indicies = np.asarray(l_indicies[:,0])
    nl_indicies = np.asarray(nl_indicies[:,0])
    #print(l_indicies.shape,nl_indicies,'indicies')

    for row in trarray: # for each trial
        tarray = stops[stops[:,2] ==row,:] # get data only for each trial
        if row in l_indicies:
            data_l = np.vstack((data_l,tarray))
        else:
            data_nl = np.vstack((data_nl,tarray))
 
    #print(data_l.shape, data_nl.shape,'split data')
    return data_l,data_nl
    




# Input: array[20,trialids], array[20,trialids]

def split_light_nolight_speed(stops, trarray):
    data_l = np.zeros((0,13))
    data_nl = np.zeros((0,13))
    indicies = get_opto_trial_indicies_speed(stops)

    l_indicies = np.delete(indicies, np.where(indicies[:,1] == 1),0) 
    nl_indicies = np.delete(indicies, np.where(indicies[:,1] == 0),0)
    l_indicies = np.asarray(l_indicies[:,0])
    nl_indicies = np.asarray(nl_indicies[:,0])
    #print(l_indicies.shape,nl_indicies,'indicies')

    for row in trarray: # for each trial
        tarray = stops[stops[:,9] ==row,:] # get data only for each trial
        if row in l_indicies:
            data_l = np.vstack((data_l,tarray))
        else:
            data_nl = np.vstack((data_nl,tarray))
 
    #print(data_l.shape, data_nl.shape,'split data')
    return data_l,data_nl
    




def get_opto_trial_indicies_speed(stops):
    trials = np.unique(stops[:,9]) # array of unique trial numbers
    rounded_total_trials = round_down(np.amax(trials), 10)
    #print(np.amax(trials), rounded_total_trials)
    iterations = rounded_total_trials/10
    iterations_array = np.arange(0,iterations,1)
    trial_indicies = []
    for x in iterations_array:
        if x == 0:
            trial_indicies = np.append(trial_indicies, np.repeat(0,20))
        elif x % 2 != 0:
            trial_indicies = np.append(trial_indicies, np.repeat(1,10))
        else:
             trial_indicies = np.append(trial_indicies, np.repeat(0,10))
    trial_indicies = np.asarray(trial_indicies)
    trials = np.arange(1,len(trial_indicies)+1, 1)
    #print(trial_indicies.shape,trials.shape)
    data = np.vstack((trials,trial_indicies))
    data = np.transpose(data)
    #print(data.shape, 'iterations')
    #print(trial_indicies)           
    return data



def find_speed_diff(data):
        differences = []
        for tcount, trial in enumerate(data[1]):
                diff = data[3,tcount] - data[9,tcount]
                differences = np.append(differences,diff)
        #difference=np.nanmean(differences)
        return differences
