import glob
import os, sys
import numpy as np
from icecube import dataio, simclasses, millipede, recclasses
from math import *
import simweights


def get_data(filename):
    f = open(filename, "r")
    
    temp = f.readline().split()
    num_files = float(temp[0])

    temp = f.readline().split()
    phi_list = np.linspace(float(temp[0]), float(temp[1]), int(temp[2]))

    temp = f.readline().split()
    cos_theta_list = np.linspace(float(temp[0]), float(temp[1]), int(temp[2]))

    box = np.zeros((len(cos_theta_list) - 1, len(phi_list) - 1))

    curr_theta = 0

    while 1:
        temp = f.readline().split()

        if len(temp) < 1:
            break

        for num in range(len(temp)):
            box[curr_theta][num] = float(temp[num])

        curr_theta += 1

    return box, phi_list, cos_theta_list, num_files


def get_place(list, target):

    if target < list[0]:
        return "nan"
    if target > list[len(list) - 1]:
        return "nan"

    if str(target) == 'nan':
        return 'nan'

    for i in range(len(list) - 1):
        if target >= list[i] and target <= list[i + 1]:
            return i

    print("should never happen")
    exit()


def classify(cos_theta_list, phi_list, theta, phi):

    theta_list = np.arccos(cos_theta_list)

    phi_index = get_place(phi_list, phi)

    theta_index = get_place(theta_list, theta)

    return theta_index, phi_index


def check_muon_filter(frame):
    try:
        return frame["FilterMask"]["MuonFilter_13"].condition_passed
    except:
        print(111)
        return False

        
def write_to_file(filename):
    
    weight_keys = [
        "CylinderLength",
        "CylinderRadius",
        "EnergyPrimaryMax",
        "EnergyPrimaryMin",
        "NEvents",
        "OverSampling",
        "ParticleType",
        "PrimaryEnergy",
        "PrimarySpectralIndex",
        "PrimaryType",
        "ThetaMax",
        "ThetaMin",
        "Weight",
    ]

    particle_keys = ["type", "energy", "zenith"]

    Muon_llh = {k: [] for k in ["zenith", "azimuth"]}


    inFile_corsika = dataio.I3File(filename)
    
    start_time = ''
    end_time = ''
    p_count = 0
    
    while inFile_corsika.more():

        try: 
            frame = inFile_corsika.pop_physics()
            #if "FilterMask" in frame:
            #if frame["FilterMask"]["MuonFilter_13"].condition_passed:
        except:
            print("empty file, jump next")
            f = open("empty_file.txt", 'a')
            f.write(filename + "\n")
            f.close()
            return
        
        if p_count == 0:
            start_time = frame["I3EventHeader"].start_time
            print("start at " + str(start_time))
            p_count += 1
            
                                
            
        if "MPEFit" in frame and check_muon_filter(frame):
            
                
            Muon_llh["zenith"].append(frame["MPEFit"].dir.zenith)
            Muon_llh["azimuth"].append(frame["MPEFit"].dir.azimuth)
            
    end_time = frame["I3EventHeader"].end_time
    print("end at " + str(end_time))
    
    duration = float((end_time - start_time))/1e9
    
    print("duration " + str(duration) + " seconds")
    
    
    print("num of LLH data after filter: " + str(len(Muon_llh['zenith'])))
    

    phi_size = 101
    theta_size = 101
    cos_theta_list = np.linspace(1, -1, theta_size)
    phi_list = np.linspace(0, 2 * pi, phi_size)
    
    """
    doing classification for llh
    """
    box_llh = np.zeros((theta_size - 1, phi_size - 1))
    
    for i in range(len(Muon_llh['azimuth'])):

        LLHAzi = Muon_llh['azimuth'][i]
        LLHZen = Muon_llh['zenith'][i]

        theta_index, phi_index = classify(cos_theta_list, phi_list, LLHZen, LLHAzi)


        if not (theta_index == "nan" or phi_index == "nan") and not (LLHAzi == 0 and LLHZen == 0):
            box_llh[theta_index][phi_index] += 1
            
            
            
            
    if not os.path.exists("./output_burn_mpe_muon.txt"):
  
        f = open("output_burn_mpe_muon.txt", 'w')
        f.write(str(duration) + "\n")
        f.write(str(phi_list[0]) + " " + str(phi_list[-1]) + " " + str(len(phi_list)) + "\n")
        f.write(str(cos_theta_list[0]) + " " + str(cos_theta_list[-1]) + " " + str(len(cos_theta_list)) + "\n")

        for i in range(len(box_llh)):
            for j in range(len(box_llh[i])):
                f.write(str(box_llh[i][j]) + " ")

            f.write("\n")

        f.close()
        
        print("\n")

        
        
    else:
        this_events = float(box_llh.sum())
        box_llh_temp, phi_list, cos_theta_list, time_passed = get_data("./output_burn_mpe_muon.txt")
        box_llh += box_llh_temp
        f = open("./output_burn_mpe_muon.txt", 'w')
        f.write(str(time_passed + duration) + "\n")
        f.write(str(phi_list[0]) + " " + str(phi_list[-1]) + " " + str(len(phi_list)) + "\n")
        f.write(str(cos_theta_list[0]) + " " + str(cos_theta_list[-1]) + " " + str(len(cos_theta_list)) + "\n")

        for i in range(len(box_llh)):
            for j in range(len(box_llh[i])):
                f.write(str(box_llh[i][j]) + " ")

            f.write("\n")

        f.close()
        
        
        total_events = float(box_llh.sum())
        print("# events here: " + str(this_events) + " rate: " + str(this_events/duration))
        print("total rate: " + str(total_events/(time_passed + duration)))
        print("\n")
        
        
        

    
    

#filename = '/data/exp/IceCube/2018/filtered/level2/0110/Run00130533/Level2_IC86.2017_data_Run00130533_Subrun00000000_00000009.i3.zst'

path = "/data/exp/IceCube/2018/filtered/level2/"

if not os.path.exists("./I3_data/"):
    os.mkdir("./I3_data/")

for index in os.listdir(path):
    
    if not os.path.exists("./I3_data/" + str(index)):
        os.mkdir("./I3_data/" + str(index))
        
    for files in os.listdir(path + str(index)):
        if str(files)[0] == "R" and str(files)[-3] != "_" and str(files)[-1] == "0":

            if os.path.exists("./I3_data/" + str(index) + "/" + str(files)):
                print("find " + str(index) + "/" + str(files))
                continue

            if not os.path.exists("./I3_data/" + str(index) + "/" + str(files)):
                os.mkdir("./I3_data/" + str(index) + "/" + str(files))

            path_2 = path+ str(index) + "/" + str(files)

            for inner_files in os.listdir(path_2):
                if inner_files.endswith(".zst") and inner_files[-8] != "T":
                    print("Reading", path_2 + "/" + str(inner_files))
                    write_to_file(path_2 + "/" + str(inner_files))
    
                
