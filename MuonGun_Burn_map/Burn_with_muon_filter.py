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

    Muon_mpe = {k: [] for k in ["zenith", "azimuth", "MuEX"]}


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
            
                                
            
        if "MPEFitMuEX" in frame and check_muon_filter(frame):
            
            Muon_mpe["zenith"].append(frame["MPEFitMuEX"].dir.zenith)
            Muon_mpe["azimuth"].append(frame["MPEFitMuEX"].dir.azimuth)
            Muon_mpe["MuEX"].append(frame["MPEFitMuEX"].energy)
            
    end_time = frame["I3EventHeader"].end_time
    print("end at " + str(end_time))
    
    duration = float((end_time - start_time))/1e9
    
    print("duration " + str(duration) + " seconds")
    
    
    print("num of LLH data after filter: " + str(len(Muon_mpe['zenith'])))
    

    phi_size = 51
    theta_size = 51
    cos_theta_list = np.linspace(1, -1, theta_size)
    phi_list = np.linspace(0, 2 * pi, phi_size)
    
    """
    doing classification for llh
    """
    box_mpe_tot = np.zeros((theta_size - 1, phi_size - 1))
    box_mpe_200 = np.zeros((theta_size - 1, phi_size - 1))
    box_mpe_500 = np.zeros((theta_size - 1, phi_size - 1))
    box_mpe_1000 = np.zeros((theta_size - 1, phi_size - 1))
    box_mpe_1000p = np.zeros((theta_size - 1, phi_size - 1))
    
    for i in range(len(Muon_mpe['azimuth'])):

        MPEAzi = Muon_mpe['azimuth'][i]
        MPEZen = Muon_mpe['zenith'][i]
        MPEMuEX = Muon_mpe['MuEX'][i]

        theta_index, phi_index = classify(cos_theta_list, phi_list, MPEZen, MPEAzi)


        if not (theta_index == "nan" or phi_index == "nan") and not (MPEAzi == 0 and MPEZen == 0):
            box_mpe_tot[theta_index][phi_index] += 1
            
            if MPEMuEX <= 200:
                box_mpe_200[theta_index][phi_index] += 1
            elif MPEMuEX > 200 and MPEMuEX <= 500:
                box_mpe_500[theta_index][phi_index] += 1
            elif MPEMuEX > 500 and MPEMuEX <= 1000:
                box_mpe_1000[theta_index][phi_index] += 1
            else:
                box_mpe_1000p[theta_index][phi_index] += 1
            
            
    box_all = [box_mpe_tot, box_mpe_200, box_mpe_500, box_mpe_1000, box_mpe_1000p]
        
    for energy_index in range(5):
        if not os.path.exists("./output_burn_mpe_muon_" + str(energy_index) + ".txt"):

            f = open("./output_burn_mpe_muon_" + str(energy_index) + ".txt", 'w')
            f.write(str(duration) + "\n")
            f.write(str(phi_list[0]) + " " + str(phi_list[-1]) + " " + str(len(phi_list)) + "\n")
            f.write(str(cos_theta_list[0]) + " " + str(cos_theta_list[-1]) + " " + str(len(cos_theta_list)) + "\n")

            for i in range(len(box_all[energy_index])):
                for j in range(len(box_all[energy_index][i])):
                    f.write(str(box_all[energy_index][i][j]) + " ")

                f.write("\n")

            f.close()

            print("\n")


        else:
            this_events = float(box_all[energy_index].sum())
            box_mpe_temp, phi_list, cos_theta_list, time_passed = get_data("./output_burn_mpe_muon_" + str(energy_index) + ".txt")
            box_all[energy_index] += box_mpe_temp
            f = open("./output_burn_mpe_muon_" + str(energy_index) + ".txt", 'w')
            f.write(str(time_passed + duration) + "\n")
            f.write(str(phi_list[0]) + " " + str(phi_list[-1]) + " " + str(len(phi_list)) + "\n")
            f.write(str(cos_theta_list[0]) + " " + str(cos_theta_list[-1]) + " " + str(len(cos_theta_list)) + "\n")

            for i in range(len(box_all[energy_index])):
                for j in range(len(box_all[energy_index][i])):
                    f.write(str(box_all[energy_index][i][j]) + " ")

                f.write("\n")

            f.close()


            total_events = float(box_all[energy_index].sum())
            print("# events in range " + str(energy_index) + ": " + str(this_events) + " rate: " + str(this_events/duration))
            print("total rate in range " + str(energy_index) + ": " + str(total_events/(time_passed + duration)))
            print("\n")
        
        
        

    
    

#filename = '/data/exp/IceCube/2018/filtered/level2/0110/Run00130533/Level2_IC86.2017_data_Run00130533_Subrun00000000_00000009.i3.zst'

path = "/data/exp/IceCube/2018/filtered/level2/"

if not os.path.exists("./I3_data/"):
    os.mkdir("./I3_data/")

for index in os.listdir(path):
    
    if not index.isdigit():
        print(index)
        continue
    
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
    
                
