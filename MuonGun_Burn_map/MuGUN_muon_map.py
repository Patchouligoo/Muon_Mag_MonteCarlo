import glob
import os, sys
import numpy as np
from icecube import dataio, MuonGun, icetray, dataclasses
from icecube.MuonGun import BundleEntry, BundleConfiguration
from math import *
from icecube.icetray.i3logging import log_info as log




def harvest_generators(infiles):
    """
    Harvest serialized generator configurations from a set of I3 files.
    """
    generator = None
    for fname in infiles:
        f = dataio.I3File(fname)
        fr = f.pop_frame(icetray.I3Frame.Stream('S'))
        f.close()
        if fr is not None:
            for k in fr.keys():
                v = fr[k]
                if isinstance(v, MuonGun.GenerationProbability):
                    log('%s: found "%s" (%s)' % (fname, k, type(v).__name__), unit="MuonGun")
                    if generator is None:
                        generator = v
                    else:
                        generator += v
    return generator



def get_data(filename):
    f = open(filename, "r")
    
    temp = f.readline().split()
    num_files = int(temp[0])

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

    Muon_mpe = {k: [] for k in ["zenith", "azimuth", "weight", "MuEX"]}


    inFile_MuGun = dataio.I3File(filename)
    
    model = MuonGun.load_model('GaisserH4a_atmod12_SIBYLL')
    
    generator = harvest_generators([filename])

    weighter = MuonGun.WeightCalculator(model, generator)
    
    
    while inFile_MuGun.more():

        try: 
            frame = inFile_MuGun.pop_physics()
        except:
            print("empty file, jump next")
            f = open("empty_file.txt", 'a')
            f.write(filename + "\n")
            f.close()
            return
                
            
        if "MPEFitMuEX" in frame and check_muon_filter(frame):
                
            Muon_mpe["zenith"].append(frame["MPEFit"].dir.zenith)
            Muon_mpe["azimuth"].append(frame["MPEFit"].dir.azimuth)
            Muon_mpe['MuEX'].append(frame["MPEFitMuEX"].energy)
            
#             x = frame['I3MCTree_preMuonProp'][1].dir.x
#             y = frame['I3MCTree_preMuonProp'][1].dir.y
#             z = frame['I3MCTree_preMuonProp'][1].dir.z
#             zenith = frame['I3MCTree_preMuonProp'][1].dir.zenith
#             azimuth = frame['I3MCTree_preMuonProp'][1].dir.azimuth
#             energy = frame['I3MCTree_preMuonProp'][1].energy
            
#             axis = dataclasses.I3Particle()
#             axis.pos = dataclasses.I3Position(x,y,z)
#             axis.dir = dataclasses.I3Direction(zenith,azimuth)
#             bundle = MuonGun.BundleConfiguration()
#             bundle.append(MuonGun.BundleEntry(float(0),float(energy)))
            
#             weights = weighter(axis, bundle)

            mctree = frame['I3MCTree_preMuonProp']   
            primary = mctree.primaries[0]
            muon = mctree.get_daughters(primary)[0]
            bundle = MuonGun.BundleConfiguration(
                [MuonGun.BundleEntry(0, muon.energy)])
            weighter = MuonGun.WeightCalculator(model, generator)
            weight   = weighter(primary,bundle)
            
            Muon_mpe["weight"].append(weight)

    print("num of MPE data after filter: " + str(len(Muon_mpe['zenith'])) )
    

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
            
            box_mpe_tot[theta_index][phi_index] += Muon_mpe["weight"][i]
            
            if MPEMuEX <= 200:
                box_mpe_200[theta_index][phi_index] += Muon_mpe["weight"][i]
            elif MPEMuEX > 200 and MPEMuEX <= 500:
                box_mpe_500[theta_index][phi_index] += Muon_mpe["weight"][i]
            elif MPEMuEX > 500 and MPEMuEX <= 1000:
                box_mpe_1000[theta_index][phi_index] += Muon_mpe["weight"][i]
            else:
                box_mpe_1000p[theta_index][phi_index] += Muon_mpe["weight"][i]
            
            
    box_all = [box_mpe_tot, box_mpe_200, box_mpe_500, box_mpe_1000, box_mpe_1000p]
            
    total_rate = 0
    
    for energy_index in range(5):       
        if not os.path.exists("./output_mugun_mpe_muon_" + str(energy_index) + ".txt"):

            f = open("./output_mugun_mpe_muon_" + str(energy_index) + ".txt", 'w')
            f.write(str(1) + "\n")
            f.write(str(phi_list[0]) + " " + str(phi_list[-1]) + " " + str(len(phi_list)) + "\n")
            f.write(str(cos_theta_list[0]) + " " + str(cos_theta_list[-1]) + " " + str(len(cos_theta_list)) + "\n")

            for i in range(len(box_all[energy_index])):
                for j in range(len(box_all[energy_index][i])):
                    f.write(str(box_all[energy_index][i][j]) + " ")

                f.write("\n")

            f.close()
            
            print("curr rate: " + str(box_all[energy_index].sum()) + "Hz")
            if energy_index == 0:
                total_rate = box_all[energy_index].sum()

        else:     

            print("curr rate: " + str(box_all[energy_index].sum()) + "Hz")

            box_mpe_temp, phi_list, cos_theta_list, num_files = get_data("./output_mugun_mpe_muon_" + str(energy_index) + ".txt")
            box_all[energy_index] += box_mpe_temp
            f = open("./output_mugun_mpe_muon_" + str(energy_index) + ".txt", 'w')
            f.write(str(num_files + 1) + "\n")
            f.write(str(phi_list[0]) + " " + str(phi_list[-1]) + " " + str(len(phi_list)) + "\n")
            f.write(str(cos_theta_list[0]) + " " + str(cos_theta_list[-1]) + " " + str(len(cos_theta_list)) + "\n")

            for i in range(len(box_all[energy_index])):
                for j in range(len(box_all[energy_index][i])):
                    f.write(str(box_all[energy_index][i][j]) + " ")

                f.write("\n")

            f.close()
            
            if energy_index == 0:
                total_rate = box_all[energy_index].sum() / (num_files + 1)
            
            
            
    print("total rate: " + str(total_rate) + "Hz")
    print("\n\n")
        
        
        

    
    

#filename = '/data/sim/IceCube/2016/filtered/level2/CORSIKA-in-ice/20904/0095000-0095999/dst/IC86.2016_corsika.020904.095802.dst.i3.zst'

path = "/data/sim/IceCube/2016/filtered/level2/MuonGun/21950/"

if not os.path.exists("./I3_data/"):
    os.mkdir("./I3_data/")


for files in os.listdir(path):
    if str(files)[0] == str(0):

        if os.path.exists("./I3_data/" + str(files)):
            print("find " + str(files))
            continue

        if not os.path.exists("./I3_data/" + str(files)):
            os.mkdir("./I3_data/" + str(files))

        path_2 = path + str(files)

        for inner_files in os.listdir(path_2):
            if inner_files.endswith(".zst"):
                print("Reading", path_2 + "/" + str(inner_files))
                write_to_file(path_2 + "/" + str(inner_files))
    