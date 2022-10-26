#!/usr/bin/env python3
import numpy as np
from glob import glob as G
import os

def Get_Molecular_Orbitals(moldenfile):
    with open(moldenfile) as f:
        lines = f.readlines()
    ene_lines = []
    for i,line in enumerate(lines):
        if "Ene=" in line:
            ene_lines.append(i)
    mol_orbs = []
    n_orbs   = len(ene_lines)
    for i in range(n_orbs-1):
        orb = []
        for line in lines[ene_lines[i]+3:ene_lines[i+1]]:
            orb.append(float(line.split()[-1]))
        mol_orbs.append(orb)
    orb = []
    for line in lines[ene_lines[n_orbs-1]+3:]:
        orb.append(float(line.split()[-1]))
    mol_orbs.append(orb)
    return mol_orbs

def Get_Sign_Array(OLD_MOL_FILE,NEW_MOL_FILE):
    old_mol_orbs = Get_Molecular_Orbitals(OLD_MOL_FILE)
    new_mol_orbs = Get_Molecular_Orbitals(NEW_MOL_FILE)
    sign_array = []
    for i in range(len(old_mol_orbs)):
        sign_array.append(np.sign(np.dot(old_mol_orbs[i],new_mol_orbs[i])))
    return sign_array

def Write_Phase_Swapped_Molden(moldenfile,sign_array):
    with open(moldenfile) as f:
        newlines = f.readlines()
    with open(moldenfile,"w") as f:
        for i,line in enumerate(newlines):
            f.write(line)
            if "[MO]" in line:
                break
        start_mo_lines = i+1
        ene_lines_encountered = -1
        for line in newlines[start_mo_lines:]:
            if "Ene=" in line:
                ene_lines_encountered += 1
                f.write(line)
            elif "Spin=" in line:
                f.write(line)
            elif "Occup=" in line:
                f.write(line)
            else:
                templine = line.split()
                if sign_array[ene_lines_encountered] == -1:
                    if "-" in templine[1]:
                        templine[1] = templine[1].replace("-"," ")
                    else:
                        templine[1] = "-"+templine[1].strip()
                newline = f"{templine[0]:>5}{float(templine[1]):>11.05f}\n"
                f.write(newline)

def VMDInitialize(inputfilename):
    with open(inputfilename,"w") as f:
        f.write("display update on\n")
        f.write("color add item Display Background white\n")
        f.write("color Display Background white\n")
        f.write("display projection perspective \n")
        f.write("display culling off\n")
        f.write("axes location off\n")
        f.write("display rendermode Normal\n")
        f.write("display depthcue off\n")
        f.write("display resize 1920 1080\n")
    return None

def Write_VMD_Command_File(ORDERED_MOLDEN_FILELIST,orbital_number,color1,color2):
    vmd_color_ids={"blue":0,"red":1,"gray":2,"orange":3,"yellow":4,
                   "tan":5,"silver":6,"green":7,"white":8,"pink":9,
                   "cyan":10,"purple":11,"lime":12,"mauve":13,"ochre":14,
                   "iceblue":15,"black":16}
    if color1 not in [key for key in vmd_color_ids.keys()]:
        print(f"Color {color1} not found in VMD dictionary.  Defaulting to blue for color 1.")
        color1 = "blue"
    if color2 not in [key for key in vmd_color_ids.keys()]:
        print(f"Color {color2} not found in VMD dictionary.  Defaulting to red for color 2.")
        color2 = "red"
    with open("orbital_trajectory.vmd","a") as f:
        num_moldens_loaded = 0
        for moldenfile in ORDERED_MOLDEN_FILELIST:
            if num_moldens_loaded == 0:
                f.write(f"mol new {moldenfile} type molden first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n")
            else:
                f.write(f"mol addfile {moldenfile} type molden first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n")
            num_moldens_loaded += 1
        f.write("mol delrep 0 top\n")
        f.write("mol representation Licorice 0.100000 100.000000 100.000000\n")
        f.write("mol color Element\n")
        f.write("mol selection {all}\n")
        f.write("mol material Glossy\n")
        f.write("mol addrep top\n")
        f.write("mol selupdate 0 top 0\n")
        f.write("mol colupdate 0 top 0\n")
        f.write("mol scaleminmax top 0 0.000000 0.000000\n")
        f.write("mol smoothrep top 0 0\n")
        f.write("mol drawframes top 0 {now}\n")
        f.write(f"mol representation Orbital 0.050000 {orbital_number} 0 0 0.050 1 6 0 0 1\n")
        f.write(f"mol color ColorID {vmd_color_ids[color1]}\n")
        f.write("mol selection {all}\n")
        f.write(f"mol material Glass3\n")
        f.write(f"mol addrep top\n")
        f.write(f"mol selupdate 1 top 0\n")
        f.write(f"mol colupdate 1 top 0\n")
        f.write(f"mol scaleminmax top 1 0.000000 0.000000\n")
        f.write(f"mol smoothrep top 1 0\n")
        f.write("mol drawframes top 1 {now}\n")
        f.write(f"mol representation Orbital -0.050000 {orbital_number} 0 0 0.050 1 6 0 0 1\n")
        f.write(f"mol color ColorID {vmd_color_ids[color2]}\n")
        f.write("mol selection {all}\n")
        f.write(f"mol material Glass3\n")
        f.write(f"mol addrep top\n")
        f.write(f"mol selupdate 1 top 0\n")
        f.write(f"mol colupdate 1 top 0\n")
        f.write(f"mol scaleminmax top 1 0.000000 0.000000\n")
        f.write(f"mol smoothrep top 1 0\n")
        f.write("mol drawframes top 1 {now}\n")
        
def VMD_Orbital_Trajectory(molden_location,orbital_number,color1="blue",color2="red"):
    ORDERED_MOLDEN_FILELIST=[]
    startdir = os.getcwd()
    os.chdir(molden_location)
    file_list = G("*.molden.*")
    file_list.sort(key=os.path.getmtime)
    n_moldens = len(file_list)
    os.makedirs("phased_molden_files/",exist_ok=True)
    os.system(f"cp *.molden.* phased_molden_files/")
    for i in range(n_moldens-1):
        OLD_MOL_FILE=file_list[i]
        NEW_MOL_FILE="phased_molden_files/"+file_list[i+1]
        if i > 0:
            OLD_MOL_FILE="phased_molden_files/"+file_list[i]
        print(OLD_MOL_FILE)
        print(NEW_MOL_FILE)
        print("----")
        sign_array = Get_Sign_Array(OLD_MOL_FILE,NEW_MOL_FILE)
        Write_Phase_Swapped_Molden(NEW_MOL_FILE,sign_array)
        ORDERED_MOLDEN_FILELIST.append(OLD_MOL_FILE)
    ORDERED_MOLDEN_FILELIST.append(NEW_MOL_FILE)
    VMDInitialize("orbital_trajectory.vmd")
    Write_VMD_Command_File(ORDERED_MOLDEN_FILELIST,orbital_number,color1,color2)
    VMDTrajectoryMovie("orbital_trajectory.vmd","TachyonLOptiXInternal")
    os.chdir(startdir)
    
def VMDTrajectoryMovie(inputfilename,renderer):
    with open(inputfilename,"a") as f:
        f.write("proc make_trajectory_movie {} {\n")
        f.write("\t# get the number of frames in the movie\n")
        f.write("\tset num [molinfo top get numframes]\n")
        f.write("\t# loop through the frames\n")
        f.write("\tfor {set i 0} {$i < $num} {incr i 1} {\n")
        f.write("\t\t# go to the given frame\n")
        f.write("\t\tanimate goto $i\n")
        f.write("\t\t\t\t# for the display to update\n")
        f.write("\t\t\t\tdisplay update\n")
        f.write("\t\t# take the picture\n")
        f.write(f"\t\tset filename temp.[format \"%05d\" [expr $i/1]].tga\n")
        f.write(f"\t\trender {renderer} $filename\n")
        f.write("\t}\n")
        f.write("}\n")
        f.write("make_trajectory_movie\n")
        f.write("quit\n")
        f.close()
    return None
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--path",dest="path",help="/path/to/directory/with/molden/files",required=True)
    parser.add_argument("-n","--number",dest="number",help="Orbital number desired for visualization",required=True)
    parser.add_argument("--color1",dest="color1",help="VMD-available color for first phase of orbitals",default="blue")
    parser.add_argument("--color2",dest="color2",help="VMD-available color for second phase of orbitals",default="red")
    args = parser.parse_args()
    VMD_Orbital_Trajectory(args.path,int(args.number),color1=args.color1,color2=args.color2)
