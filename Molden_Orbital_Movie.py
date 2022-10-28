#!/usr/bin/env python3
##################################################
# Program Name:         Molden_Orbital_Movie
# Original Author:      Mark A. Hix
# Date Last Updated:    2022 October 26
# Version:              1.0
# Software Required:    Python3,VMD
# -- Special thanks to Prof. Todd Martinez and Dr. Jacqueline R. Shea 
# for help on the math and Molden files.
##################################################

## Working from standard python libraries only
import numpy as np
from glob import glob as G
import os


def Get_Molecular_Orbitals(moldenfile):
    # Process Molden file to obtain the coefficients for all molecular orbitals
    # Probably could do this in a way that only iterates through the file once, 
    # but this works fine as-is and it's not worth the time to optimize.
    
    ## Read Moldenfile into memory
    with open(moldenfile) as f:
        lines = f.readlines()

    ## Get number of orbitals by counting lines with "Ene=" in them, 
    ## since each orbital chunk begins with that information.
    ene_lines = []
    for i,line in enumerate(lines):
        if "Ene=" in line:
            ene_lines.append(i)
    n_orbs   = len(ene_lines)

    ## Run through orbital chunks, obtain coefficients for each atomic orbital in it.
    ## Add orbital chunk to array of orbital chunks
    mol_orbs = []
    for i in range(n_orbs-1):
        orb = []
        for line in lines[ene_lines[i]+3:ene_lines[i+1]]:
            orb.append(float(line.split()[-1]))
        mol_orbs.append(orb)
    orb = []
    ## Get last orbital chunk because the slicing is different and
    ## I don't care to make it more elegant at the moment.  This works fine.
    for line in lines[ene_lines[n_orbs-1]+3:]:
        orb.append(float(line.split()[-1]))
    mol_orbs.append(orb)
    return mol_orbs

def Get_Sign_Array(OLD_MOL_FILE,NEW_MOL_FILE):
    # Compare two orbital sets, assumed to be from sequential QM/MM-MD
    # steps, and take the dot-product of each orbital to determine if the
    # phase has been swapped from one to the next.
    
    ## Process both molden files
    old_mol_orbs = Get_Molecular_Orbitals(OLD_MOL_FILE)
    new_mol_orbs = Get_Molecular_Orbitals(NEW_MOL_FILE)
    
    ## initialize list of orbital-phase-signs
    sign_array = []
    
    ## iterate through all orbital-pairs, get the sign, append to the list.
    for i in range(len(old_mol_orbs)):
        sign_array.append(np.sign(np.dot(old_mol_orbs[i],new_mol_orbs[i])))
    return sign_array

def Write_Phase_Swapped_Molden(moldenfile,sign_array):
    # Rewrite Molden file with signs swapped on coefficients of designated MOs
    # based on value of sign_array
    
    ## Get the original file information
    with open(moldenfile) as f:
        newlines = f.readlines()
    
    ## Rewrite that file with the modified information
    with open(moldenfile,"w") as f:
        ## Same stuff above the MOs as before.  No changes here.
        for i,line in enumerate(newlines):
            f.write(line)
            if "[MO]" in line:
                break
        ## Now that we're done with pre-MO stuff, we can add our modified coefficients.
        start_mo_lines = i+1
        ## a lazy counter, I know, but it makes indexing *way* easier for me.
        ene_lines_encountered = -1
        ## run through remaining lines from the file, checking the contents and adjusting.
        for line in newlines[start_mo_lines:]:
            if "Ene=" in line:
                ene_lines_encountered += 1
                f.write(line)
            elif "Spin=" in line:
                f.write(line)
            elif "Occup=" in line:
                f.write(line)
            else:
                ## Here's where we actually change signs on any orbital with a negative.
                ## split into AO index and coefficient
                templine = line.split()
                if sign_array[ene_lines_encountered] == -1:
                    ## change sign of coefficients
                    if "-" in templine[1]:
                        templine[1] = templine[1].replace("-"," ")
                    else:
                        templine[1] = "-"+templine[1].strip()
                ## rebuild the line based on (un)modified coefficient values, maintain formatting.
                newline = f"{templine[0]:>5}{float(templine[1]):>11.05f}\n"
                ## output (un)modified line to Molden file.
                f.write(newline)

def VMDInitialize(inputfilename):
    # Initialize VMD environment for some standard behaviors.  
    # End Users can modify this as much as they want, but Tcl isn't fun.
    # Values below should be fairly self-explanatory, see VMD manual for details.
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

def Write_VMD_Command_File(vmdcommandfile,ORDERED_MOLDEN_FILELIST,orbital_number,color1,color2):
    # With the list of Molden files, a specified orbital number, and two color keywords, 
    # we generate the VMD input commands.
    
    ## VMD colorID dictionary
    vmd_color_ids={"blue":0,"red":1,"gray":2,"orange":3,"yellow":4,
                   "tan":5,"silver":6,"green":7,"white":8,"pink":9,
                   "cyan":10,"purple":11,"lime":12,"mauve":13,"ochre":14,
                   "iceblue":15,"black":16}
    ## check if color1 and color2 values will be recognized by VMD.  If not, you get defaults.
    if color1 not in [key for key in vmd_color_ids.keys()]:
        print(f"Color {color1} not found in VMD dictionary.  Defaulting to blue for color 1.")
        color1 = "blue"
    if color2 not in [key for key in vmd_color_ids.keys()]:
        print(f"Color {color2} not found in VMD dictionary.  Defaulting to red for color 2.")
        color2 = "red"
        
    with open(vmdcommandfile,"a") as f:
        num_moldens_loaded = 0
        ## Iterate through ordered molden filelist
        for moldenfile in ORDERED_MOLDEN_FILELIST:
            if num_moldens_loaded == 0:
                ## The first datafile in the molecule has to use 'mol new', while the rest can use 'mol addfile'
                f.write(f"mol new {moldenfile} type molden first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n")
            else:
                f.write(f"mol addfile {moldenfile} type molden first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n")
            num_moldens_loaded += 1
        
        ## Add representation for the molecule contained in the Molden file.
        ## I like licorice because it's clear but doesn't overwhelm the viewer 
        ## when the focus is going to be on the orbitals.
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
        
        ## Add representation for + phase of orbital, using color1.
        ## Using "Glass3" material because it offers good transparency
        ## while still showing the color choice
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
        
        ## Add representation for - phase of orbital, using color2.
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
        
def VMD_Orbital_Trajectory(molden_location,orbital_number,vmdcommandfile,x,y,z,enlarge,color1="blue",color2="red"):
    # Main Job Function.
    # Takes the location of the molden files (assumed to be in the same directory, such as a scratch folder).
    # desired orbital number, and color choices (defaulted to blue/red).
    ORDERED_MOLDEN_FILELIST=[]
    startdir = os.getcwd()
    
    ## Move into molden directory
    os.chdir(molden_location)
    
    ## Get list of molden files and sort by time, get total count.
    ## (assumes the files are output sequentially during QM/MM-MD)
    file_list = G("*.molden.*")
    file_list.sort(key=os.path.getmtime)
    n_moldens = len(file_list)
    
    ## Probably should make this a function parameter, but eh.
    ## Make a folder for the modified Molden files so we don't 
    ## actually overwrite the originals.
    output_folder = "phased_molden_files/" 
    os.makedirs(output_folder,exist_ok=True)
    
    ## Copy original moldens into new folder for safer processing.
    ## I feel like this is overly paranoid, but I'd rather be too 
    ## careful than not careful enough.
    os.system(f"cp *.molden.* {output_folder}")
    for i in range(n_moldens-1):
        ## Run through each pair of molden files, writing modified files as we go.
        OLD_MOL_FILE=output_folder+file_list[i]
        NEW_MOL_FILE=output_folder+file_list[i+1]
        sign_array = Get_Sign_Array(OLD_MOL_FILE,NEW_MOL_FILE)
        Write_Phase_Swapped_Molden(NEW_MOL_FILE,sign_array)
        ORDERED_MOLDEN_FILELIST.append(OLD_MOL_FILE)
    ## add the final file to the list.
    ORDERED_MOLDEN_FILELIST.append(NEW_MOL_FILE)
    
    ## Initialize the VMD input file with the display settings.
    VMDInitialize(vmdcommandfile)
    
    ## Add the molden files and representations
    Write_VMD_Command_File(vmdcommandfile,ORDERED_MOLDEN_FILELIST,orbital_number,color1,color2)
    
    ## Add the moviemaker VMD function and call it so that VMD 
    ## automatically renders the video on completion of file I/O.
    ## Probably not ideal since the orientation of the molecule in 
    ## space won't necessarily be what you want...
    VMDRotateStructure(vmdcommandfile,x,y,z,enlarge)
    VMDTrajectoryMovie(vmdcommandfile,"TachyonLOptiXInternal")
    os.chdir(startdir)

def VMDRotateStructure(inputfilename,x,y,z,enlarge):
    with open(inputfilename,"a") as f:
        f.write(f"rotate x by {x}\n")
        f.write(f"rotate y by {y}\n")
        f.write(f"rotate z by {z}\n")
        f.write(f"scale by {enlarge}\n")

def VMDTrajectoryMovie(inputfilename,renderer):
    # Tcl code for VMD to generate a set of frames from a trajectory.
    with open(inputfilename,"a") as f:
        f.write("## Uncomment lines below to change initial rotation/orientation of molecule before rendering\n")
        f.write(f"## then rerun this file with 'vmd -e {inputfilename}'\n")
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
        f.write(f"\t\tset filename VMD_Images/temp.[format \"%05d\" [expr $i/1]].tga\n")
        f.write(f"\t\trender {renderer} $filename\n")
        f.write("\t}\n")
        f.write("}\n")
        f.write("make_trajectory_movie\n")
        f.close()
    return None
    
if __name__ == "__main__":
    # Since argparse is beefy, only load it up if it's actually getting used, 
    # such as when this script is called from the command line.
    import argparse
    import subprocess
    ## construct the argument parser and add arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("-p","--path",dest="path",help="/path/to/directory/with/molden/files",required=True)
    parser.add_argument("-n","--number",dest="number",help="Orbital number desired for visualization",required=True)
    parser.add_argument("--color1",dest="color1",help="VMD-available color for first phase of orbitals",default="blue")
    parser.add_argument("--color2",dest="color2",help="VMD-available color for second phase of orbitals",default="red")
    parser.add_argument("-v","--vmd",dest="vmd",help="VMD command file output",default="orbital_trajectory.vmd")
    parser.add_argument("-x",dest="x_rotate",help="angle in degrees to rotate structure about x-axis before rendering.",default=0)
    parser.add_argument("-y",dest="y_rotate",help="angle in degrees to rotate structure about y-axis before rendering.",default=0)
    parser.add_argument("-z",dest="z_rotate",help="angle in degrees to rotate structure about z-axis before rendering.",default=0)
    parser.add_argument("-e","--enlarge",dest="enlarge",help="Scaling multiplier (makes molecule larger/smaller in view)",default=1.0)
    parser.add_argument("-o","--output",dest="output",help="OutputMovie.mp4",default="./orbitals.mp4")

    ## Parse those args.
    args = parser.parse_args()
    
    ## Call the main job
    VMD_Orbital_Trajectory(os.path.abspath(args.path),int(args.number),args.vmd,args.x_rotate,args.y_rotate,args.z_rotate,args.enlarge,color1=args.color1,color2=args.color2)
    with open(args.vmd,"a") as f:
        f.write("quit\n")
    subprocess.call("mkdir -p VMD_Images/",shell=True)
    subprocess.call(f"vmd -e {args.vmd}",shell=True)
    outputmovie = os.path.abspath(args.output)
    os.chdir("VMD_Images/")
    subprocess.call("for i in $(ls temp.*.tga); do convert $i ${i%.tga}.png; rm $i; done",shell=True)
    subprocess.call(f"ffmpeg -y -r 30 -f image2 -s 1920x1080 -i temp.%05d.png -vf \"scale=trunc(iw/2)*2:trunc(ih/2)*2\" -vcodec libx264 -crf 5 -pix_fmt yuv420p {outputmovie}",shell=True)
    os.chdir("../")
