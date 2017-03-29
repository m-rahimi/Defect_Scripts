#! /usr/bin/env python

from MDAnalysis import *
from math import *
import numpy as np
import time
import sys
import argparse
import os
import matplotlib.pyplot as plt
from subprocess import *
from pylab import *


parser = argparse.ArgumentParser(description="Read trajectory and coordinate file")

parser.add_argument("-t", 
                    action="store", nargs='?', default="",
                    required=False, dest="traj", 
                    help="specifies a dcd file created using NAMD package")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-d",
                    action="store",nargs='?',
                    required=True, dest="dz", type=float,
                    help="Grid thickness")
parser.add_argument("-ti",
                    action="store",nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")                                        
parser.add_argument("-o",
                    action="store", nargs='?', 
                    required=True, dest="output1",
                    help="output filename for ScalarZ")
parser.add_argument("-b",
                    action="store", nargs='?', default=int('0'), 
                    required=False, dest="start_time", type=int,
                    help="First time to start reading trajectory")
parser.add_argument("-e",
                    action="store", nargs='?', default=float('inf'), 
                    required=False, dest="end_time", type=int,
                    help="Last time to stop readining trajectory")
parser.add_argument("-f", 
                    action="store", nargs='?',
                    required=False, dest="firstFileName", type=int,
                    help="The name of the first trajectory file before dot")
parser.add_argument("-l", 
                    action="store", nargs='?',
                    required=False, dest="lastFileName", type=int,
                    help="The name of the last trajectory file before dot")
parser.add_argument("-s",
                    action="store", nargs='?', 
                    required=True, dest="select",
                    help="atoms selection argument")


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
dz = args.dz
output_filename = args.output1
start_frame = args.start_time
end_frame = args.end_time
firstFile = args.firstFileName
lastFile = args.lastFileName
atom_selection = args.select


time1 = time.time()
output_file = open(output_filename,"w")

# Read the trj and psf file
if traj_filename == "" :
	list = [str(firstFile)+".dcd"]
	for i in xrange(firstFile+1,lastFile+1) :
		list.append(str(i)+".dcd")
	print "List of Trajectory files:" + str(list)
	u = Universe(psf_filename, list)
else :
	u = Universe(psf_filename, traj_filename)

# Obtain initial information form psf and trj files
natoms = len(u.selectAtoms(atom_selection))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)

# Define Slabs
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
  end_frame = num_frames
print "reading frame " + str(start_frame) + " to " + str(end_frame)

# Deine arrays

dens_f = np.empty((0),dtype='float')



pvector = np.zeros((natoms),dtype='float')
nframes = 0
for curr_frame in xrange(start_frame, end_frame) :
    if curr_frame != 0 :
        trj = u.trajectory.next()
    else :
        trj = u.trajectory[0]
    curr_time = (curr_frame + 1) * trj_time
    if curr_frame > end_frame :
        break
    else :
        nframes += 1
        print "Reading frame " + str(curr_frame) + " at " + str(curr_time) + " ps"
        coor = u.selectAtoms(atom_selection).coordinates()
 
        # Applay PBC 
        #box = u.dimensions[0:3]
        
        #Calculate the position of each slab to have almost same number of atoms in each slab
        if curr_frame == start_frame :
            for i in xrange(0,len(coor)) :
                loc = coor[i]
                pvector[i] = loc[2]
            
            pvector = np.sort(pvector)
            highZ = int(pvector[natoms-1]) + 5
            lowZ = int(pvector[0])
            nslab = int((highZ - lowZ)/dz)
            bins = np.zeros((nslab+1),dtype='float')
            bins[nslab] = highZ
            bins[0] = lowZ
            slab_ind = nslab - 1
            dens = np.zeros((nslab),dtype='int')
            for i in xrange(1,nslab):
		        bins[i] = bins[i-1] + dz

        dens.fill(0)
        for i in xrange(0,len(coor)) :
            loc = coor[i]
            p = loc[2]
            x = min(int((p-bins[0])/dz),slab_ind)
            dens[x] += 1

        dens_f = np.append(dens_f,dens,0)

		

dens_f = np.reshape(dens_f,( len(dens_f)/nslab,nslab))
dens_average =  np.average(dens_f,0)
dens_deviation =  np.std(dens_f,0)       
dens_average = dens_average/len(coor) 
p0 = (bins[0]+bins[1])/2
output_file.write('%5.3f %10.7f\n' % ( bins[0] ,0))
for x in range(nslab):
    p = (bins[x]+bins[x+1])/2
    d = (bins[x+1]-bins[x])
    output_file.write('%5.3f %10.7f\n' % ( p ,dens_average[x]))
output_file.write('%5.3f %10.7f\n' % ( bins[nslab] ,0))


time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file.close()

##Plot the distribution results
#outFile = open(str(output_filename), 'r')
#binCenter=[]
#P = []
#for line in outFile:
#	line = line.strip()
#	if line != "":
#		columns = line.split()
#		if '#' in columns[0]:
#			check='yes'
#		else:
#			firstCol = float(columns[0])
#			binCenter.append(firstCol)
#			secCol = float(columns[1])
#			P.append(secCol) 
##Plot the cosZ data
#plt.rcParams['axes.linewidth'] = 1.5
##plt.xticks(rotation=25)
#plt.ylabel(r'$\mathrm{P(z)}$', fontsize=27)
#plt.xlabel(r'$\mathrm{z\/(\AA)}$', fontsize=27)	
#plt.tick_params(which='both', width=2)		
##plt.ylim((0,0.1))
##plt.ylim((binCenter[0],binCenter[len(binCenter)-1]))       
##plt.yticks(weight='bold')
##plt.xticks(weight='bold')
#plt.yticks(fontsize=17)
#plt.xticks(fontsize=17)
##plt.gcf().set_size_inches(5,10)
#plt.plot(binCenter, P, marker='o', markeredgecolor='b')
#plt.savefig('../Analysis/plot_dist'+atom_selection+'_.png', dpi=600, facecolor='w', edgecolor='w',
#        orientation='portrait', papertype='letter', format='png',
#        transparent=True, bbox_inches='tight', pad_inches=0.1,
#        frameon=None)
