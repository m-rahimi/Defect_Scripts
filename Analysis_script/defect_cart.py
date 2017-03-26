#! /usr/bin/env python
# Mohammad Rahimi 4 Dec 2014
# Calculate scaler order, nematic field and local density in cartesian coordinate
# input psf and dcd from NAMD
# for more than one trj line 74 should be modified
# output file can be visulazie by gnuplot
# the gnupolt file was copied at the end of the code

from MDAnalysis import *
from math import *
import numpy as np
from scipy import interpolate
import time
import sys
import argparse
import os

parser = argparse.ArgumentParser(description="Read trajectory and coordinate file")

parser.add_argument("-t", 
                    action="store", nargs='?', 
                    required=True, dest="traj", 
                    help="specifies an .dcd file created using the '-pbc mol' option")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-nx",
                    action="store",nargs='?',
                    required=True, dest="nslabx", type=int,
                    help="Number of slabs")
parser.add_argument("-ny",
                    action="store",nargs='?',
                    required=True, dest="nslaby", type=int,
                    help="Number of slabs")
parser.add_argument("-vx",
                    action="store",nargs='?',
                    required=True, dest="vslabx", type=int,
                    help="Number of slabs")
parser.add_argument("-vy",
                    action="store",nargs='?',
                    required=True, dest="vslaby", type=int,
                    help="Number of slabs")
parser.add_argument("-ti",
                    action="store",nargs='?', default=int('20'),
                    required=False, dest="time", type=int,
                    help="time interval (ps)")
parser.add_argument("-o",
                    action="store", nargs='?', default="Output", 
                    required=False, dest="output",
                    help="output filename")
parser.add_argument("-b",
                    action="store", nargs='?', default=int('0'), 
                    required=False, dest="start_frame", type=int,
                    help="First frame to start reading trajectory")
parser.add_argument("-e",
                    action="store", nargs='?', default=float('inf'), 
                    required=False, dest="end_frame", type=int,
                    help="Last frame to stop readining trajectory")


args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
trj_time = args.time
nslabx = args.nslabx
nslaby = args.nslaby
vslabx = args.vslabx
vslaby = args.vslaby
output_filename = args.output
start_frame = args.start_frame
end_frame = args.end_frame

atom_selection_N = "name NY1"
atom_selection_C = "name CA12"
atom_selection_N5 = "resname 5CB and name NY1"
atom_selection_N8 = "resname 8CB and name NY1"

time1 = time.time()
output_file = open(output_filename + "_scaler","w")
output_file2 = open(output_filename + "_vector" ,"w")

def vector(v1, v2):
    return [[v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]]
def position(v1,v2):
    return (v1+v2)/2.0
def anint(a):
    if a > 0 :
       an = int(a+0.5)
    else :
       an = int(a-0.5)
    return int(an)
# Read the trj and gro file
#u = Universe("mix.psf", ["65.dcd", "66.dcd", "67.dcd", "68.dcd"])
u = Universe(psf_filename, traj_filename)

# Obtain initial information form gro and trj files
natoms = len(u.selectAtoms(atom_selection_N))
natoms5 = len(u.selectAtoms(atom_selection_N5))
natoms8 = len(u.selectAtoms(atom_selection_N8))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of selected atoms = " + str(natoms)
print "Number of 5CB molecules = " + str(natoms5)
print "Number of 8CB molecules = " + str(natoms8)

delta = (natoms/nslabx/nslaby)
deltav = (natoms/vslabx/vslaby)
print "Number of atoms in each slab is around " + str(delta)
print "Number of atoms in each slab is around " + str(deltav)
num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
   end_frame = num_frames
print "read frame " + str(start_frame) + " to " + str(end_frame)

# Define arrays
binsx = np.zeros((nslabx+1),dtype='float')
binsy = np.zeros((nslaby+1),dtype='float')
vbinsx = np.zeros((vslabx+1),dtype='float')
vbinsy = np.zeros((vslaby+1),dtype='float')

Q_slab = np.empty((nslabx,nslaby,9),dtype='float')
dens = np.zeros((nslabx,nslaby),dtype='int')
dens5 = np.zeros((nslabx,nslaby),dtype='int')
dens8 = np.zeros((nslabx,nslaby),dtype='int')
dens_f = np.zeros((nslabx,nslaby),dtype='float')
Scalar = np.empty((nslabx,nslaby),dtype='float')
Biaxial = np.empty((nslabx,nslaby),dtype='float')

Q_slabv = np.empty((vslabx,vslaby,9),dtype='float')
densv = np.zeros((vslabx,vslaby),dtype='int')
densv_f = np.zeros((vslabx,vslaby),dtype='int')
Scalarv = np.empty((vslabx,vslaby),dtype='float')

vector_x = np.zeros((vslabx,vslaby),dtype='float')
vector_y = np.zeros((vslabx,vslaby),dtype='float')
vector_z = np.zeros((vslabx,vslaby),dtype='float')

Qt = np.zeros((3,3),dtype='float')
I = np.matrix([[1.,0,0],[0,1.,0],[0,0,1.]])

Scalar.fill(0)
Biaxial.fill(0)
Scalarv.fill(0)
Q_slab.fill(0)
Q_slabv.fill(0)
nframes = 0

# Loops in the trj
for curr_frame in xrange(0, num_frames) :
    if curr_frame != 0 :
        trj = u.trajectory.next()
    else :
        trj = u.trajectory[0]
    curr_time = (curr_frame + 1) * trj_time
    if curr_frame < start_frame :
        continue
    if curr_frame > end_frame :
        break
    else :
        nframes += 1
        print "Reading frame " + str(curr_frame) + " at " + str(curr_time) + " ps"
        coor_N = u.selectAtoms(atom_selection_N).coordinates()
        coor_C = u.selectAtoms(atom_selection_C).coordinates()
        
        coor_N5 = u.selectAtoms(atom_selection_N5).coordinates()
        coor_N8 = u.selectAtoms(atom_selection_N8).coordinates()
        
        if len(coor_N) != len(coor_C) :
           print >>sys.stderr, "Error: number of atoms in each group does not match"
           sys.exit(1)
        
        #Calculate the opsition of each slab to have almost same number of atoms in each slab
        #The method is correct if the falctuation of volume is very small
        box = u.dimensions[0:3]
        box2 = box/2
        dx = box[0]/nslabx
        dy = box[1]/nslaby
        volume = (dx * dy * box[2]) / 1000  # nm^3
        for i in range(nslabx+1):
            binsx[i] = i*dx - box2[0]
        for i in range(nslaby+1):
            binsy[i] = i*dy - box2[1]
 
        dxv = box[0]/vslabx
        dyv = box[1]/vslaby
        volumev = (dxv * dyv * box[2]) / 1000  # nm^3
        for i in range(vslabx+1):
            vbinsx[i] = i*dxv - box2[0]
        for i in range(vslaby+1):
            vbinsy[i] = i*dyv - box2[1]
       
#        # Calculate Q tensor for atoms
#        # Calculate Q tensor for each slab by summing Q of atomes located in the slab
#        # Calculate global Q 
#        # Calculate number of atoms in each slab
        j_5CB = 0
        j_8CB = 0
        Qt.fill(0)
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]

            px = position(loc_C[0],loc_N[0])
            py = position(loc_C[1],loc_N[1])

            # apply periodic boundary conditions
            px = px - box[0]*anint(px/box[0]) + box2[0]
            py = py - box[1]*anint(py/box[1]) + box2[1]
            
            vv = vector(loc_C,loc_N)
            vv = vv/np.linalg.norm(vv)
            Q = (3*np.dot(vv.T,vv) - I)/2.0 # factor 3/2 means that biggest eigvalu = nematic order
            Qt += Q
            Q = np.reshape(Q,(9))
           
            nx = int(px/dx) 
            ny = int(py/dy)
            if nx == nslabx : nx -= 1
            if nx < 0 :       nx += 1
            if ny == nslaby : ny -= 1
            if ny < 0 :       ny += 1
            dens[nx,ny] += 1
            Q_slab[nx,ny,:] = Q + Q_slab[nx,ny,:]
            
            if (j_5CB<natoms5 and np.allclose(loc_N, coor_N5[j_5CB])) : 
                j_5CB += 1
                dens5[nx,ny] += 1
            if (j_8CB<natoms8 and np.allclose(loc_N, coor_N8[j_8CB])) : 
                j_8CB += 1
                dens8[nx,ny] += 1

            nx = int(px/dxv) 
            ny = int(py/dyv)
            if nx == vslabx : nx -= 1
            if nx < 0 :       nx += 1
            if ny == vslaby : ny -= 1
            if ny < 0 :       ny += 1
            densv[nx,ny] += 1
            Q_slabv[nx,ny,:] = Q + Q_slabv[nx,ny,:]

        Qt = Qt/natoms  # Normalize global Q
        e_value, e_vector = np.linalg.eigh(Qt)  # Calculate global order parameter
        print "time " + str(curr_time) + " Scalar " + str(np.amax(e_value))

print "Total number of frames " + str(nframes)

for x in range(nslabx):
  for y in range(nslaby):
      dens_f[x,y] = float(dens[x,y]) / nframes
      if dens_f[x,y]>delta*0.4 :
         Q = Q_slab[x,y,:] / dens[x,y]  
         Q = np.reshape(Q,(3,3))
         e_value, e_vector = np.linalg.eigh(Q) # Calculate local order parameter
         e_value = e_value[np.argsort(e_value)] 
         Scalar[x,y] = e_value[2]   # Larger eigenvalue is S
         Biaxial[x,y] = 2*e_value[1] + e_value[2]
      else : 
         Scalar[x,y] = 0.56

for x in range(vslabx):
  for y in range(vslaby):
      densv_f[x,y] = float(densv[x,y]) / nframes
      if densv_f[x,y]>deltav*0.4 :
         Q = Q_slabv[x,y,:] / densv[x,y]  
         Q = np.reshape(Q,(3,3))
         e_value, e_vector = np.linalg.eigh(Q) # Calculate local order parameter
         Scalarv[x,y] = np.amax(e_value)   # Larger eigenvalue is S
         vector_x[x,y] = e_vector[0,np.argmax(e_value)]
         vector_y[x,y] = e_vector[1,np.argmax(e_value)]
         vector_z[x,y] = e_vector[2,np.argmax(e_value)]
      else : 
         Scalarv[x,y] = 0.56

##write Scalar
output_file.write("# X   Y   Scaler   Biaxial    Ndens    dens(nm-3)      dens5   dens8 \n")
for x in range(nslabx):
  for y in range(nslaby):
    if dens[x,y]!=0 :
        ro85 = (dens8[x,y]-dens5[x,y])/(dens8[x,y]+dens5[x,y])
    else :
        ro85 = 0
    output_file.write('%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n' % ((binsx[x]+binsx[x+1])/2,(binsy[y]+binsy[y+1])/2,Scalar[x,y],Biaxial[x,y],dens[x,y]/nframes,dens[x,y]/nframes/volume,dens5[x,y]/nframes/volume,dens8[x,y]/nframes/volume,ro85))
  output_file.write('\n')
##write vectors
output_file2.write("# X   Y   Z  Scaler   Ndens   VX   Vy   Vz\n")
zz = 0.0
for x in range(vslabx):
  for y in range(vslaby):
      output_file2.write('%5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f \n' % (((vbinsx[x+1]+vbinsx[x])/2.0),((vbinsy[y+1]+vbinsy[y])/2.0),zz,Scalarv[x,y],densv_f[x,y],vector_x[x,y],vector_y[x,y],vector_z[x,y]))
  output_file2.write('\n')
#
#
time2 = time.time()
S_time = time2-time1
print "Simulation time = " + str(S_time)
output_file.close()
output_file2.close()


##################################################
## Set terminal and output
#set terminal pngcairo size 1000,1000 enhanced font 'Verdana,15'
##set terminal  postscript enhanced 18
#set output 'P.png'
## Set various features of the plot
#set pm3d
#unset surface  # don't need surfaces
#set view map
##set contour
##set cntrparam bspline #cubicspline  # smooth out the lines
##set cntrparam levels -2    # sets the num of contour lines
#set pm3d interpolate 20,20 # interpolate the color
#
## Set a nice color palette
#set palette defined ( 0 "green", 1 "blue", 2 "red") 
#
##set title 'Free Energy Profile ({/Symbol D}F)'
## Axes
#set xlabel 'x (A)' font 'Verdana,20' #rotate by 0
#set ylabel 'y (A)' rotate by 90 font 'Verdana,20'
#set cblabel 'Density' rotate by 0 offset -8.,14.5,0 font 'Verdana,20'
#
#set tics nomirror
#set tics out
#set xtics 26.121 font 'Verdana,13'
#set ytics 26.121 font 'Verdana,13'
#set mxtics 3
#set mytics 3
#set style line 12 lc rgb 'blue' lt 1 lw 2
#set style line 13 lc rgb 'yellow' lt 1 lw 1
#set grid xtics ytics mxtics mytics #ls 12, ls 13
##set grid mxtics mytics
#
#
#set xrange [-78.363:78.363]
#set yrange [-78.363:78.363]
#
##set cbrange [30:80]
## Now plot
#set multiplot
#splot 'Output_scaler' using 1:2:4 notitle with lines lt 1, 'Output_vector' u ($1-5*$6):($2-5*$7):3:(10*$6):(10*$7):8 notitle with vector nohead lt 1 lc 0 
#
#r = 21
#zz = 0
#set parametric
#set vrange [0:2*pi]
#fx(v) = r*cos(v)
#fy(v) = r*sin(v)
#unset pm3d
#set surface
#
#splot fx(v),fy(v),zz
#
#unset multiplot
#pause -1
