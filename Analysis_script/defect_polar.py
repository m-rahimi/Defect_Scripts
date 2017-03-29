#! /usr/bin/env python
# 9 Dec 2014
# calculate density number, scaler order and nematic director 
# input dcd, psf and input files
# input file :: number of bin; number of molecule in a bin; angle
# output :: output.vtk
# Triangle program must be run at line 257
# http://www.cs.cmu.edu/~quake/triangle.html
#############################################
# input
# Number of bins, Number atoms in each bin, angle
# For example: 20 30 10
# Number of layer = number of line in input file

from MDAnalysis import *
from math import *
import numpy as np
import time
import sys
import argparse
import os
import random
from subprocess import *

parser = argparse.ArgumentParser(description="Make triangle")

parser.add_argument("-t", 
                    action="store", nargs='?', 
                    required=True, dest="traj", 
                    help="specifies an .dcd file created using the '-pbc mol' option")
parser.add_argument("-p",
                    action="store", nargs='?',
                    required=True, dest="psf",
                    help="path of a sample psf file")
parser.add_argument("-l", 
                    action="store", nargs='?', type=int, 
                    required=True, dest="layer", 
                    help="")
parser.add_argument("-s",
                    action="store", nargs='?', default=int('0'), 
                    required=False, dest="start_frame", type=int,
                    help="First frame to start reading trajectory")
parser.add_argument("-e",
                    action="store", nargs='?', default=float('inf'), 
                    required=False, dest="end_frame", type=int,
                    help="Last frame to stop readining trajectory")
parser.add_argument("-time",
                    action="store", nargs='?', default=float('20'), 
                    required=False, dest="time", type=float,
                    help="Time between to frames (ps)")

args = parser.parse_args()
traj_filename = args.traj
psf_filename = args.psf
Nlayer = args.layer
start_frame = args.start_frame
end_frame = args.end_frame
trj_time = args.time
atom_selection_N = "name NY1"
atom_selection_C = "name CA12"
atom_selection_N5 = "resname 5CB and name NY1"
atom_selection_N8 = "resname 8CB and name NY1"

time1 = time.time()

Nbin = np.zeros((Nlayer),dtype='int')
Nb_atoms = np.zeros((Nlayer),dtype='int')
angle = np.zeros((Nlayer),dtype='float')
Rbin = np.zeros((Nlayer+1),dtype='float')
Nl_atoms = np.zeros((Nlayer),dtype='float')
dtheta = np.zeros((Nlayer),dtype='float')
output_file = open("85CB.vtk","w")
f1=open('input', 'r')
i=0
for i in range(0,Nlayer):
    line = f1.readline()
    line = line.strip()
    columns = line.split()
    Nbin[i] = float(columns[0])
    Nb_atoms[i] = float(columns[1])
    angle[i] = float(columns[2])   #*pi/180
    Nl_atoms[i] = Nbin[i] * Nb_atoms[i]
    dtheta[i] = 360 / Nbin[i] 
f1.close()

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
u = Universe(psf_filename, traj_filename)
#u = Universe("mix.psf", ["65.dcd", "66.dcd", "67.dcd", "68.dcd"])

# Obtain initial information form gro and trj files
natoms = len(u.selectAtoms(atom_selection_N))
natoms5 = len(u.selectAtoms(atom_selection_N5))
natoms8 = len(u.selectAtoms(atom_selection_N8))
print "Total number of atoms = " + str(len(u.atoms))
print "Number of total molecules = " + str(natoms)
print "Number of 5CB molecules = " + str(natoms5)
print "Number of 8CB molecules = " + str(natoms8)

num_frames = u.trajectory.numframes
print "number of frames " + str(num_frames) 
if end_frame > num_frames :
   end_frame = num_frames
print "read frame " + str(start_frame) + " to " + str(end_frame)

print "Number of layers = " + str(Nlayer)
print "Number of bins = " + str(Nbin)


# calculte density and Q
Tbin = 0
#Max_atoms = 0
Sbin = np.zeros((Nlayer),dtype='int')
for i in range(0,Nlayer):
    Sbin[i] = Tbin
    Tbin += Nbin[i]
#    if Max_atoms < Nbin[i]*Natoms[i] : Max_atoms = Nbin[i]*Natoms[i]
#Max_bin = Nbin[Nlayer-1]
dens = np.zeros((Tbin),dtype='int')
dens5 = np.zeros((Tbin),dtype='int')
dens8 = np.zeros((Tbin),dtype='int')
Q_slab = np.zeros((Tbin,9),dtype='float')
volume = np.zeros((Tbin),dtype='float')
Scalar = np.zeros((Tbin),dtype='float')
Biaxial = np.zeros((Tbin),dtype='float')
vector_x = np.zeros((Tbin),dtype='float')
vector_y = np.zeros((Tbin),dtype='float')
vector_z = np.zeros((Tbin),dtype='float')

pos = np.zeros((natoms),dtype='float') 
I = np.matrix([[1.,0,0],[0,1.,0],[0,0,1.]])

nframes = 0
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
        if (curr_frame%10)==0 :
           print "Reading frame " + str(curr_frame) + " at " + str(curr_time) + " ps"
           print "calculating frame " + str(nframes) + " at " + str(nframes*trj_time) + " ps"
        coor_N = u.selectAtoms(atom_selection_N).coordinates()
        coor_C = u.selectAtoms(atom_selection_C).coordinates()

        coor_N5 = u.selectAtoms(atom_selection_N5).coordinates()
        coor_N8 = u.selectAtoms(atom_selection_N8).coordinates()
        if len(coor_N) != len(coor_C) :
           print "Error: number of atoms in each group does not match"
           sys.exit(1)

        box = u.dimensions[0:3]
        if nframes == 1 :
           for i in xrange(0,len(coor_N)) :
               loc_N = coor_N[i]
               loc_C = coor_C[i]
               px = position(loc_C[0],loc_N[0])
               py = position(loc_C[1],loc_N[1])
               # apply periodic boundary conditions
               px = px - box[0]*anint(px/box[0])
               py = py - box[1]*anint(py/box[1])

               pos[i] = sqrt(px*px + py*py)
           pos=pos[np.argsort(pos)]
           
           jj = 0
           ll = 1
           Rbin[0] = pos[0]
           for i in xrange(0,len(coor_N)) :
               if ll <= Nlayer :
                  if jj+1 == Nl_atoms[ll-1] :
                     Rbin[ll] = (pos[i] + pos[i+1])/2
                     jj = 0
                     ll += 1
                  else :
                     jj += 1
               else :
                  break  
           print str(Rbin) 
        j_5CB = 0
        j_8CB = 0
        for i in xrange(0,len(coor_N)) :
            loc_N = coor_N[i]
            loc_C = coor_C[i]
            px = position(loc_C[0],loc_N[0])
            py = position(loc_C[1],loc_N[1])
            # apply periodic boundary conditions
            px = px - box[0]*anint(px/box[0])
            py = py - box[1]*anint(py/box[1])

            RR = sqrt(px*px + py*py)
            if RR < Rbin[Nlayer] : 
               for j in range(Nlayer,0,-1) : 
                   if RR < Rbin[j] and RR > Rbin[j-1] :
                      NN = j - 1
                      break
               Theta = atan2(py,px) 
               if Theta<0 : Theta += 2*pi
               Theta = Theta*180/pi - angle[NN]
               if Theta<0 : Theta += 360
               
               indx = Sbin[NN] + int(Theta/dtheta[NN])
               
               dens[indx] += 1
               if (j_5CB<natoms5 and np.allclose(loc_N, coor_N5[j_5CB])) :
                  j_5CB += 1
                  dens5[indx] += 1
               if (j_8CB<natoms8 and np.allclose(loc_N, coor_N8[j_8CB])) :
                  j_8CB += 1
                  dens8[indx] += 1
               
               vv = vector(loc_C,loc_N)
               vv = vv/np.linalg.norm(vv)
               Q = (3*np.dot(vv.T,vv) - I)/2.0 # factor 3/2 means that biggest eigvalu = nematic order
               Q = np.reshape(Q,(9))
               Q_slab[indx,:] = Q + Q_slab[indx,:]

            else :
               if (j_5CB<natoms5 and np.allclose(loc_N, coor_N5[j_5CB])) :  j_5CB += 1
               if (j_8CB<natoms8 and np.allclose(loc_N, coor_N8[j_8CB])) :  j_8CB += 1

print "Total number of frames " + str(nframes)

for i in range(0,Nlayer) :
    VV = (pi*(Rbin[i+1]*Rbin[i+1]-Rbin[i]*Rbin[i]) * box[2]) / 1000
    for j in range(0,Nbin[i]):
        indx = Sbin[i] + j
        volume[indx] = VV / Nbin[i]
        Q = Q_slab[indx,:] / dens[indx]
        Q = np.reshape(Q,(3,3))
        e_value, e_vector = np.linalg.eigh(Q) # Calculate local order parameter
        e_value = e_value[np.argsort(e_value)]
        Scalar[indx] = e_value[2]
        Biaxial[indx] = 2*e_value[1] + e_value[2]
        vector_x[indx] = e_vector[0,np.argmax(e_value)]
        vector_y[indx] = e_vector[1,np.argmax(e_value)]
        vector_z[indx] = e_vector[2,np.argmax(e_value)]

### write coordinates    
f2=open('in.poly', 'w+')
f2.write('%5i 2 0 0 \n' % (Tbin))
zz = 0.0
i=0
output_file.write('# vtk DataFile Version 2.0 \n')
output_file.write('triangle \n')
output_file.write('ASCII \n')
output_file.write('DATASET POLYDATA \n')
output_file.write('POINTS %5i float \n'% (Tbin))
output_file.write('\n')

for i in xrange(0,Nlayer):
    Raverage = (Rbin[i]+Rbin[i+1])/2
    for j in xrange(0,Nbin[i]):
        theta = (j*dtheta[i]+angle[i])*pi/180 
        xx = Raverage*cos(theta)
        yy = Raverage*sin(theta)
        output_file.write('%5.3f %5.3f %5.3f \n' % ( xx, yy, zz))
        indx = Sbin[i] + j
        f2.write('%5i %5.3f %5.3f \n' % (indx+1, xx, yy))

f2.write('%5i 0 \n' % (Nbin[0]))
for i in range(1,Nbin[0]):
    f2.write('%5i %5i %5i \n' % (i, i, i+1))
f2.write('%5i %5i %5i \n' % (Nbin[0], 1, Nbin[0]))
f2.write('%5i 0 \n' % (1))
f2.write('%5i %5.3f %5.3f \n' % (1, 0.0, 0.))
f2.close()
#
##check_output("rm in.1.*", shell=True)
os.system("rm in.1.*")
##os.system("/home/amin/Downloads/triangle/triangle -pcBev in.poly")
check_output("/home/amin/Software/triangle/triangle -pcBev in.poly", shell=True)
##check_output("/home/amin/Downloads/triangle/showme in.poly", shell=True)
#    
#
## makes triangle element
f3=open('in.1.ele', 'r')
line = f3.readline()
columns = line.split()
Nelement = int(columns[0])
output_file.write('\n')
output_file.write('TRIANGLE_STRIPS %5i %5i \n'% (Nelement,4*Nelement))
output_file.write('\n')

for i in range(0,Nelement):
    line = f3.readline()
    columns = line.split()
    output_file.write('%1i %5i %5i %5i \n'% (3,int(columns[1])-1,int(columns[2])-1,int(columns[3])-1))


output_file.write('\n')
output_file.write('POINT_DATA %5i \n'% (Tbin))
output_file.write('SCALARS dens float \n')
output_file.write('LOOKUP_TABLE default \n')
output_file.write('\n')

for i in xrange(0,Tbin):
    output_file.write('%5.3f \n' % (float(dens[i])/volume[i]/float(nframes)))

output_file.write('\n')
output_file.write('SCALARS dens5 float \n')
output_file.write('LOOKUP_TABLE default \n')
output_file.write('\n')

for i in xrange(0,Tbin):
    output_file.write('%5.3f \n' % (float(dens5[i])/volume[i]/float(nframes)))

output_file.write('\n')
output_file.write('SCALARS dens8 float \n')
output_file.write('LOOKUP_TABLE default \n')
output_file.write('\n')

for i in xrange(0,Tbin):
    output_file.write('%5.3f \n' % (float(dens8[i])/volume[i]/float(nframes)))

output_file.write('\n')
output_file.write('SCALARS dens58 float \n')
output_file.write('LOOKUP_TABLE default \n')
output_file.write('\n')

for i in xrange(0,Tbin):
    output_file.write('%5.3f \n' % (float(dens8[i]-dens5[i])/float(dens[i])))

output_file.write('\n')
output_file.write('SCALARS scalar float \n')
output_file.write('LOOKUP_TABLE default \n')
output_file.write('\n')

for i in xrange(0,Tbin):
    output_file.write('%5.3f \n' % (Scalar[i]))
 
output_file.write('\n')
output_file.write('SCALARS biaxial float \n')
output_file.write('LOOKUP_TABLE default \n')
output_file.write('\n')

for i in xrange(0,Tbin):
    output_file.write('%5.3f \n' % (Biaxial[i]))
   
output_file.write('\n')
output_file.write('VECTORS directors float \n')
output_file.write('\n')

for i in xrange(0,Tbin):
    output_file.write('%5.3f %5.3f %5.3f \n' % (vector_x[i],vector_y[i],vector_z[i]))
#for ll in xrange(0,Nlayer):
#    for bb in xrange(0,Nbin[ll]):
#        NN = Sbin[ll]+bb
##        if (ll%10)==0 and ll>0 :
#        output_file.write('%5.3f %5.3f %5.3f \n' % (vector_x[NN],vector_y[NN],vector_z[NN]))
##        else :
##           output_file.write('%5.3f %5.3f %5.3f \n' % (0.0,0.0,0.0))
#
#


