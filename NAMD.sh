#!/bin/bash
#Automaticaly submit jobs to continue the simulation

if [ $1 == "start" ]; then
   echo "Start"
   timeout 85500s mpirun -n 16 namd2 sample.conf > log
#   timeout 85500s mpirun -n 32 /project/depablo/mrahimi/NAMD/NAMD_src/Source/Linux-x86_64-g++/namd2 sample.conf > log
   cp sample.conf sample.conf.bck
   cp job_midway_depablo job_midway_depablo.bck
elif [ $1 == "restart" ]; then
   echo "Restart"
#modify sample file
   sed -i '19s/.*/if {1} {/' sample.conf
   sed -i '39s/.*/#temperature         $temperature/' sample.conf
   sed -i '45s/.*/ if {0} {/' sample.conf
#   sed -i '144s/.*/ if {0} {/' sample.conf

   if [ $2 == "step" ]; then
      echo "Continue by matching timestep"
      timestep=$(awk '{if(NR==3) print $1}' myinput.restart.xsc )
      sed -i '26s/.*/firsttimestep      '$timestep'/' sample.conf
      echo $timestep
   fi
   timeout 85500s mpirun -n 16 namd2 sample.conf > log
#   timeout 85500s mpirun -n 32 /project/depablo/mrahimi/NAMD/NAMD_src/Source/Linux-x86_64-g++/namd2 sample.conf > log
elif [ $1 == "clean" ]; then
  rm -r DCD/ LOG/ Restart/ myoutput.* 
   mv sample.conf.bck sample.conf
   mv job_midway_depablo.bck job_midway_depablo
  exit
else
   echo "Incorrect input"
   exit
fi

if ! [ -d DCD ]
then
  mkdir DCD
  mkdir LOG
  mkdir Restart 
fi

cd DCD 
newfile=0
for file in `ls *.dcd|sort -g`
do

    filename=$(basename "$file")
    name=${filename%.*}
    extension=${filename##*.}

    echo $filename
    echo $name
    echo $extension

    if [ $name -ge $newfile ]
    then
      newfile=$name
      #mv "$file" "$(($filename + 1))".$extension
    fi  
done

newfile=$((newfile+1))
echo $newfile
cd ..

#Back up file
file=myoutput
cp "$file".dcd  DCD/"$newfile".dcd
cp log          LOG/"$newfile".log
cp "$file".restart.coor Restart/"$newfile".coor
cp "$file".restart.vel Restart/"$newfile".vel
cp "$file".restart.xsc Restart/"$newfile".xsc

#Copy output to input for restart
cp myoutput.restart.coor myinput.restart.coor
cp myoutput.restart.vel myinput.restart.vel
cp myoutput.restart.xsc myinput.restart.xsc

timestep=$(awk '{if(NR==3) print $1}' myinput.restart.xsc )
echo "Last timestep" $timestep
Ttimestep=$(awk '{if($1=="run") print $2}' sample.conf )
echo "Target timestep" $Ttimestep

if [ $timestep -lt $Ttimestep ]; then
   echo "Submit new job"
   awk '{if(NR==18) print $1, "restart step"; else print $0}' job_midway_depablo > ttt
   mv ttt job_midway_depablo
   sbatch job_midway_depablo
fi
sleep 60
echo "END"

