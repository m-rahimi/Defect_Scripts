#! /bin/bash

lcDir="/home/ramezani/midway/lc"

##
####################################################################################################
##5CB_vac @ 16c
#first=10
#last=80
##Navigate to the directory
#cd $lcDir/5cb/vac/16c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -T 289.15 #-b 1 -e 5
##
####################################################################################################
##5CB_vac @ 20c
#first=10
#last=80
##Navigate to the directory
#cd $lcDir/5cb/vac/20c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -T 293.15 #-b 1 -e 5
##
####################################################################################################
##5CB_vac @ 25c
#first=10
#last=80
##Navigate to the directory
#cd $lcDir/5cb/vac/25c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -T 298.15 #-b 1 -e 5
##
####################################################################################################
##5CB_vac @ 38c
#first=10
#last=80
##Navigate to the directory
#cd $lcDir/5cb/vac/38c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w5cb_vac.psf -f $first -l $last -n 100 -T 311.15 #-b 1 -e 5
##
####################################################################################################
##8CB_vac @ 38c
#first=10
#last=80
##Navigate to the directory
#cd $lcDir/8cb/vac/38c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -T 310.85 #-b 1 -e 5
##
####################################################################################################
##8CB_vac @ 43c
#first=10
#last=80
##Navigate to the directory
#cd $lcDir/8cb/vac/43c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_vac.psf -f $first -l $last -n 100 -T 316.15 #-b 1 -e 5
##
###################################################################################################
#5CB_wat @ 25c
first=60
last=78
#Navigate to the directory
cd $lcDir/5cb/wat/25c/DCD
#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
$lcDir/anchoring.py  -p ../w5cb_wat.psf -f $first -l $last -n 100 -T 298.15 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" #-b 1 -e 5
#
###################################################################################################
#5CB_wat @ 40c
first=60
last=97
#Navigate to the directory
cd $lcDir/5cb/wat/40c/DCD
#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
$lcDir/anchoring.py  -p ../w5cb_wat.psf -f $first -l $last -n 100 -T 313.15 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" #-b 1 -e 5
#
##################################################################################################
#5CB_wat_nacl @ 25c
first=60
last=71
#Navigate to the directory
cd $lcDir/5cb/nacl/2m/25c/DCD
#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
$lcDir/anchoring.py  -p ../w5cb_wat_nacl.psf -f $first -l $last -n 100 -T 298.15 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" #-b 1 -e 5
#
###################################################################################################
#5CB_wat_nacl @ 40c
first=60
last=104
#Navigate to the directory
cd $lcDir/5cb/nacl/2m/40c/DCD
#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
$lcDir/anchoring.py  -p ../w5cb_wat_nacl.psf -f $first -l $last -n 100 -T 313.15 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" #-b 1 -e 5
#
###################################################################################################
#5CB_wat_nai @ 25c
first=60
last=78
#Navigate to the directory
cd $lcDir/5cb/nai/2m/25c/DCD
#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
$lcDir/anchoring.py  -p ../w5cb_wat_nai.psf -f $first -l $last -n 100 -T 298.15 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" #-b 1 -e 5
#
###################################################################################################
#5CB_wat_nai @ 40c
first=60
last=101
#Navigate to the directory
cd $lcDir/5cb/nai/2m/40c/DCD
#Excecute the python script to calculate anchoring strength and distribution of cos(beta)
$lcDir/anchoring.py  -p ../w5cb_wat_nai.psf -f $first -l $last -n 100 -T 313.15 -sc "name CA12 and prop z < 320" -sn "name NY1 and prop z < 320" #-b 1 -e 5
##
####################################################################################################
##8CB_wat @ 25c
#first=10
#last=38
##Navigate to the directory
#cd $lcDir/8cb/wat/25c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -T 298.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat @ 37c
#first=10
#last=37
##Navigate to the directory
#cd $lcDir/8cb/wat/37c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -T 310.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat @ 47c
#first=10
#last=74
##Navigate to the directory
#cd $lcDir/8cb/wat/47c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat.psf -f $first -l $last -n 100 -T 320.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat_nacl @ 25c
#first=10
#last=46
##Navigate to the directory
#cd $lcDir/8cb/nacl/2m/25c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -T 298.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat_nacl @ 37c
#first=10
#last=46
##Navigate to the directory
#cd $lcDir/8cb/nacl/2m/37c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -T 310.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat_nacl @ 47c
#first=10
#last=82
##Navigate to the directory
#cd $lcDir/8cb/nacl/2m/47c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat_nacl.psf -f $first -l $last -n 100 -T 320.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat_nai @ 25c
#first=10
#last=41
##Navigate to the directory
#cd $lcDir/8cb/nai/2m/25c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -T 298.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat_nai @ 37c
#first=10
#last=41
##Navigate to the directory
#cd $lcDir/8cb/nai/2m/37c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -T 310.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##
####################################################################################################
##8CB_wat_nai @ 47c
#first=10
#last=76
##Navigate to the directory
#cd $lcDir/8cb/nai/2m/47c/DCD
##Excecute the python script to calculate anchoring strength and distribution of cos(beta)
#$lcDir/anchoring.py  -p ../w8cb_wat_nai.psf -f $first -l $last -n 100 -T 320.15 -sc "name CA12 and prop z < 450" -sn "name NY1 and prop z < 450" #-b 1 -e 5
##