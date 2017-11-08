#!/usr/bin/python

import sys
import subprocess
import time
import numpy as np
import os

ncpu=int(sys.argv[1])

proc=list(range(ncpu))

dist_kappa_dyn =['step30','step45']
kappa_min_dyn  =[0.5,1.0]

dist_kappa_w   =['step30','step45']
kappa_min_w    =[0.1,0.5]
kappa_max_w    =[1.0,5.0]

kappa_min_vis  =[1.0,5.0,10.,30.]

input_str=[]
input_folder=[]
counter=0

for i in range(len(dist_kappa_dyn)):
  for j in range(len(kappa_min_dyn)):
    for k in range(len(dist_kappa_w)):
      for l in range(len(kappa_min_w)):
        for m in range(len(kappa_min_vis)):

          input_par="'"+dist_kappa_dyn[i]+"'"+' '+str(kappa_min_dyn[j])+' '+"'"+dist_kappa_w[k]+"'"+' '+str(kappa_min_w[l])+' '+str(kappa_max_w[l])+' '+str(kappa_min_vis[m])
          input_folder.append(dist_kappa_dyn[i]+'_'+str(kappa_min_dyn[j])+'_'+dist_kappa_w[k]+'_'+str(kappa_min_w[l])+'_'+str(kappa_max_w[l])+'_'+str(kappa_min_vis[m]))

          input_str.append(input_par)
          print(counter,input_str[counter])
          counter=counter+1

icount = 0
for icpu in range(ncpu):

# create folder
  folder_name='lc_'+str(input_folder[icount])
  os.makedirs(folder_name)
  os.chdir(folder_name)
  os.makedirs("./output")
  os.system('cp -r ../source/lightcurve_multicomp_multiang_v4.py ./')

#  cmd = "python /scratch/snx1600/albino/mass_disk_david/scripts/mass_disc_serial.py "+str(icpu)
  cmd = "python ./lightcurve_multicomp_multiang_v4.py "+input_str[icount]
  proc[icpu]=subprocess.Popen(cmd,shell=True)
  icount = icount + 1
  os.chdir("..")
  time.sleep(30)

working = list(range(ncpu))
while (working):
  working=False
  for icpu in range(ncpu):
    if ( proc[icpu].poll() == None ): working = True
  time.sleep(30)

print('done')

