#!/usr/bin/env python2
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
import numpy as np 
import math
from lib import scripts 
from lib import analysis

def split_list(a_list):
    half = len(a_list)//2
    return a_list[:half], a_list[half:]

# Initiate some arrays

ligand_list = []
lambdas = []  
weights = []  
mols = []

# Defaults

e_steps1 = 0
e_steps2 = 0
x_steps1 = 0
x_steps2 = 0
ti_points = 0

add_sc = '1'

# Read arguments that define input file and stage
if len(sys.argv) < 5:
  scripts.help_message()
  sys.exit(0)
for i in [1, 3]:
  if '-i' == sys.argv[i].lower():
    input_file = sys.argv[i + 1]
  elif '-s' == sys.argv[i].lower():
    stage = sys.argv[i + 1]
  else:
    scripts.help_message()
    sys.exit(1)

# Open input file
with open(input_file) as f_in:       
    # Remove spaces and tabs
    lines = (line.strip(' \t\n\r') for line in f_in)
    lines = list(line for line in lines if line)  # Non-blank lines in a list

for i in range(0, len(lines)):
    # split line using the equal sign, and remove text after #
    if not lines[i][0] == '#':
        lines[i] = lines[i].split('#')[0].split('=')

# Read parameters from input file 
for i in range(0, len(lines)):
    if not lines[i][0] == '#':
        lines[i][0] = lines[i][0].strip().lower()
        lines[i][1] = lines[i][1].strip()

        if lines[i][0] == 'temperature':
            temperature = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'e_steps1':
            e_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'e_steps2':
            e_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'x_steps1':
            x_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'x_steps2':
            x_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'ti_points':
            ti_points = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'gti_add_sc':
            add_sc = lines[i][1]
        elif lines[i][0] == 'ligand_list':
            newline = lines[i][1].strip('\'\"-,.:;#()][').split(',')
            for j in range(0, len(newline)):
                ligand_list.append(newline[j])
        elif lines[i][0] == 'dec_int':
            if lines[i][1].lower() == 'mbar':
                dec_int = lines[i][1].lower()
            elif lines[i][1].lower() == 'ti':
                dec_int = lines[i][1].lower()
            else:
                print('Integration method not recognized, please choose ti or mbar')
                sys.exit(1)
        elif lines[i][0] == 'blocks':
            blocks = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'hmr':
            if lines[i][1].lower() == 'yes':
                hmr = 'yes'
            elif lines[i][1].lower() == 'no':
                hmr = 'no'
            else:
                print('Wrong input! Please use yes or no to indicate whether hydrogen mass repartitioning '
                      'will be used.')
                sys.exit(1)
        elif lines[i][0] == 'lambdas':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                lambdas.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))


# Obtain all ligand names
mols = []
for i in range(0, len(ligand_list)):
  with open('./all-poses/%s.pdb' %ligand_list[i].lower()) as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   
    for j in range(0, len(lines)):
      if (lines[j][0:6].strip() == 'ATOM') or (lines[j][0:6].strip() == 'HETATM'):  
        lig_name = (lines[j][17:20].strip())
        mols.append(lig_name)
        break 

# Obtain Gaussian Quadrature lambdas and weights

if dec_int == 'ti': 
  if ti_points != 0:
    lambdas = []
    weights = []
    x,y = np.polynomial.legendre.leggauss(ti_points)
    # Adjust Gaussian lambdas
    for i in range(0, len(x)):
      lambdas.append(float((x[i]+1)/2))
    # Adjust Gaussian weights
    for i in range(0, len(y)):
      weights.append(float(y[i]/2))
  else:
    print('Wrong input! Please choose a positive integer for the ti_points variable when using the TI-GQ method')
    sys.exit(1)
  print('lambda values:', lambdas) 
  print('Gaussian weights:', weights) 
elif dec_int == 'mbar': 
  if lambdas == []:
    print('Wrong input! Please choose a set of lambda values when using the MBAR method')
    sys.exit(1)
  if ti_points != 0:
    print('Wrong input! Do not define the ti_points variable when applying the MBAR method, instead choose a set of lambda values')
    sys.exit(1)
  print('lambda values:', lambdas) 




if stage == 'rbfe':

  # Call the original BAT program to create initial systems

  sp.call('python BAT.py -i %s -s fe' % input_file, shell=True)

  # Find the common and soft-core atoms for the two ligands

  print('Ligand names:')
  print(ligand_list)
  print('Ligand codes:')
  print(mols)

  cc_atoms = []
  sc1_atoms = []
  sc2_atoms = []
  atoms_lig = []
  ref_atoms = []
  target_atoms = []


  # Open reference ligand file and save atoms:
  ligand = ligand_list[0]
  with open('fe/%s/sdr/x00/vac_ligand.pdb' %ligand)  as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

  for i in range(0, len(lines)):
    if (lines[i][17:20].strip() == mols[0]):
      ref_atoms.append(lines[i][12:16].strip())

  # Open target ligand file and save atoms:
  ligand = ligand_list[1]
  with open('fe/%s/sdr/x00/vac_ligand.pdb' %ligand)  as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

  for i in range(0, len(lines)):
    if (lines[i][17:20].strip() == mols[1]):
      target_atoms.append(lines[i][12:16].strip())

  # Find common atoms
  for i in ref_atoms:
    if i in target_atoms:
      cc_atoms.append(i)

  for i in ref_atoms:
    if i not in cc_atoms:
      sc2_atoms.append(i)

  for i in target_atoms:
    if i not in cc_atoms:
      sc1_atoms.append(i)

  print('Reference ligand atoms:') 
  print(ref_atoms) 
  print('Total: %s' %(str(len(ref_atoms)))) 
  print('Target ligand atoms:') 
  print(target_atoms)
  print('Total: %s' %(str(len(target_atoms)))) 
  print('Common atoms:') 
  print(cc_atoms)
  print('Total: %s' %(str(len(cc_atoms)))) 
  print('Unique (soft-core) reference ligand atoms:') 
  print(sc2_atoms)
  print('Total: %s' %(str(len(sc2_atoms)))) 
  print('Unique (soft-core) target ligand atoms:') 
  print(sc1_atoms)
  print('Total: %s' %(str(len(sc1_atoms)))) 
  print('') 

  # Build new system for LJ calculation
  
  if not os.path.exists('rbfe'):
    os.makedirs('rbfe')
  os.chdir('rbfe')
  os.makedirs('%s_%s' % (ligand_list[0],ligand_list[1]))
  os.chdir('%s_%s' % (ligand_list[0],ligand_list[1]))
  os.makedirs('LJ-sc-exchange')
  os.chdir('LJ-sc-exchange')
  os.makedirs('build_files')
  os.chdir('build_files')
  shutil.copy('../../../../fe/%s/sdr/x00/tleap_solvate.in' %ligand_list[1], './')
  shutil.copy('../../../../fe/%s/sdr/x00/tleap_vac.in' %ligand_list[1], './')
  shutil.copy('../../../../fe/%s/sdr/x00/%s.mol2' %(ligand_list[1],mols[0].lower()), './')
  shutil.copy('../../../../fe/%s/sdr/x00/%s.mol2' %(ligand_list[1],mols[1].lower()), './')
  shutil.copy('../../../../fe/%s/sdr/x00/%s.frcmod' %(ligand_list[1],mols[0].lower()), './')
  shutil.copy('../../../../fe/%s/sdr/x00/%s.frcmod' %(ligand_list[1],mols[1].lower()), './')
  shutil.copy('../../../../fe/%s/sdr/x00/dum.mol2' %ligand_list[1], './')
  shutil.copy('../../../../fe/%s/sdr/x00/dum.frcmod' %ligand_list[1], './')
  shutil.copy('../../../../fe/%s/sdr/x00/build.pdb' %ligand_list[1], './build-ini.pdb')
  shutil.copy('../../../../fe/%s/sdr/x00/vac.pdb' %ligand_list[1], './vac-ini.pdb')
  shutil.copy('../../../../fe/%s/sdr/x00/cv.in' %ligand_list[1], './')
  shutil.copy('../../../../fe/%s/sdr/x00/run-local.bash' %ligand_list[1], './')
  shutil.copy('../../../../fe/%s/sdr/x00/disang.rest' %ligand_list[1], './')
  shutil.copy('../../../../fe/%s/sdr/x00/parmed-hmr.in' %ligand_list[1], './')

  # Clean restraint file
  file = open("disang.rest", "r")
  lines=file.readlines()
  file.close()
  f = open("disang.rest", "w")
  j = 0
  for line in lines:
    f.writelines(line)
    j += 1
    if j == 4:
      f.writelines('\n')
      break
  f.close()

  # Read original build file

  with open('build-ini.pdb') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

  with open('./build.pdb', 'w') as outfile:
    for i in range(0,len(lines)):
      if lines[i][17:20].strip() == '%s' %mols[1]:
        break
      outfile.write(lines[i]+'\n')
  outfile.close()

  lines_after = []

  k = 0
  for i in range(0,len(lines)):
    if lines[i][17:20].strip() == '%s' %mols[1]:
      k = 1
    if k == 1 and  lines[i][17:20].strip() != '%s' %mols[0] and  lines[i][17:20].strip() != '%s' %mols[1]:
      lines_after.append(lines[i]+'\n')

  with open('vac-ini.pdb') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

  lines_ligand = []
  
  for i in range(0,len(lines)):
    if lines[i][17:20].strip() == '%s' %mols[1]:
      lines_ligand.append(lines[i]+'\n')

  lines1, lines2 = split_list(lines_ligand)

  # Write final build file

  build_file = open('build.pdb', 'a')
  for i in range(0,len(lines1)):
    if lines1[i][12:16].strip() in cc_atoms:
      build_file.write(lines1[i])        
  build_file.write('TER\n')       
  for i in range(0,len(lines2)):
    if lines2[i][12:16].strip() in cc_atoms:
      build_file.write(lines2[i].replace(mols[1],mols[0]))        
  build_file.write('TER\n')       
  for i in range(0,len(lines1)):
    if lines1[i][12:16].strip() in cc_atoms:
      build_file.write(lines1[i].replace(mols[1],mols[0]))        
  build_file.write('TER\n')       
  for i in range(0,len(lines2)):
    if lines2[i][12:16].strip() in cc_atoms:
      build_file.write(lines2[i])        
  for i in range(0,len(lines_after)):
    if i >= 3:
      build_file.write(lines_after[i])        
  build_file.close()       

  # Write dry build file

  with open('build.pdb') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   
  with open('./build-dry.pdb', 'w') as outfile:
    for i in range(0,len(lines)):
      if lines[i][17:20].strip() == 'WAT':
        break
      outfile.write(lines[i]+'\n')

  outfile.close()

  # Solvate system with existing parameters
  
  p = sp.call('tleap -s -f tleap_vac.in > tleap_vac.log', shell=True)

  p = sp.call('tleap -s -f tleap_solvate.in > tleap_solvate.log', shell=True)

  p = sp.call('parmed -O -n -i parmed-hmr.in > parmed-hmr.log', shell=True)

  # Configure simulation files 

  # First ligand residue number
  with open('vac.pdb') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

  for i in range(0,len(lines)):
    if lines[i][17:20].strip() == '%s' %mols[1]:
      lig_resid = lines[i][22:26].strip()
      break

  # Convert atom arrays to variables
  sc1_atoms_list = ",".join(sc1_atoms)
  #print(sc1_atoms_list)
  sc2_atoms_list = ",".join(sc2_atoms)
  #print(sc2_atoms_list)

  os.chdir('../')

  # Create windows for LJ calculations

  print('LJ soft-core exchange windows:')
  for k in range(0, len(lambdas)):
    weight = lambdas[k]
    win = k
    print('window: x%02d lambda: %s' %(int(win), str(weight)))
    shutil.copytree('./build_files', './x%02d' % int(win))
    os.chdir('./x%02d' % int(win))
    shutil.copy('../../../../fe/%s/sdr/x%02d/mini.in' %(ligand_list[1],int(win)), './')
    shutil.copy('../../../../fe/%s/sdr/x%02d/heat.in' %(ligand_list[1],int(win)), './')
    shutil.copy('../../../../fe/%s/sdr/x%02d/eqnpt.in' %(ligand_list[1],int(win)), './')
    shutil.copy('../../../../fe/%s/sdr/x%02d/mdin-00' %(ligand_list[1],int(win)), './')
    shutil.copy('../../../../fe/%s/sdr/x%02d/SLURMM-run' %(ligand_list[1],int(win)), './')
    # Edit simulation files
    file = open("mini.in", "r")
    lines=file.readlines()
    file.close()
    f = open("mini.in", "w")
    newlines = []
    for line in lines:
        if "ifsc" in line:
            line = '  ifsc=1, scmask1=\':'+lig_resid+'@'+sc1_atoms_list+' | :'+str(int(lig_resid)+1)+'@'+sc2_atoms_list+'\', scmask2=\':'+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\', crgmask=\':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+lig_resid+','+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\'\n' 
        if "gti_chg_keep" in line:
            line = '  gti_chg_keep = 1,\n'
        if "gti_add_sc" in line:
            line = '  gti_add_sc = '+add_sc+',\n'
        newlines.append(line)
    f.writelines(newlines)
    f.close()
    file = open("heat.in", "r")
    lines=file.readlines()
    file.close()
    f = open("heat.in", "w")
    newlines = []
    for line in lines:
      if "ifsc" in line:
          line = '   ifsc=1, scmask1=\':'+lig_resid+'@'+sc1_atoms_list+' | :'+str(int(lig_resid)+1)+'@'+sc2_atoms_list+'\', scmask2=\':'+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\', crgmask=\':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+lig_resid+','+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\'\n' 
      if "gti_chg_keep" in line:
          line = '   gti_chg_keep = 1,\n'
      if "gti_add_sc" in line:
          line = '   gti_add_sc = '+add_sc+',\n'
      newlines.append(line)
    f.writelines(newlines)
    f.close()
    file = open("eqnpt.in", "r")
    lines=file.readlines()
    file.close()
    f = open("eqnpt.in", "w")
    newlines = []
    for line in lines:
      if "ifsc" in line:
          line = '  ifsc=1, scmask1=\':'+lig_resid+'@'+sc1_atoms_list+' | :'+str(int(lig_resid)+1)+'@'+sc2_atoms_list+'\', scmask2=\':'+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\', crgmask=\':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+lig_resid+','+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\'\n' 
      if "gti_chg_keep" in line:
          line = '  gti_chg_keep = 1,\n'
      if "gti_add_sc" in line:
          line = '  gti_add_sc = '+add_sc+',\n'
      newlines.append(line)
    f.writelines(newlines)
    f.close()
    for i in range(0,3):
      file = open("mdin-00", "r")
      lines=file.readlines()
      file.close()
      f = open("mdin-%02d" %i, "w")
      newlines = []
      for line in lines:
        if "ifsc" in line:
            line = '  ifsc=1, scmask1=\':'+lig_resid+'@'+sc1_atoms_list+' | :'+str(int(lig_resid)+1)+'@'+sc2_atoms_list+'\', scmask2=\':'+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\', crgmask=\':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+2)+'@'+sc2_atoms_list+' | :'+lig_resid+','+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\'\n' 
        if "gti_chg_keep" in line:
            line = '  gti_chg_keep = 1,\n'
        if "gti_add_sc" in line:
            line = '  gti_add_sc = '+add_sc+',\n'
        if "nstlim" in line:
          if i != 2:
            line = '  nstlim = %s\n' %(str(round(x_steps1/2)))
          else:
            line = '  nstlim = %s\n' %(str(x_steps2))
        newlines.append(line)
      f.writelines(newlines)
      f.close()
    os.chdir('../')

  os.chdir('../')

  # Create windows for charged calculations

  s = 0
  for i in ligand_list:
    print('Partial charging windows for %s:' %i)
    os.makedirs('./charge-%s' %i)
    os.chdir('./charge-%s' %i)
    for k in range(0, len(lambdas)):
      weight = lambdas[k]
      win = k
      print('window: e%02d lambda: %s' %(int(win), str(weight)))
      os.makedirs('./e%02d' % int(win))
      os.chdir('./e%02d' % int(win))
      for file in glob.glob('../../../../fe/%s/sdr/e%02d/vac.*' %(i,int(win))):
        shutil.copy(file, './')
      for file in glob.glob('../../../../fe/%s/sdr/e%02d/full.*' %(i,int(win))):
        shutil.copy(file, './')
      shutil.copy('../../../../fe/%s/sdr/e%02d/disang.rest' %(i,int(win)), './')
      shutil.copy('../../../../fe/%s/sdr/e%02d/cv.in' %(i,int(win)), './')
#      shutil.copy('../../../../fe/%s/sdr/e%02d/mini.in' %(i,int(win)), './')
      shutil.copy('../../../../fe/%s/sdr/e%02d/heat.in' %(i,int(win)), './')
      shutil.copy('../../../../fe/%s/sdr/e%02d/eqnpt.in' %(i,int(win)), './')
      shutil.copy('../../../../fe/%s/sdr/e%02d/mdin-00' %(i,int(win)), './')
      shutil.copy('../../../../fe/%s/sdr/e%02d/SLURMM-run' %(i,int(win)), './')
      shutil.copy('../../../../fe/%s/sdr/e%02d/run-local.bash' %(i,int(win)), './')

      # Clean restraint file
      file = open("disang.rest", "r")
      lines=file.readlines()
      file.close()
      f = open("disang.rest", "w")
      j = 0
      for line in lines:
        f.writelines(line)
        j += 1
        if j == 4:
          f.writelines('\n')
          break
      f.close()

      # Edit simulation files
      file = open("heat.in", "r")
      lines=file.readlines()
      file.close()
      f = open("heat.in", "w")
      newlines = []
      for line in lines:
        if "ifsc" in line:
          if s == 0:
            line = '   ifsc=0, crgmask = \':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+3)+'@'+sc2_atoms_list+'\',\n' 
          elif s == 1:
            line = '   ifsc=0, crgmask = \':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\',\n' 
        newlines.append(line)
      f.writelines(newlines)
      f.close()
      file = open("eqnpt.in", "r")
      lines=file.readlines()
      file.close()
      f = open("eqnpt.in", "w")
      newlines = []
      for line in lines:
        if "ifsc" in line:
          if s == 0:
            line = '  ifsc=0, crgmask = \':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+3)+'@'+sc2_atoms_list+'\',\n' 
          elif s == 1:
            line = '  ifsc=0, crgmask = \':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\',\n' 
        newlines.append(line)
      f.writelines(newlines)
      f.close()
      for b in range(0,3):
        file = open("mdin-00", "r")
        lines=file.readlines()
        file.close()
        f = open("mdin-%02d" %b, "w")
        newlines = []
        for line in lines:
          if "ifsc" in line:
            if s == 0:
              line = '  ifsc=0, crgmask = \':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+3)+'@'+sc2_atoms_list+'\',\n' 
            elif s == 1:
              line = '  ifsc=0, crgmask = \':'+str(int(lig_resid)+1)+','+str(int(lig_resid)+3)+'@'+sc1_atoms_list+'\',\n' 
          if "nstlim" in line:
            if b != 2:
              line = '  nstlim = %s\n' %(str(round(e_steps1/2)))
            else:
              line = '  nstlim = %s\n' %(str(e_steps2))
          newlines.append(line)
        f.writelines(newlines)
        f.close()
      os.chdir('../') 
    s += 1
    os.chdir('../')
  os.chdir('../../')  
  shutil.rmtree('./fe')

elif stage == 'analysis':
  if dec_int == 'ti':
    # Get LJ exchange free energies
    comp = 'x'
    os.chdir('rbfe/%s_%s/LJ-sc-exchange' % (ligand_list[0],ligand_list[1]))
    # Get dvdl values from output file
    for j in range(0, len(lambdas)):
      data = []
      win = j
      os.chdir('%s%02d' %(comp, int(win)))
      dvdl = open('dvdl.dat', "w")
      with open("md-02.out", "r") as fin:
        s = 0
        n = 0
        for line in fin:
          if 'TI region  1' in line:             
            s = 1
          if 'DV/DL  = ' in line and s == 1:
            n = n+1
            splitdata = line.split()
            data.append(splitdata[2])
            dvdl.write('%5d   %9.4f\n' % (n, float(splitdata[2])))
            s = 0
          if 'A V E' in line:
            break
        dvdl.close()
      # Separate in blocks
      for k in range(0, blocks):
        fout = open('dvdl%02d.dat' % (k+1), "w")
        for t in range(k*int(round(len(data)//blocks)), (k+1)*int(round(len(data)//blocks))):
          fout.write('%5d   %9.4f\n' % (t+1, float(data[t])))
        fout.close()
      os.chdir('../')
    os.chdir('../../../')
    pose = "_".join(ligand_list)
    rest_file = 'dvdl.dat'
    mode = 'all'
    dec_method = 'LJ-sc-exchange'
    analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

    for k in range(0, blocks):
      rest_file = 'dvdl%02d.dat' % (k+1) 
      mode = 'b%02d' % (k+1) 
      analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)
  
    sys.stdout = sys.__stdout__

    # Get free energy value
    os.chdir('rbfe/%s_%s/LJ-sc-exchange' % (ligand_list[0],ligand_list[1]))
    with open('./data/ti-'+comp+'-all.dat', "r") as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line)
      data = lines[-1]
      splitdata = data.split()
      fe_x = float(splitdata[1])

    # Get error
    b_data = [] 
    for k in range(0, blocks):
      with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line)
        data = lines[-1]
        splitdata = data.split()
        b_data.append(float(splitdata[1]))
    sd_x = np.std(b_data)
    os.chdir('../../../')

    print(fe_x)
    print(sd_x)

    # Get electrostatic free energies
    fe_e = []
    sd_e = []
    for i in ligand_list:
      comp = 'e'
      os.chdir('rbfe/%s_%s/charge-%s' % (ligand_list[0],ligand_list[1],i))
      # Get dvdl values from output file
      for j in range(0, len(lambdas)):
        data = []
        win = j
        os.chdir('%s%02d' %(comp, int(win)))
        dvdl = open('dvdl.dat', "w")
        with open("md-02.out", "r") as fin:
          s = 0
          n = 0
          for line in fin:
            if 'TI region  1' in line:             
              s = 1
            if 'DV/DL  = ' in line and s == 1:
              n = n+1
              splitdata = line.split()
              data.append(splitdata[2])
              dvdl.write('%5d   %9.4f\n' % (n, float(splitdata[2])))
              s = 0
            if 'A V E' in line:
              break
          dvdl.close()
        # Separate in blocks
        for k in range(0, blocks):
          fout = open('dvdl%02d.dat' % (k+1), "w")
          for t in range(k*int(round(len(data)//blocks)), (k+1)*int(round(len(data)//blocks))):
            fout.write('%5d   %9.4f\n' % (t+1, float(data[t])))
          fout.close()
        os.chdir('../')
      os.chdir('../../../')
      pose = "_".join(ligand_list)
      rest_file = 'dvdl.dat'
      mode = 'all'
      prefix = "charge"
      dec_method = f"{prefix}-{i}"
      analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

      for k in range(0, blocks):
        rest_file = 'dvdl%02d.dat' % (k+1) 
        mode = 'b%02d' % (k+1) 
        analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

      sys.stdout = sys.__stdout__

      # Get free energy value
      os.chdir('rbfe/%s_%s/%s' % (ligand_list[0],ligand_list[1],dec_method))
      with open('./data/ti-'+comp+'-all.dat', "r") as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line)
        data = lines[-1]
        splitdata = data.split()
        fe_e.append(float(splitdata[1]))

      # Get error
      b_data = [] 
      for k in range(0, blocks):
        with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
          lines = (line.rstrip() for line in f_in)
          lines = list(line for line in lines if line)
          data = lines[-1]
          splitdata = data.split()
          b_data.append(float(splitdata[1]))
        stdev = np.std(b_data)
      sd_e.append(stdev)
      os.chdir('../../../')

    print(fe_e)
    print(sd_e)
  elif dec_int == 'mbar':
    # Get LJ exchange free energies
    comp = 'x'
    os.chdir('rbfe/%s_%s/LJ-sc-exchange' % (ligand_list[0],ligand_list[1]))
    # Get potential energy values from output file
    for j in range(0, len(lambdas)):
      data = []
      win = j
      os.chdir('%s%02d' %(comp, int(win)))
      potl = open('energies.dat', "w")
      with open("md-02.out", "r") as fin:
        n = 0
        for line in fin:
          cols = line.split()
          if 'MBAR Energy analysis' in line:
            if n != 0:
              potl.write('\n')
            n = n+1
          if len(cols) >= 2 and cols[0] == 'Energy' and cols[1] == 'at':
            potl.write('%5d  %6s   %10s\n' % (n, cols[2], cols[4]))
      potl.write('\n')
      potl.close()
      # Separate in blocks
      for k in range(0, blocks):
        s = 0
        fout = open('ener%02d.dat' % (k+1), "w")
        with open("energies.dat", "r") as fin:
          for line in fin:
            cols = line.split()
            low = int(k*int(round(n/blocks)))+1
            high = int((k+1)*int(round(n/blocks)))+1
            if len(cols) >= 1 and int(cols[0]) == low:
              s = 1
            if len(cols) >= 1 and int(cols[0]) == high:
              s = 0
            if s == 1:
              fout.write(line)
        fout.close()
      os.chdir('../')
    os.chdir('../../../')
    pose = "_".join(ligand_list)
    rest_file = 'energies.dat'
    mode = 'all'
    dec_method = 'LJ-sc-exchange'
    analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

    for k in range(0, blocks):
      rest_file = 'ener%02d.dat' % (k+1) 
      mode = 'b%02d' % (k+1) 
      analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)
  
    sys.stdout = sys.__stdout__

    # Get free energy value
    os.chdir('rbfe/%s_%s/LJ-sc-exchange' % (ligand_list[0],ligand_list[1]))
    with open('./data/mbar-'+comp+'-all.dat', "r") as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line)
      data = lines[-1]
      splitdata = data.split()
      fe_x = float(splitdata[1])

    # Get error
    b_data = [] 
    for k in range(0, blocks):
      with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line)
        data = lines[-1]
        splitdata = data.split()
        b_data.append(float(splitdata[1]))
    sd_x = np.std(b_data)
    os.chdir('../../../')

    print(fe_x)
    print(sd_x)

    # Get electrostatic free energies
    fe_e = []
    sd_e = []
    for i in ligand_list:
      comp = 'e'
      os.chdir('rbfe/%s_%s/charge-%s' % (ligand_list[0],ligand_list[1],i))
      # Get potential energy values from output file
      for j in range(0, len(lambdas)):
        data = []
        win = j
        os.chdir('%s%02d' %(comp, int(win)))
        potl = open('energies.dat', "w")
        with open("md-02.out", "r") as fin:
          n = 0
          for line in fin:
            cols = line.split()
            if 'MBAR Energy analysis' in line:
              if n != 0:
                potl.write('\n')
              n = n+1
            if len(cols) >= 2 and cols[0] == 'Energy' and cols[1] == 'at':
              potl.write('%5d  %6s   %10s\n' % (n, cols[2], cols[4]))
        potl.write('\n')
        potl.close()
        # Separate in blocks
        for k in range(0, blocks):
          s = 0
          fout = open('ener%02d.dat' % (k+1), "w")
          with open("energies.dat", "r") as fin:
            for line in fin:
              cols = line.split()
              low = int(k*int(round(n/blocks)))+1
              high = int((k+1)*int(round(n/blocks)))+1
              if len(cols) >= 1 and int(cols[0]) == low:
                s = 1
              if len(cols) >= 1 and int(cols[0]) == high:
                s = 0
              if s == 1:
                fout.write(line)
          fout.close()
        os.chdir('../')
      os.chdir('../../../')
      pose = "_".join(ligand_list)
      rest_file = 'energies.dat'
      mode = 'all'
      prefix = "charge"
      dec_method = f"{prefix}-{i}"
      analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

      for k in range(0, blocks):
        rest_file = 'ener%02d.dat' % (k+1) 
        mode = 'b%02d' % (k+1) 
        analysis.fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

      sys.stdout = sys.__stdout__

      # Get free energy value
      os.chdir('rbfe/%s_%s/%s' % (ligand_list[0],ligand_list[1],dec_method))
      with open('./data/mbar-'+comp+'-all.dat', "r") as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line)
        data = lines[-1]
        splitdata = data.split()
        fe_e.append(float(splitdata[1]))

      # Get error
      b_data = [] 
      for k in range(0, blocks):
        with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
          lines = (line.rstrip() for line in f_in)
          lines = list(line for line in lines if line)
          data = lines[-1]
          splitdata = data.split()
          b_data.append(float(splitdata[1]))
        stdev = np.std(b_data)
      sd_e.append(stdev)
      os.chdir('../../../')

    print(fe_e)
    print(sd_e)


  # Write final results
  os.chdir('rbfe/%s_%s/' % (ligand_list[0],ligand_list[1]))
  # Create Results folder
  if not os.path.exists('Results'):
    os.makedirs('Results')
  total_fe = fe_x + fe_e[1] - fe_e[0]
  total_sd2 = sd_x**2 + sd_e[0]**2 + sd_e[1]**2
  resfile = open('./Results/Results.dat', 'w')
  resfile.write('\n---------------------------------------------------\n')
  resfile.write('RBFE from molecules %s to %s' %(ligand_list[0],ligand_list[1]))
  resfile.write('\n---------------------------------------------------\n\n')
  resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
  resfile.write('%-20s %8.2f;    %3.2f\n' % ('Electro 2 ('+dec_int.upper()+');', fe_e[1], sd_e[1]))
  resfile.write('%-20s %8.2f;    %3.2f\n' % ('LJ exchange ('+dec_int.upper()+');', fe_x, sd_x))
  resfile.write('%-20s %8.2f;    %3.2f\n' % ('-Electro 1 ('+dec_int.upper()+');', -1*fe_e[0], sd_e[0]))
  resfile.write('\n%-20s %8.2f;    %3.2f\n' % ('Total RBFE value;', -1.0*total_fe, math.sqrt(total_sd2)))
  resfile.write('\n---------------------------------------------------\n\n')
  resfile.write('Energies in kcal/mol\n\n')
  resfile.close()



  # Write results for the blocks
  for k in range(0, blocks):
    fb_x = 0
    fb_e = []
    with open('./LJ-sc-exchange/data/'+dec_int+'-x-b%02d.dat' %(k+1), "r") as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line)
      data = lines[-1]
      splitdata = data.split()
      fb_x = float(splitdata[1])
    with open('./charge-'+ligand_list[0]+'/data/'+dec_int+'-e-b%02d.dat' %(k+1), "r") as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line)
      data = lines[-1]
      splitdata = data.split()
      fb_e.append(float(splitdata[1]))
    with open('./charge-'+ligand_list[1]+'/data/'+dec_int+'-e-b%02d.dat' %(k+1), "r") as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line)
      data = lines[-1]
      splitdata = data.split()
      fb_e.append(float(splitdata[1]))
    total_fb = fb_x + fb_e[1] - fb_e[0]
    resfile = open('./Results/Res-b%02d.dat' %(k+1), 'w')
    resfile.write('\n---------------------------------------------------\n')
    resfile.write('RBFE from molecules %s to %s' %(ligand_list[0],ligand_list[1]))
    resfile.write('\n---------------------------------------------------\n\n')
    resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
    resfile.write('%-20s %8.2f\n' % ('Electro 2 ('+dec_int.upper()+');', fb_e[1]))
    resfile.write('%-20s %8.2f\n' % ('LJ exchange ('+dec_int.upper()+');', fb_x))
    resfile.write('%-20s %8.2f\n' % ('-Electro 1 ('+dec_int.upper()+');', fb_e[0]))
    resfile.write('\n%-20s %8.2f\n' % ('Total RBFE value;', -1.0*total_fb))
    resfile.write('\n---------------------------------------------------\n\n')
    resfile.write('Energies in kcal/mol\n\n')
    resfile.close()

