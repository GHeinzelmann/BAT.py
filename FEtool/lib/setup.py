#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
import scripts as scripts

def restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp):

    rst = []
    atm_num = []
    mlines = []
    msk = []
    pdb_file = ('vac.pdb')
    ligand_pdb_file = ('vac_ligand.pdb')

    # Restraint identifiers
    recep_tr = '#Rec_TR'
    recep_c = '#Rec_C'
    recep_d = '#Rec_D'
    lign_tr = '#Lig_TR'
    lign_c = '#Lig_C'
    lign_d = '#Lig_D'

    # Change to simulation directory
    if stage != 'fe':
      os.chdir(pose)

    # Find anchors
    with open(stage+'-%s.pdb' % mol.lower(), 'r') as f:
      data = f.readline().split()    
      P1 = data[2].strip()   
      P2 = data[3].strip()   
      P3 = data[4].strip()   
      L1 = data[5].strip()   
      L2 = data[6].strip()   
      L3 = data[7].strip()   
      lig_res = L1.split('@')[0][1:]
      first_res = data[8].strip()        
      p1_res = P1.split('@')[0][1:]  
      p2_res = P2.split('@')[0][1:]    
      p3_res = P3.split('@')[0][1:]    
    
    # Adjust anchors for ligand only
    if (comp == 'c' or comp == 'w' or comp == 'f'):
      L1 = L1.replace(':'+lig_res, ':1') 
      L2 = L2.replace(':'+lig_res, ':1') 
      L3 = L3.replace(':'+lig_res, ':1') 

    # Adjust anchors for protein only
    if (comp == 'r'):
      p1_new = int(p1_res) - 3    
      p2_new = int(p2_res) - 3    
      p3_new = int(p3_res) - 3    
      P1 = P1.replace(':'+p1_res, ':'+str(p1_new)) 
      P2 = P2.replace(':'+p2_res, ':'+str(p2_new)) 
      P3 = P3.replace(':'+p3_res, ':'+str(p3_new)) 
    

    # Get a relation between atom number and masks
    atm_num = scripts.num_to_mask(pdb_file)
    ligand_atm_num = scripts.num_to_mask(ligand_pdb_file)

    # Define translational/rotational and anchor atom distance restraints on the protein

    rst.append(':2@Pb '+P1+'')
    rst.append(':1@Pb :2@Pb '+P1+'') 
    rst.append(':3@Pb :1@Pb :2@Pb '+P1+'') 
    rst.append(':2@Pb '+P1+' '+P2+'') 
    rst.append(':1@Pb :2@Pb '+P1+' '+P2+'') 
    rst.append(':2@Pb '+P1+' '+P2+' '+P3+'') 
    rst.append(''+P1+' '+P2+'') 
    rst.append(''+P2+' '+P3+'') 
    rst.append(''+P3+' '+P1+'') 

    # Define protein dihedral restraints in the given range

    beg = bb_start - int(first_res) + 4 
    end = bb_end - int(first_res) + 4 
    nd = 0
    for i in range(beg, end+1):
      j = i+1
      psi1 = ':'+str(i)+'@N' 
      psi2 = ':'+str(i)+'@CA' 
      psi3 = ':'+str(i)+'@C' 
      psi4 = ':'+str(j)+'@N' 
      psit = '%s %s %s %s' % (psi1, psi2, psi3, psi4)
      rst.append(psit)
      nd += 1  
      phi1 = ':'+str(i)+'@C' 
      phi2 = ':'+str(j)+'@N' 
      phi3 = ':'+str(j)+'@CA' 
      phi4 = ':'+str(j)+'@C' 
      phit = '%s %s %s %s' % (phi1, phi2, phi3, phi4)
      rst.append(phit)
      nd += 1  

    # Define translational/rotational and anchor atom distance restraints on the ligand

    rst.append(':1@Pb '+L1+'')
    rst.append(':2@Pb :1@Pb '+L1+'') 
    rst.append(':3@Pb :2@Pb :1@Pb '+L1+'') 
    rst.append(':1@Pb '+L1+' '+L2+'') 
    rst.append(':2@Pb :1@Pb '+L1+' '+L2+'') 
    rst.append(':1@Pb '+L1+' '+L2+' '+L3+'') 
    rst.append(''+L1+' '+L2+'') 
    rst.append(''+L2+' '+L3+'') 
    rst.append(''+L3+' '+L1+'') 

    # New restraints for ligand only
    if (comp == 'c' or comp == 'w' or comp == 'f'):
      rst = []
      rst.append(''+L1+' '+L2+'') 
      rst.append(''+L2+' '+L3+'') 
      rst.append(''+L3+' '+L1+'') 

    # Get ligand dihedral restraints from ligand parameter/pdb file

    spool = 0
    with open('./vac_ligand.prmtop') as fin:
	lines = (line.rstrip() for line in fin)
	lines = list(line for line in lines if line) # Non-blank lines in a list   
	for line in lines:
	  if 'FLAG DIHEDRALS_WITHOUT_HYDROGEN' in line:
	    spool=1
	  elif 'FLAG EXCLUDED_ATOMS_LIST' in line:
	    spool=0
	  if spool != 0 and (len(line.split()) > 3):
	    mlines.append(line)


    for i in range(0, len(mlines)):
      data = mlines[i].split()
      if int(data[3]) > 0:
	anum = []
	for j in range (0, len(data)): 
	  anum.append(abs(int(data[j])/3)+1)      
	msk.append('%s %s %s %s' %(ligand_atm_num[anum[0]], ligand_atm_num[anum[1]], ligand_atm_num[anum[2]], ligand_atm_num[anum[3]]))   
	   
    for i in range(0, len(mlines)):
      data = mlines[i].split()
      if len(data) > 7:
	if int(data[8]) > 0:
	  anum = []
	  for j in range (0, len(data)): 
	    anum.append(abs(int(data[j])/3)+1)      
	  msk.append('%s %s %s %s' %(ligand_atm_num[anum[5]], ligand_atm_num[anum[6]], ligand_atm_num[anum[7]], ligand_atm_num[anum[8]]))   
	   
    excl = msk[:]
    ind = 0
    mat = []
    for i in range(0, len(excl)):
       data = excl[i].split()
       for j in range(0, len(excl)):   
	 if j == i:
	   break 
	 data2 = excl[j].split()
	 if (data[1] == data2[1] and data[2] == data2[2]) or (data[1] == data2[2] and data[2] == data2[1]):
	   ind = 0
	   for k in range(0, len(mat)):
	     if mat[k] == j:
	       ind = 1
	   if ind == 0:
	     mat.append(j) 

    for i in range(0, len(mat)):
      msk[mat[i]]= ''

    if (comp != 'c' and comp != 'w' and comp != 'f'):
      msk = filter(None, msk) 
      msk = [m.replace(':1',':'+lig_res) for m in msk]

    for i in range(0, len(msk)):
      rst.append(msk[i])

    # New restraints for protein only
    if (comp == 'r'):
      rst = []
      rst.append(''+P1+' '+P2+'') 
      rst.append(''+P2+' '+P3+'') 
      rst.append(''+P3+' '+P1+'') 
      beg = bb_start - int(first_res) + 1 
      end = bb_end - int(first_res) + 1
      nd = 0
      for i in range(beg, end+1):
	j = i+1
	psi1 = ':'+str(i)+'@N' 
	psi2 = ':'+str(i)+'@CA' 
	psi3 = ':'+str(i)+'@C' 
	psi4 = ':'+str(j)+'@N' 
	psit = '%s %s %s %s' % (psi1, psi2, psi3, psi4)
	rst.append(psit)
	nd += 1  
	phi1 = ':'+str(i)+'@C' 
	phi2 = ':'+str(j)+'@N' 
	phi3 = ':'+str(j)+'@CA' 
	phi4 = ':'+str(j)+'@C' 
	phit = '%s %s %s %s' % (phi1, phi2, phi3, phi4)
	rst.append(phit)
	nd += 1  


    # Get initial restraint values for references

    assign_file = open('assign.in', 'w')
    assign_file.write('%s  %s  %s  %s  %s  %s  %s\n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3))
    assign_file.write('parm full.hmr.prmtop\n')
    assign_file.write('trajin full.inpcrd\n')
    for i in range(0, len(rst)):
	arr = rst[i].split()
	if len(arr) == 2:
	  assign_file.write('%s %s %s'%('distance r'+str(i), rst[i], 'noimage out assign.dat\n'))
	if len(arr) == 3:
	  assign_file.write('%s %s %s'%('angle r'+str(i), rst[i], 'out assign.dat\n'))
	if len(arr) == 4:
	  assign_file.write('%s %s %s'%('dihedral r'+str(i), rst[i], 'out assign.dat\n'))

    assign_file.close() 
    sp.call('cpptraj -i assign.in > assign.log', shell=True)

    # Assign reference values for restraints
    with open('./assign.dat') as fin:
	lines = (line.rstrip() for line in fin)
	lines = list(line for line in lines if line) # Non-blank lines in a list   
	vals = lines[1].split()
	vals.append(vals.pop(0))
	del vals[-1]

    # For preparation or free energy, get pull distance and equilibrium values for protein backbone restraints
    if (stage == 'prep'):
      vals[9+nd] = float(vals[9+nd]) + trans_dist
      shutil.copy('../../equil/'+pose+'/assign.dat', './assign-eq.dat')
      with open('./assign-eq.dat') as fin:
	lines = (line.rstrip() for line in fin)
	lines = list(line for line in lines if line) # Non-blank lines in a list   
	valse = lines[1].split()
	valse.append(valse.pop(0))
	del valse[-1]

    if (stage == 'fe' and comp != 'c' and comp != 'w' and comp != 'f'):
      if comp != 'r':
        vals[9+nd] = float(vals[9+nd]) + trans_dist
      shutil.copy('../../../../prep/'+pose+'/assign-eq.dat', './')
      with open('./assign-eq.dat') as fin:
	lines = (line.rstrip() for line in fin)
	lines = list(line for line in lines if line) # Non-blank lines in a list   
	valse = lines[1].split()
	valse.append(valse.pop(0))
	del valse[-1]

    # Define spring constants based on stage and weight 
    if stage == 'equil':
      rdf = rest[0]
      raf = rest[1]
      rdhf = rest[2]
      rdsf = rest[3]
      ldf = 2*weight*rest[4]/100
      laf = weight*rest[5]/100
      ldhf = 0    
      ldsf = weight*rest[7]/100
    elif comp == 's' or comp == 'w' or comp == 'f':
      rdf = rest[0]
      raf = rest[1]
      rdhf = rest[2]
      rdsf = rest[3]
      ldf = rest[4]
      laf = rest[5]
      ldhf = rest[6]   
      ldsf = rest[7]
    elif comp == 'a' or comp == 'r':
      rdf = rest[0]
      raf = rest[1]
      rdhf = weight*rest[2]/100
      rdsf = weight*rest[3]/100
      ldf = 0
      laf = 0
      ldhf = 0    
      ldsf = 0
    elif comp == 'l' or comp == 'c':
      rdf = rest[0]
      raf = rest[1]
      rdhf = rest[2]
      rdsf = rest[3]
      ldf = 0
      laf = 0
      ldhf = weight*rest[6]/100   
      ldsf = weight*rest[7]/100
    elif comp == 't':
      rdf = rest[0]
      raf = rest[1]
      rdhf = rest[2]
      rdsf = rest[3]
      ldf = weight*rest[4]/100
      laf = weight*rest[5]/100
      ldhf = rest[6]   
      ldsf = rest[7]
      
    # Write AMBER restraint file for the full system 
    if (comp != 'c' and comp != 'r' and comp != 'w' and comp != 'f'): 
      disang_file = open('disang.rest', 'w')
      disang_file.write('%s  %s  %s  %s  %s  %s  %s  %s  %s \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage, 'weight = '+str(weight)))
      for i in range(0, len(rst)):
	data = rst[i].split()
	# Protein translational/rotational restraints
	if i < 6 and comp != 'a':
	  if len(data) == 2:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','   
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), rdf, rdf, recep_tr))
	  elif len(data) == 3:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','  
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(180.0), raf, raf, recep_tr))
	  elif len(data) == 4:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, raf, raf, recep_tr))
	# Protein conformation (P1-P3 distance restraints)
	elif i >= 6 and i < 9: 
	  if (stage != 'equil'):
	    if len(data) == 2:
	      nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','   
	      disang_file.write('%s %-23s '%('&rst iat=', nums))
	      disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(valse[i]), float(valse[i]), float(999.0), rdsf, rdsf, recep_c))
	  else:
	    if len(data) == 2:
	      nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','   
	      disang_file.write('%s %-23s '%('&rst iat=', nums))
	      disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), rdsf, rdsf, recep_c))
	# Protein conformation (backbone restraints)
	elif i >= 9 and i < 9+nd: 
	  if (stage != 'equil'):
	    if len(data) == 4:
	      nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
	      disang_file.write('%s %-23s '%('&rst iat=', nums))
	      disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(valse[i]) - 180, float(valse[i]), float(valse[i]), float(valse[i]) + 180 , rdhf, rdhf, recep_d))
	  else:
	    if len(data) == 4:
	      nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
	      disang_file.write('%s %-23s '%('&rst iat=', nums))
	      disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180 , rdhf, rdhf, recep_d))
	# Ligand translational/rotational restraints
	elif i >= 9+nd and i < 15+nd and comp != 'a': 
	  if len(data) == 2:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','   
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldf, ldf, lign_tr))
	  elif len(data) == 3:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','  
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(180.0), laf, laf, lign_tr))
	  elif len(data) == 4:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf, laf, lign_tr))      
	# Ligand conformation (L1-L3 distance restraints)
	elif i >= 15+nd and i < 18+nd and comp != 'a': 
	  if len(data) == 2:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','   
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldsf, ldsf, lign_c))
	# Ligand conformation (non-hydrogen dihedrals)
	elif i >= 18+nd and comp != 'a': 
	  if len(data) == 4:
	    nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
	    disang_file.write('%s %-23s '%('&rst iat=', nums))
	    disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))

      # Analysis of simulations

      if (comp != 'l' and comp != 'a'):
	restraints_file = open('restraints.in', 'w')
	restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage))
	restraints_file.write('noexitonerror\n')
	restraints_file.write('parm vac.prmtop\n')
	for i in range(2,11):
	  restraints_file.write('trajin md%02.0f.nc\n' % i)
	for i in range(9+nd, 15+nd):
	  arr = rst[i].split()
	  if len(arr) == 2:
	    restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
	  if len(arr) == 3:
	    restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
	  if len(arr) == 4:
	    restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
      elif (comp == 'a'):
	restraints_file = open('restraints.in', 'w')
	restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage))
	restraints_file.write('noexitonerror\n')
	restraints_file.write('parm vac.prmtop\n')
	for i in range(2,11):
	  restraints_file.write('trajin md%02.0f.nc\n' % i)
	for i in range(6, 9+nd):
	  arr = rst[i].split()
	  if len(arr) == 2:
	    restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
	  if len(arr) == 3:
	    restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
	  if len(arr) == 4:
	    restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
      elif (comp == 'l'):
	restraints_file = open('restraints.in', 'w')
	restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage))
	restraints_file.write('noexitonerror\n')
	restraints_file.write('parm vac.prmtop\n')
	for i in range(2,11):
	  restraints_file.write('trajin md%02.0f.nc\n' % i)
	for i in range(15+nd, len(rst)):
	  arr = rst[i].split()
	  if len(arr) == 2:
	    restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
	  if len(arr) == 3:
	    restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
	  if len(arr) == 4:
	    restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
    elif comp == 'c' or comp == 'w' or comp == 'f':
      while '' in rst:
        rst.remove('')
    # Write restraint file for ligand system 
      disang_file = open('disang.rest', 'w')
      disang_file.write('%s  %s  %s  %s  %s  %s  %s  %s  %s \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage, 'weight = '+str(weight)))
      for i in range(0, len(rst)):
	data = rst[i].split()
	# Ligand conformational restraints
	if len(data) == 2:
	  nums = str(ligand_atm_num.index(data[0]))+','+str(ligand_atm_num.index(data[1]))+','   
	  disang_file.write('%s %-23s '%('&rst iat=', nums))
	  disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldsf, ldsf, lign_c))
	elif len(data) == 4:
	  nums = str(ligand_atm_num.index(data[0]))+','+str(ligand_atm_num.index(data[1]))+','+str(ligand_atm_num.index(data[2]))+','+str(ligand_atm_num.index(data[3]))+','  
	  disang_file.write('%s %-23s '%('&rst iat=', nums))
	  disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
      # Analysis of simulations
      restraints_file = open('restraints.in', 'w')
      restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage))
      restraints_file.write('noexitonerror\n')
      restraints_file.write('parm vac.prmtop\n')
      for i in range(2,11):
	restraints_file.write('trajin md%02.0f.nc\n' % i)
      for i in range(0, len(rst)):
	arr = rst[i].split()
	if len(arr) == 2:
	  restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
	if len(arr) == 3:
	  restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
	if len(arr) == 4:
	  restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
    elif comp == 'r':
      while '' in rst:
        rst.remove('')
    # Write restraint file for protein system 
      disang_file = open('disang.rest', 'w')
      disang_file.write('%s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', P1, P2, P3, 'stage = '+stage, 'weight = '+str(weight)))
      for i in range(0, len(rst)):
	data = rst[i].split()
	# Protein conformational restraints
	if len(data) == 2:
	  nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','
	  disang_file.write('%s %-23s '%('&rst iat=', nums))
	  disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(valse[i+6]), float(valse[i+6]), float(999.0), rdsf, rdsf, recep_c))
	if len(data) == 4:
	  nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','
	  disang_file.write('%s %-23s '%('&rst iat=', nums))
	  disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(valse[i+6]) - 180, float(valse[i+6]), float(valse[i+6]), float(valse[i+6]) + 180 , rdhf, rdhf, recep_d))
      # Analysis of simulations
      restraints_file = open('restraints.in', 'w')
      restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage))
      restraints_file.write('noexitonerror\n')
      restraints_file.write('parm vac.prmtop\n')
      for i in range(2,11):
	restraints_file.write('trajin md%02.0f.nc\n' % i)
      for i in range(0, len(rst)):
	arr = rst[i].split()
	if len(arr) == 2:
	  restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
	if len(arr) == 3:
	  restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
	if len(arr) == 4:
	  restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))


    if stage != 'fe':
      os.chdir('../')

def sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, rng):

    if stage != 'fe':
      # Copy folder for running simulations
      if not os.path.exists('run_files'):
	try:
	  shutil.copytree('../run_files', './run_files')
	# Directories are the same
	except shutil.Error as e:
	  print('Directory not copied. Error: %s' % e)
	# Any error saying that the directory doesn't exist
	except OSError as e:
	  print('Directory not copied. Error: %s' % e)
      if hmr == 'no': 
        replacement = 'full.prmtop'
        for dname, dirs, files in os.walk('./run_files'):
          for fname in files:
            fpath = os.path.join(dname, fname)
            with open(fpath) as f:
              s = f.read()
              s = s.replace('full.hmr.prmtop', replacement)
            with open(fpath, "w") as f:
              f.write(s)
      elif hmr == 'yes': 
        replacement = 'full.hmr.prmtop'
        for dname, dirs, files in os.walk('./run_files'):
          for fname in files:
            fpath = os.path.join(dname, fname)
            with open(fpath) as f:
              s = f.read()
              s = s.replace('full.prmtop', replacement)
            with open(fpath, "w") as f:
              f.write(s)
      # Change to simulation directory
      os.chdir(pose)

    # Find anchors
    with open('disang.rest', 'r') as f:
	data = f.readline().split()    
	L1 = data[6].strip()   
	L2 = data[7].strip()   
	L3 = data[8].strip()   


    # Get number of atoms in vacuum
    with open('./vac.pdb') as myfile:
	data = myfile.readlines()
	vac_atoms = data[-3][6:11].strip()

    # Create minimization and NPT equilibration files for big box and small ligand box (APR)
    if comp != 'c' and comp != 'r':
      with open("../amber_files/mini.in", "rt") as fin:
        with open("./mini.in", "wt") as fout:
	  for line in fin:
	    fout.write(line.replace('_L1_', L1).replace('_L2_', L2).replace('_L3_', L3))
      with open("../amber_files/therm1.in", "rt") as fin:
        with open("./therm1.in", "wt") as fout:
	  for line in fin:
	    fout.write(line.replace('_L1_', L1).replace('_L2_', L2).replace('_L3_', L3))
      with open("../amber_files/therm2.in", "rt") as fin:
        with open("./therm2.in", "wt") as fout:
	  for line in fin:
	    fout.write(line.replace('_L1_', L1).replace('_L2_', L2).replace('_L3_', L3).replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
	  for line in fin:
	    fout.write(line.replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
	  for line in fin:
	    fout.write(line.replace('_temperature_', str(temperature)))
    else:
      with open("../amber_files/mini.in", "rt") as fin:
        with open("./mini.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line)
      with open("../amber_files/therm1.in", "rt") as fin:
        with open("./therm1.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line)
      with open("../amber_files/therm2.in", "rt") as fin:
        with open("./therm2.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line.replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt-lig.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line.replace('_temperature_', str(temperature)))

    # Create gradual release files for equilibrium
    if (stage == 'equil'):
      for i in range(0, num_sim):
	with open('../amber_files/mdin-equil', "rt") as fin:
	  with open("./mdin-%02d" %int(i), "wt") as fout:
	    if i == (num_sim-1):
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang%02d' %int(i)))
	    else:
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('disang_file', 'disang%02d' %int(i)))


    # Create preparation files
    if (stage == 'prep'):
      for i in range(0, num_sim):
	with open('../amber_files/mdin-prep', "rt") as fin:
	  with open("./mdin-%03d" %int(i), "wt") as fout:
	    if i == 0:
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('disang_file', 'disang%03d' %int(i)))
	    else:
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang%03d' %int(i)))


    # Create free energy files
    if (stage == 'fe'):
      if (comp != 'c' and comp != 'r'):
	for i in range(0, num_sim+1):
	  with open('../amber_files/mdin-apr', "rt") as fin:
	    with open("./mdin-%02d" %int(i), "wt") as fout:
	      if i == 1 or i == 0:
		for line in fin:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('disang_file', 'disang'))
	      else:
		for line in fin:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang'))
      else:
	for i in range(0, num_sim+1):
	  with open('../amber_files/mdin-lig', "rt") as fin:
	    with open("./mdin-%02d" %int(i), "wt") as fout:
	      if i == 1 or i == 0:
		for line in fin:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('disang_file', 'disang'))
	      else:
		for line in fin:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang'))

    # Create running scripts for local and server
    if (stage == 'fe'):
      if (comp != 'c' and comp != 'r'): 
	with open('../run_files/local-'+stage+'.bash', "rt") as fin:
	  with open("./run-local.bash", "wt") as fout:
	    for line in fin:
	      fout.write(line)
	with open('../run_files/PBS-'+stage, "rt") as fin:
	  with open("./PBS-run", "wt") as fout:
	    for line in fin:
    	      fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
      else:
	with open('../run_files/local-lig.bash', "rt") as fin:
	  with open("./run-local.bash", "wt") as fout:
	    for line in fin:
	      fout.write(line)
	with open('../run_files/PBS-lig', "rt") as fin:
	  with open("./PBS-run", "wt") as fout:
	    for line in fin:
    	      fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
    else:
      with open('../run_files/local-'+stage+'.bash', "rt") as fin:
	with open("./run-local.bash", "wt") as fout:
	  for line in fin:
	    fout.write(line.replace('RANGE', str(rng)))
      with open('../run_files/PBS-'+stage, "rt") as fin:
	with open("./PBS-run", "wt") as fout:
	  for line in fin:
	    fout.write(line.replace('STAGE', stage).replace('POSE', pose).replace('RANGE', str(rng)))

    os.chdir('../')


def dec_files(temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, weight, lambdas):

    # Find anchors
    with open('disang.rest', 'r') as f:
	data = f.readline().split()    
	L1 = data[6].strip()   
	L2 = data[7].strip()   
	L3 = data[8].strip()   

    # Get number of atoms in vacuum
    with open('./vac.pdb') as myfile:
	data = myfile.readlines()
	vac_atoms = data[-3][6:11].strip()


    if (comp == 'v'):
      for i in range(0, num_sim+1):
	with open('./vac.pdb') as myfile:
	  data = myfile.readlines()
	  mask = int(data[-3][22:26].strip())
	with open('../amber_files/mdin-lj', "rt") as fin:
	  with open("./mdin-%02d" %int(i), "wt") as fout:
	    if i == 1 or i == 0:
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mask)))
	    else:
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mask)))
	mdin = open("./mdin-%02d" %int(i), 'a')
	mdin.write('  mbar_states = %02d\n' %len(lambdas))
	mdin.write('  mbar_lambda = ')
	for i in range(0, len(lambdas)):
	  mdin.write(' %6.5f,' %(lambdas[i]))
	mdin.write('\n')
	mdin.write(' /\n')
	mdin.write(' &wt type = \'END\' , /\n')
	mdin.write('DISANG=disang.rest\n')
	mdin.write('LISTOUT=POUT\n')

    if (comp == 'e'):
      for i in range(0, num_sim+1):
	with open('./vac.pdb') as myfile:
	  data = myfile.readlines()
	  mk2 = int(data[-3][22:26].strip())
          mk1 = int(mk2 - 1)
	with open('../amber_files/mdin-ch', "rt") as fin:
	  with open("./mdin-%02d" %int(i), "wt") as fout:
	    if i == 1 or i == 0:
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
	    else:
	      for line in fin:
		fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
	mdin = open("./mdin-%02d" %int(i), 'a')
	mdin.write('  mbar_states = %02d\n' %len(lambdas))
	mdin.write('  mbar_lambda = ')
	for i in range(0, len(lambdas)):
	  mdin.write(' %6.5f,' %(lambdas[i]))
	mdin.write('\n')
	mdin.write(' /\n')
	mdin.write(' &wt type = \'END\' , /\n')
	mdin.write('DISANG=disang.rest\n')
	mdin.write('LISTOUT=POUT\n')

    if (comp == 'w'):
      for i in range(0, num_sim+1):
	mask = '1'
	with open('../amber_files/mdin-lj', "rt") as fin:
	  with open("./mdin-%02d" %int(i), "wt") as fout:
	    if i == 1 or i == 0:
	      for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mask)))
	    else:
	      for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mask)))
	mdin = open("./mdin-%02d" %int(i), 'a')
	mdin.write('  mbar_states = %02d\n' %len(lambdas))
	mdin.write('  mbar_lambda = ')
	for i in range(0, len(lambdas)):
	  mdin.write(' %6.5f,' %(lambdas[i]))
	mdin.write('\n')
	mdin.write(' /\n')
	mdin.write(' &wt type = \'END\' , /\n')
	mdin.write('DISANG=disang.rest\n')
	mdin.write('LISTOUT=POUT\n')

    if (comp == 'f'):
      for i in range(0, num_sim+1):
	mk1 = '1'
	mk2 = '2'
	with open('../amber_files/mdin-ch', "rt") as fin:
	  with open("./mdin-%02d" %int(i), "wt") as fout:
	    if i == 1 or i == 0:
	      for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps1)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
	    else:
	      for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line:
		  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
	mdin = open("./mdin-%02d" %int(i), 'a')
	mdin.write('  mbar_states = %02d\n' %len(lambdas))
	mdin.write('  mbar_lambda = ')
	for i in range(0, len(lambdas)):
	  mdin.write(' %6.5f,' %(lambdas[i]))
	mdin.write('\n')
	mdin.write(' /\n')
	mdin.write(' &wt type = \'END\' , /\n')
	mdin.write('DISANG=disang.rest\n')
	mdin.write('LISTOUT=POUT\n')

    if (comp == 'v'):
      # Create running scripts for local and server
      with open('../run_files/local-'+stage+'.bash', "rt") as fin:
	with open("./run-local.bash", "wt") as fout:
	  for line in fin:
	    fout.write(line)
      with open('../run_files/PBS-'+stage, "rt") as fin:
	with open("./PBS-run", "wt") as fout:
	  for line in fin:
    	    fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))

    if (comp == 'e'):
      # Create heating and NPT equilibration files for charge decoupling
      with open("../amber_files/eqnpt-ch.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
	  for line in fin: 
	    fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
      with open("../amber_files/heat-ch.in", "rt") as fin:
        with open("./heat.in", "wt") as fout:
	  for line in fin: 
	    fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
      # Create running scripts for local and server
      with open('../run_files/local-ch.bash', "rt") as fin:
	with open("./run-local.bash", "wt") as fout:
	  for line in fin:
	    fout.write(line)
      with open('../run_files/PBS-ch', "rt") as fin:
	with open("./PBS-run", "wt") as fout:
	  for line in fin:
    	    fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))


    if (comp == 'w'):
      # Create minimization and NPT equilibration files for big box and small ligand box (APR)
      with open("../amber_files/mini.in", "rt") as fin:
        with open("./mini.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line)
      with open("../amber_files/therm1.in", "rt") as fin:
        with open("./therm1.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line)
      with open("../amber_files/therm2.in", "rt") as fin:
        with open("./therm2.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line.replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt-lig.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line.replace('_temperature_', str(temperature)))
      # Create running scripts for local and server
      with open('../run_files/local-lig.bash', "rt") as fin:
	with open("./run-local.bash", "wt") as fout:
	  for line in fin:
	    fout.write(line)
      with open('../run_files/PBS-lig', "rt") as fin:
	with open("./PBS-run", "wt") as fout:
	  for line in fin:
    	    fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))

    if (comp == 'f'):
      # Create minimization and NPT equilibration files for big box and small ligand box (APR)
      mk1 = '1'
      mk2 = '2'
      with open("../amber_files/heat-ch.in", "rt") as fin:
        with open("./heat.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
      with open("../amber_files/eqnpt-ch.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
	  for line in fin: 
            if not 'restraint' in line and not 'ntr = 1' in line:
	      fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
      # Create running scripts for local and server
      with open('../run_files/local-ch.bash', "rt") as fin:
	with open("./run-local.bash", "wt") as fout:
	  for line in fin:
	    fout.write(line)
      with open('../run_files/PBS-ch', "rt") as fin:
	with open("./PBS-run", "wt") as fout:
	  for line in fin:
    	    fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))

    os.chdir('../')


