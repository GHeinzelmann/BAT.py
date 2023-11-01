#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
from lib import scripts as scripts

def restraints(pose, rest, bb_start, bb_end, weight, stage, mol, comp, bb_equil, sdr_dist, dec_method, other_mol):

    rst = []
    atm_num = []
    mlines = []
    hvy_h = []
    hvy_g = []
    msk = []
    pdb_file = ('vac.pdb')
    ligand_pdb_file = ('vac_ligand.pdb')

    if comp == 'n':
      dec_method == 'sdr'

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
      p1_res = P1.split('@')[0][1:]  
      p2_res = P2.split('@')[0][1:]    
      p3_res = P3.split('@')[0][1:]    
      p1_atom = P1.split('@')[1]
      p2_atom = P2.split('@')[1]
      p3_atom = P3.split('@')[1]
      L1 = data[5].strip()   
      L2 = data[6].strip()   
      L3 = data[7].strip()   
      l1_atom = L1.split('@')[1]
      l2_atom = L2.split('@')[1]
      l3_atom = L3.split('@')[1]
      lig_res = L1.split('@')[0][1:]
      first_res = data[8].strip()        
      recep_last = data[9].strip()
 
    # Get backbone atoms and adjust anchors 
    if (comp != 'c' and comp != 'r' and comp != 'f' and comp != 'w'):

      # Get protein backbone atoms
      with open('./vac.pdb') as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line) # Non-blank lines in a list   
        for i in range(0, len(lines)):
          if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
            if int(lines[i][22:26].strip()) >= 2 and int(lines[i][22:26].strip()) < int(lig_res):
              data = lines[i][12:16].strip()
              if data == 'CA' or data == 'N' or data == 'C' or data == 'O':
                hvy_h.append(lines[i][6:11].strip())


      if dec_method == 'sdr':
        if (comp == 'e' or comp == 'v' or comp == 'n'):

          rec_res = int(recep_last) + 2
          lig_res = str((int(lig_res) + 1))
          L1 = ':'+lig_res+'@'+l1_atom
          L2 = ':'+lig_res+'@'+l2_atom
          L3 = ':'+lig_res+'@'+l3_atom
          hvy_h = []
          hvy_g = []

          # Adjust anchors

          p1_resid = str(int(p1_res) + 1)
          p2_resid = str(int(p2_res) + 1)
          p3_resid = str(int(p3_res) + 1)

          P1 = ":"+p1_resid+"@"+p1_atom
          P2 = ":"+p2_resid+"@"+p2_atom
          P3 = ":"+p3_resid+"@"+p3_atom

          # Get receptor heavy atoms
          with open('./vac.pdb') as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if line) # Non-blank lines in a list   
            for i in range(0, len(lines)):
              if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                if int(lines[i][22:26].strip()) >= 3 and int(lines[i][22:26].strip()) <= rec_res:
                  data = lines[i][12:16].strip()
                  if data == 'CA' or data == 'N' or data == 'C' or data == 'O':
                    hvy_h.append(lines[i][6:11].strip())

          # Get bulk ligand heavy atoms
          with open('./vac.pdb') as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if line) # Non-blank lines in a list   
            if comp == 'e':
              for i in range(0, len(lines)):
                if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                  if lines[i][22:26].strip() == str(int(lig_res) + 2):
                    data = lines[i][12:16].strip()
                    if data[0] != 'H':
                      hvy_g.append(lines[i][6:11].strip())
            if comp == 'v':
              for i in range(0, len(lines)):
                if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                  if lines[i][22:26].strip() == str(int(lig_res) + 1):
                    data = lines[i][12:16].strip()
                    if data[0] != 'H':
                      hvy_g.append(lines[i][6:11].strip())
            if comp == 'n':
              for i in range(0, len(lines)):
                if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                  if lines[i][22:26].strip() == str(int(lig_res)):
                    data = lines[i][12:16].strip()
                    if data[0] != 'H':
                      hvy_g.append(lines[i][6:11].strip())


    # Adjust anchors for ligand only
    if (comp == 'c'  or comp == 'w' or comp == 'f'):
      L1 = L1.replace(':'+lig_res, ':1') 
      L2 = L2.replace(':'+lig_res, ':1') 
      L3 = L3.replace(':'+lig_res, ':1') 

    # Get a relation between atom number and masks
    atm_num = scripts.num_to_mask(pdb_file)
    ligand_atm_num = scripts.num_to_mask(ligand_pdb_file)

    # Get number of ligand atoms
    with open('./vac_ligand.pdb') as myfile:
        data = myfile.readlines()
        vac_atoms = int(data[-3][6:11].strip())

    # Define anchor atom distance restraints on the protein

    rst.append(''+P1+' '+P2+'') 
    rst.append(''+P2+' '+P3+'') 
    rst.append(''+P3+' '+P1+'') 

    # Define protein dihedral restraints in the given range
    nd = 0
    for i in range(0, len(bb_start)):
      beg = bb_start[i] - int(first_res) + 2 
      end = bb_end[i] - int(first_res) + 2
      if dec_method == 'sdr':
        if (comp == 'e' or comp == 'v' or comp == 'n'):
          beg = bb_start[i] - int(first_res) + 3
          end = bb_end[i] - int(first_res) + 3
      for i in range(beg, end):
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

    rst.append(''+P1+' '+L1+'')
    rst.append(''+P2+' '+P1+' '+L1+'') 
    rst.append(''+P3+' '+P2+' '+P1+' '+L1+'') 
    rst.append(''+P1+' '+L1+' '+L2+'') 
    rst.append(''+P2+' '+P1+' '+L1+' '+L2+'') 
    rst.append(''+P1+' '+L1+' '+L2+' '+L3+'') 

    # New restraints for ligand only
    if (comp == 'c'  or comp == 'w' or comp == 'f'):
      rst = []

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
          anum.append(abs(int(data[j])//3)+1)      
        msk.append('%s %s %s %s' %(ligand_atm_num[anum[0]], ligand_atm_num[anum[1]], ligand_atm_num[anum[2]], ligand_atm_num[anum[3]]))   
           
    for i in range(0, len(mlines)):
      data = mlines[i].split()
      if len(data) > 7:
        if int(data[8]) > 0:
          anum = []
          for j in range (0, len(data)): 
            anum.append(abs(int(data[j])//3)+1)      
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
      msk = list(filter(None, msk)) 
      msk = [m.replace(':1',':'+lig_res) for m in msk]
    else:
      msk = list(filter(None, msk))

    # Remove dihedral restraints on sp carbons to avoid crashes 
    sp_carb = []
    with open('./'+mol.lower()+'.mol2') as fin:
      lines = (line.rstrip() for line in fin)
      lines = list(line for line in lines if line) # Non-blank lines in a list   
      for line in lines:
        data = line.split()
        if len(data) > 6:
          if data[5] == 'cg' or data[5] == 'c1':
            sp_carb.append(data[1])
    for i in range(0, len(msk)):
      rem_dih = 0
      data = msk[i].split()
      for j in range(0, len(sp_carb)):
        atom_name1 = data[1].split('@')[1]
        atom_name2 = data[2].split('@')[1]
        if atom_name1 == sp_carb[j] or atom_name2 == sp_carb[j]:
          rem_dih = 1
          break
      if rem_dih == 0:
        rst.append(msk[i])

    # New restraints for protein only
    if (comp == 'r'):
      rst = []
      rst.append(''+P1+' '+P2+'') 
      rst.append(''+P2+' '+P3+'') 
      rst.append(''+P3+' '+P1+'') 
      nd = 0
      for i in range(beg, end):
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

    # If chosen, apply initial reference for the protein backbone restraints
    if (stage == 'fe' and comp != 'c' and comp != 'w' and comp != 'f'):
      if (bb_equil == 'yes'):
        shutil.copy('../../../../equil/'+pose+'/assign.dat', './assign-eq.dat')
      else:
        shutil.copy('./assign.dat', './assign-eq.dat')
      with open('./assign-eq.dat') as fin:
        lines = (line.rstrip() for line in fin)
        lines = list(line for line in lines if line) # Non-blank lines in a list   
        valse = lines[1].split()
        valse.append(valse.pop(0))
        del valse[-1]

    # Define spring constants based on stage and weight 
    if stage == 'equil':
      if bb_equil == 'yes':
        rdhf = rest[0]
      else:
        rdhf = 0
      rdsf = rest[1]
      ldf =  weight*rest[2]/100
      laf =  weight*rest[3]/100
      ldhf = weight*rest[4]/100
      rcom = rest[5]
    elif comp == 'l' or comp == 'c':
      rdhf = rest[0]
      rdsf = rest[1]
      ldf = 0
      laf = 0
      ldhf = weight*rest[4]/100
      rcom = rest[5]
    elif comp == 'a' or comp == 'r':
      rdhf = weight*rest[0]/100
      rdsf = weight*rest[1]/100
      ldf = 0
      laf = 0
      ldhf = 0
      rcom = rest[5]
    elif comp == 't':
      rdhf = rest[0]
      rdsf = rest[1]
      ldf = weight*rest[2]/100
      laf = weight*rest[3]/100
      ldhf = rest[4]
      rcom = rest[5]
    elif comp == 'm':
      rdhf = weight*rest[0]/100
      rdsf = weight*rest[1]/100
      ldf = weight*rest[2]/100
      laf = weight*rest[3]/100
      ldhf = weight*rest[4]/100
      rcom = rest[5]
    elif comp == 'n':
      rdhf = weight*rest[0]/100
      rdsf = weight*rest[1]/100
      ldf = 0
      laf = 0
      ldhf = weight*rest[4]/100
      rcom = rest[5]
      lcom = rest[6]
    elif comp == 'v' or comp == 'e' or comp == 'w' or comp == 'f':
      rdhf = rest[0]
      rdsf = rest[1]
      ldf = rest[2]
      laf = rest[3]
      ldhf = rest[4]
      rcom = rest[5]
      lcom = rest[6]

      
    # Write AMBER restraint file for the full system 
    if (comp != 'c' and comp != 'r' and comp != 'w' and comp != 'f'): 
      disang_file = open('disang.rest', 'w')
      disang_file.write('%s  %s  %s  %s  %s  %s  %s  %s  %s \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage, 'weight = '+str(weight)))
      for i in range(0, len(rst)):
        data = rst[i].split()
        # Protein conformation (P1-P3 distance restraints)
        if i < 3: 
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
        elif i >= 3 and i < 3+nd: 
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
        elif i >= 3+nd and i < 9+nd and comp != 'a': 
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
          if comp == 'e':
            if i == (3+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldf, ldf, lign_tr))
            if i == (4+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2])+vac_atoms)+','
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(180.0), laf, laf, lign_tr))
            if i == (5+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3])+vac_atoms)+','  
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf, laf, lign_tr))      
            if i == (6+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','  
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(180.0), laf, laf, lign_tr))
            if i == (7+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','  
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf, laf, lign_tr))      
            if i == (8+nd):
              nums2 = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','  
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, laf, laf, lign_tr))      
        # Ligand conformation (non-hydrogen dihedrals)
        elif i >= 9+nd and comp != 'a': 
          if len(data) == 4:
            nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','  
            disang_file.write('%s %-23s '%('&rst iat=', nums))
            disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
            if comp == 'v' and dec_method == 'sdr':
              nums2 = str(atm_num.index(data[0])+vac_atoms)+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','  
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
            if comp == 'e':
              nums2 = str(atm_num.index(data[0])+vac_atoms)+','+str(atm_num.index(data[1])+vac_atoms)+','+str(atm_num.index(data[2])+vac_atoms)+','+str(atm_num.index(data[3])+vac_atoms)+','  
              disang_file.write('%s %-23s '%('&rst iat=', nums2))
              disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
              if (dec_method == 'sdr'):
                nums3 = str(atm_num.index(data[0])+2*vac_atoms)+','+str(atm_num.index(data[1])+2*vac_atoms)+','+str(atm_num.index(data[2])+2*vac_atoms)+','+str(atm_num.index(data[3])+2*vac_atoms)+','
                disang_file.write('%s %-23s '%('&rst iat=', nums3))
                disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
                nums4 = str(atm_num.index(data[0])+3*vac_atoms)+','+str(atm_num.index(data[1])+3*vac_atoms)+','+str(atm_num.index(data[2])+3*vac_atoms)+','+str(atm_num.index(data[3])+3*vac_atoms)+','
                disang_file.write('%s %-23s '%('&rst iat=', nums4))
                disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))

      # COM restraints
      cv_file = open('cv.in', 'w')
      cv_file.write('cv_file \n')
      cv_file.write('&colvar \n')
      cv_file.write(' cv_type = \'COM_DISTANCE\' \n')
      cv_file.write(' cv_ni = %s, cv_i = 1,0,' % str(len(hvy_h)+2))
      for i in range(0, len(hvy_h)):
        cv_file.write(hvy_h[i])
        cv_file.write(',')
      cv_file.write('\n')
      cv_file.write(' anchor_position = %10.4f, %10.4f, %10.4f, %10.4f \n' % (float(0.0), float(0.0), float(0.0), float(999.0)))
      cv_file.write(' anchor_strength = %10.4f, %10.4f, \n' % (rcom, rcom))
      cv_file.write('/ \n')
      if dec_method == 'sdr':
        if comp == 'e' or comp == 'v' or comp == 'n':
          cv_file.write('&colvar \n')
          cv_file.write(' cv_type = \'COM_DISTANCE\' \n')
          cv_file.write(' cv_ni = %s, cv_i = 2,0,' % str(len(hvy_g)+2))
          for i in range(0, len(hvy_g)):
            cv_file.write(hvy_g[i])
            cv_file.write(',')
          cv_file.write('\n')
          cv_file.write(' anchor_position = %10.4f, %10.4f, %10.4f, %10.4f \n' % (float(0.0), float(0.0), float(0.0), float(999.0)))
          cv_file.write(' anchor_strength = %10.4f, %10.4f, \n' % (lcom, lcom))
          cv_file.write('/ \n')
        cv_file.close()



      # Analysis of simulations

      if (comp != 'l' and comp != 'a' and comp != 'm' and comp != 'n'):  
        restraints_file = open('restraints.in', 'w')
        restraints_file.write('%s  %s  %s  %s  %s  %s  %s  %s  \n'%('# Anchor atoms', P1, P2, P3, L1, L2, L3, 'stage = '+stage))
        restraints_file.write('noexitonerror\n')
        restraints_file.write('parm vac.prmtop\n')
        for i in range(2,11):
          restraints_file.write('trajin md%02.0f.nc\n' % i)
        for i in range(3+nd, 9+nd):
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
        for i in range(0, 3+nd):
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
        for i in range(9+nd, len(rst)):
          arr = rst[i].split()
          if len(arr) == 2:
            restraints_file.write('%s %s %s'%('distance d'+str(i), rst[i], 'noimage out restraints.dat\n'))
          if len(arr) == 3:
            restraints_file.write('%s %s %s'%('angle a'+str(i), rst[i], 'out restraints.dat\n'))
          if len(arr) == 4:
            restraints_file.write('%s %s %s'%('dihedral a'+str(i), rst[i], 'out restraints.dat\n'))
      elif (comp == 'm' or comp == 'n'):
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
    elif comp == 'f':
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
          nums2 = str(ligand_atm_num.index(data[0])+vac_atoms)+','+str(ligand_atm_num.index(data[1])+vac_atoms)+','
          disang_file.write('%s %-23s '%('&rst iat=', nums))
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldsf, ldsf, lign_c))
          disang_file.write('%s %-23s '%('&rst iat=', nums2))
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(vals[i]), float(vals[i]), float(999.0), ldsf, ldsf, lign_c))
        elif len(data) == 4:
          nums = str(ligand_atm_num.index(data[0]))+','+str(ligand_atm_num.index(data[1]))+','+str(ligand_atm_num.index(data[2]))+','+str(ligand_atm_num.index(data[3]))+','
          nums2 = str(ligand_atm_num.index(data[0])+vac_atoms)+','+str(ligand_atm_num.index(data[1])+vac_atoms)+','+str(ligand_atm_num.index(data[2])+vac_atoms)+','+str(ligand_atm_num.index(data[3])+vac_atoms)+','
          disang_file.write('%s %-23s '%('&rst iat=', nums))
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(vals[i]) - 180, float(vals[i]), float(vals[i]), float(vals[i]) + 180, ldhf, ldhf, lign_d))
          disang_file.write('%s %-23s '%('&rst iat=', nums2))
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
    elif comp == 'c'  or comp == 'w':
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
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(0.0), float(valse[i]), float(valse[i]), float(999.0), rdsf, rdsf, recep_c))
        if len(data) == 4:
          nums = str(atm_num.index(data[0]))+','+str(atm_num.index(data[1]))+','+str(atm_num.index(data[2]))+','+str(atm_num.index(data[3]))+','
          disang_file.write('%s %-23s '%('&rst iat=', nums))
          disang_file.write('r1= %10.4f, r2= %10.4f, r3= %10.4f, r4= %10.4f, rk2= %11.7f, rk3= %11.7f, &end %s \n' % (float(valse[i]) - 180, float(valse[i]), float(valse[i]), float(valse[i]) + 180 , rdhf, rdhf, recep_d))
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
 
    disang_file.write('\n')


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

    # Create minimization and NPT equilibration files for big box and small ligand box 
    if comp != 'c' and comp != 'r' and comp != 'n':
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
    elif (comp == 'r' or comp == 'c'):
      with open("../amber_files/mini-lig.in", "rt") as fin:
        with open("./mini.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line)
      with open("../amber_files/therm1-lig.in", "rt") as fin:
        with open("./therm1.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line)
      with open("../amber_files/therm2-lig.in", "rt") as fin:
        with open("./therm2.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line.replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt-lig.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
            if not 'restraint' in line and not 'ntr = 1' in line:
              fout.write(line.replace('_temperature_', str(temperature)))
    else: # n component
      with open("../amber_files/mini-sim.in", "rt") as fin:
        with open("./mini.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_L1_', L1).replace('_L2_', L2).replace('_L3_', L3))
      with open("../amber_files/therm1-sim.in", "rt") as fin:
        with open("./therm1.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_L1_', L1).replace('_L2_', L2).replace('_L3_', L3))
      with open("../amber_files/therm2-sim.in", "rt") as fin:
        with open("./therm2.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_L1_', L1).replace('_L2_', L2).replace('_L3_', L3).replace('_temperature_', str(temperature)))
      with open("../amber_files/eqnpt-sim.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
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


    # Create preparation files (not used anymore)
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
      if (comp != 'c' and comp != 'r' and comp != 'n'):
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-rest', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('disang_file', 'disang'))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang'))
      elif (comp == 'r' or comp == 'c'):
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-lig', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('disang_file', 'disang'))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang'))
      else: # n
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-sim', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('disang_file', 'disang'))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('disang_file', 'disang'))

    # Create running scripts for local and server
    if (stage == 'fe'):
      with open('../run_files/local-lig.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open('../run_files/PBS-Am', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
              fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
      with open('../run_files/SLURMM-Am', "rt") as fin:
        with open("./SLURMM-run", "wt") as fout:
          for line in fin:
              fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
    else:
      with open('../run_files/local-'+stage+'.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line.replace('RANGE', str(rng)))
      with open('../run_files/PBS-Am', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', stage).replace('POSE', pose))
      with open('../run_files/SLURMM-Am', "rt") as fin:
        with open("./SLURMM-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', stage).replace('POSE', pose))

    os.chdir('../')


def dec_files(temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, weight, lambdas, dec_method, ntwx):

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

    # Get last ligand residue number
    with open('./vac.pdb') as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line) # Non-blank lines in a list   
        for i in range(0, len(lines)):
          if lines[i][17:20].strip() == mol.upper():
            last_lig = lines[i][22:26].strip()


    if (comp == 'v'):
      # Create simulation files for vdw decoupling
      if (dec_method == 'sdr'): 
      # Simulation files for simultaneous decoupling
        with open('./vac.pdb') as myfile:
          data = myfile.readlines()
          mk2 = int(last_lig)
          mk1 = int(mk2 - 1)
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-lj', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
          mdin = open("./mdin-%02d" %int(i), 'a')
          mdin.write('  mbar_states = %02d\n' %len(lambdas))
          mdin.write('  mbar_lambda = ')
          for i in range(0, len(lambdas)):
            mdin.write(' %6.5f,' %(lambdas[i]))
          mdin.write('\n')
          mdin.write('  infe = 1,\n')
          mdin.write(' /\n')
          mdin.write(' &pmd \n')
          mdin.write(' output_file = \'cmass.txt\' \n')
          mdin.write(' output_freq = %02d \n' % int(ntwx))
          mdin.write(' cv_file = \'cv.in\' \n')
          mdin.write(' /\n')
          mdin.write(' &wt type = \'END\' , /\n')
          mdin.write('DISANG=disang.rest\n')
          mdin.write('LISTOUT=POUT\n')

        with open("../amber_files/eqnpt-lj.in", "rt") as fin:
          with open("./eqnpt.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
        with open("../amber_files/heat-lj.in", "rt") as fin:
          with open("./heat.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))

      # Simulation files for double decoupling
      elif (dec_method == 'dd'): 
        with open('./vac.pdb') as myfile:
          data = myfile.readlines()
          mk1 = int(last_lig)
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-lj-dd', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))
          mdin = open("./mdin-%02d" %int(i), 'a')
          mdin.write('  mbar_states = %02d\n' %len(lambdas))
          mdin.write('  mbar_lambda = ')
          for i in range(0, len(lambdas)):
            mdin.write(' %6.5f,' %(lambdas[i]))
          mdin.write('\n')
          mdin.write('  infe = 1,\n')
          mdin.write(' /\n')
          mdin.write(' &pmd \n')
          mdin.write(' output_file = \'cmass.txt\' \n')
          mdin.write(' output_freq = %02d \n' % int(ntwx))
          mdin.write(' cv_file = \'cv.in\' \n')
          mdin.write(' /\n')
          mdin.write(' &wt type = \'END\' , /\n')
          mdin.write('DISANG=disang.rest\n')
          mdin.write('LISTOUT=POUT\n')

        with open("../amber_files/eqnpt-lj-dd.in", "rt") as fin:
          with open("./eqnpt.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))
        with open("../amber_files/heat-lj-dd.in", "rt") as fin:
          with open("./heat.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))

      # Create running scripts for local and server
      with open('../run_files/local-dd.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open('../run_files/PBS-Am', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
      with open('../run_files/SLURMM-Am', "rt") as fin:
        with open("./SLURMM-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))

    if (comp == 'e'):
      # Create simulation files for charge decoupling
      if (dec_method == 'sdr'): 
        # Simulation files for simultaneous decoupling
        with open('./vac.pdb') as myfile:
          data = myfile.readlines()
          mk4 = int(last_lig)
          mk3 = int(mk4 - 1)
          mk2 = int(mk4 - 2)
          mk1 = int(mk4 - 3)
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-ch', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))
          mdin = open("./mdin-%02d" %int(i), 'a')
          mdin.write('  mbar_states = %02d\n' %len(lambdas))
          mdin.write('  mbar_lambda = ')
          for i in range(0, len(lambdas)):
            mdin.write(' %6.5f,' %(lambdas[i]))
          mdin.write('\n')
          mdin.write('  infe = 1,\n')
          mdin.write(' /\n')
          mdin.write(' &pmd \n')
          mdin.write(' output_file = \'cmass.txt\' \n')
          mdin.write(' output_freq = %02d \n' % int(ntwx))
          mdin.write(' cv_file = \'cv.in\' \n')
          mdin.write(' /\n')
          mdin.write(' &wt type = \'END\' , /\n')
          mdin.write('DISANG=disang.rest\n')
          mdin.write('LISTOUT=POUT\n')

        with open("../amber_files/eqnpt-ch.in", "rt") as fin:
          with open("./eqnpt.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))
        with open("../amber_files/heat-ch.in", "rt") as fin:
          with open("./heat.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)).replace('mk3',str(mk3)).replace('mk4',str(mk4)))

      elif (dec_method == 'dd'): 
        with open('./vac.pdb') as myfile:
          # Simulation files for double decoupling
          data = myfile.readlines()
          mk2 = int(last_lig)
          mk1 = int(mk2 - 1)
        for i in range(0, num_sim+1):
          with open('../amber_files/mdin-ch-dd', "rt") as fin:
            with open("./mdin-%02d" %int(i), "wt") as fout:
              if i == 1 or i == 0:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
              else:
                for line in fin:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
          mdin = open("./mdin-%02d" %int(i), 'a')
          mdin.write('  mbar_states = %02d\n' %len(lambdas))
          mdin.write('  mbar_lambda = ')
          for i in range(0, len(lambdas)):
            mdin.write(' %6.5f,' %(lambdas[i]))
          mdin.write('\n')
          mdin.write('  infe = 1,\n')
          mdin.write(' /\n')
          mdin.write(' &pmd \n')
          mdin.write(' output_file = \'cmass.txt\' \n')
          mdin.write(' output_freq = %02d \n' % int(ntwx))
          mdin.write(' cv_file = \'cv.in\' \n')
          mdin.write(' /\n')
          mdin.write(' &wt type = \'END\' , /\n')
          mdin.write('DISANG=disang.rest\n')
          mdin.write('LISTOUT=POUT\n')

        with open("../amber_files/eqnpt-ch-dd.in", "rt") as fin:
          with open("./eqnpt.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
        with open("../amber_files/heat-ch-dd.in", "rt") as fin:
          with open("./heat.in", "wt") as fout:
            for line in fin: 
              fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))

      # Create running scripts for local and server
      with open('../run_files/local-dd.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open('../run_files/PBS-Am', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
      with open('../run_files/SLURMM-Am', "rt") as fin:
        with open("./SLURMM-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))

    if (comp == 'f'):
      mk1 = '1'
      mk2 = '2'
      for i in range(0, num_sim+1):
        with open('../amber_files/mdin-ch-dd', "rt") as fin:
          with open("./mdin-%02d" %int(i), "wt") as fout:
            if i == 1 or i == 0:
              for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line and not 'infe' in line:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
            else:
              for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line and not 'infe' in line:
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

      with open("../amber_files/heat-ch-lig.in", "rt") as fin:
        with open("./heat.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))
      with open("../amber_files/eqnpt-ch-lig.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)).replace('mk2',str(mk2)))

      # Create running scripts for local and server
      with open('../run_files/local-dd.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open('../run_files/PBS-Am', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
      with open('../run_files/SLURMM-Am', "rt") as fin:
        with open("./SLURMM-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))

    if (comp == 'w'):
      for i in range(0, num_sim+1):
        mk1 = '1'
        with open('../amber_files/mdin-lj-dd', "rt") as fin:
          with open("./mdin-%02d" %int(i), "wt") as fout:
            if i == 1 or i == 0:
              for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(round(steps1/2))).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))
            else:
              for line in fin:
                if not 'restraint' in line and not 'ntr = 1' in line and not 'ntwprt' in line:
                  fout.write(line.replace('_temperature_', str(temperature)).replace('_num-atoms_', str(vac_atoms)).replace('_num-steps_', str(steps2)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))
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

      with open("../amber_files/heat-lj-lig.in", "rt") as fin:
        with open("./heat.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))
      with open("../amber_files/eqnpt-lj-lig.in", "rt") as fin:
        with open("./eqnpt.in", "wt") as fout:
          for line in fin:
            fout.write(line.replace('_temperature_', str(temperature)).replace('lbd_val', '%6.5f' %float(weight)).replace('mk1',str(mk1)))

      # Create running scripts for local and server
      with open('../run_files/local-dd.bash', "rt") as fin:
        with open("./run-local.bash", "wt") as fout:
          for line in fin:
            fout.write(line)
      with open('../run_files/PBS-Am', "rt") as fin:
        with open("./PBS-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))
      with open('../run_files/SLURMM-Am', "rt") as fin:
        with open("./SLURMM-run", "wt") as fout:
          for line in fin:
            fout.write(line.replace('STAGE', pose).replace('POSE', '%s%02d' %(comp, int(win))))

    os.chdir('../')


