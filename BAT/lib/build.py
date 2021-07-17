#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re as re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
from lib import scripts as scripts

def build_equil(pose, celp_st, mol, H1, H2, H3, calc_type, l1_x, l1_y, l1_z, l1_range, min_adis, max_adis, ligand_ff, ligand_ph):

    # Not apply SDR distance when equilibrating
    sdr_dist = 0

    # Create equilibrium directory
    if not os.path.exists('equil'):
      os.makedirs('equil')
    os.chdir('equil')
    if os.path.exists('./build_files'):
      shutil.rmtree('./build_files')
    try:
      shutil.copytree('../build_files', './build_files')
    # Directories are the same
    except shutil.Error as e:
      print('Directory not copied. Error: %s' % e)
    # Any error saying that the directory doesn't exist
    except OSError as e:
      print('Directory not copied. Error: %s' % e)
    os.chdir('build_files')

    if calc_type == 'dock':
      shutil.copy('../../all-poses/%s_docked.pdb' %(celp_st), './protein.pdb')
      shutil.copy('../../all-poses/%s.pdb' %(pose), './')
    elif calc_type == 'crystal':    
      shutil.copy('../../all-poses/%s.pdb' %(pose), './')
      # Replace names and run initial VMD script
      with open("prep-crystal.tcl", "rt") as fin:
        with open("prep.tcl", "wt") as fout:
          for line in fin:
            fout.write(line.replace('MMM', mol).replace('mmm', mol.lower()).replace('CCCC', pose))
      sp.call('vmd -dispdev text -e prep.tcl', shell=True)
  
 
    # Get beginning and end of protein and save first residue as global variable
    with open('./protein.pdb') as myfile:
        data = myfile.readlines()
        first_res = data[1][22:26].strip()
        last_res = data[-2][22:26].strip()
        print('Receptor first residue: %s' %first_res)
        print('Receptor last residue: %s' %last_res)

    # Get original protein length
    recep_resid_num = int(last_res) - int(first_res) + 1

    # Adjust protein anchors to the new residue numbering
    h1_resid = H1.split('@')[0][1:]
    h2_resid = H2.split('@')[0][1:]
    h3_resid = H3.split('@')[0][1:]

    h1_atom = H1.split('@')[1]
    h2_atom = H2.split('@')[1]
    h3_atom = H3.split('@')[1]
   
    p1_resid = str(int(h1_resid) - int(first_res) + 2)
    p2_resid = str(int(h2_resid) - int(first_res) + 2)
    p3_resid = str(int(h3_resid) - int(first_res) + 2)

    p1_vmd = str(int(h1_resid) - int(first_res) + 1)

    P1 = ":"+p1_resid+"@"+h1_atom
    P2 = ":"+p2_resid+"@"+h2_atom
    P3 = ":"+p3_resid+"@"+h3_atom

    print('Receptor anchors:')
    print(P1)
    print(P2)
    print(P3)

    # Replace names in initial files and VMD scripts
    with open("prep-ini.tcl", "rt") as fin:
      with open("prep.tcl", "wt") as fout:
        for line in fin:
          fout.write(line.replace('MMM', mol).replace('mmm', mol.lower()).replace('NN', h1_atom).replace('P1A', p1_vmd).replace('FIRST','1').replace('LAST',str(recep_resid_num)).replace('STAGE','equil').replace('XDIS','%4.2f' %l1_x).replace('YDIS','%4.2f' %l1_y).replace('ZDIS','%4.2f' %l1_z).replace('RANG','%4.2f' %l1_range).replace('DMAX','%4.2f' %max_adis).replace('DMIN','%4.2f' %min_adis).replace('SDRD','%4.2f' %sdr_dist))
#    with open('%s.pdb' %pose) as f:
#       data=f.read().replace('LIG','%s' %mol)
#    with open('%s.pdb' %pose, "w") as f:
#       f.write(data)


    # Save parameters in ff folder
    if not os.path.exists('../ff/'):
      os.makedirs('../ff/')
    shutil.copy('./dum.mol2', '../ff/')
    shutil.copy('./dum.frcmod', '../ff/')

    # Get parameters and adjust files
    if calc_type == 'dock':
      sp.call('babel -i pdb '+pose+'.pdb -o pdb '+mol.lower()+'.pdb -d', shell=True)
    sp.call('babel -i pdb '+mol.lower()+'.pdb -o pdb '+mol.lower()+'-h.pdb -p %4.2f' %ligand_ph, shell=True)
    sp.call('babel -i pdb '+mol.lower()+'.pdb -o mol2 '+mol.lower()+'-crg.mol2 -p %4.2f' %ligand_ph, shell=True)
    # Get ligand net charge from babel
    lig_crg = 0
    with open('%s-crg.mol2' %mol.lower()) as f_in:
      for line in f_in:
        splitdata = line.split()
        if len(splitdata) > 8:
          lig_crg = lig_crg + float(splitdata[8].strip())
    ligand_charge = round(lig_crg)
    print('The protonation of the ligand is for pH %4.2f' %ligand_ph)
    print('The net charge of the ligand is %d' %ligand_charge)
    if not os.path.exists('../ff/%s.mol2' %mol.lower()):
      print('Antechamber parameters command: antechamber -i '+mol.lower()+'-h.pdb -fi pdb -o '+mol.lower()+'.mol2 -fo mol2 -c bcc -s 2 -at '+ligand_ff.lower()+' -nc %d' % ligand_charge)
      sp.call('antechamber -i '+mol.lower()+'-h.pdb -fi pdb -o '+mol.lower()+'.mol2 -fo mol2 -c bcc -s 2 -at '+ligand_ff.lower()+' -nc %d' % ligand_charge, shell=True)
      shutil.copy('./%s.mol2' %(mol.lower()), '../ff/')
    if not os.path.exists('../ff/%s.frcmod' %mol.lower()):
      if ligand_ff == 'gaff':
        sp.call('parmchk2 -i '+mol.lower()+'.mol2 -f mol2 -o '+mol.lower()+'.frcmod -s 1', shell=True)
      elif ligand_ff == 'gaff2':
        sp.call('parmchk2 -i '+mol.lower()+'.mol2 -f mol2 -o '+mol.lower()+'.frcmod -s 2', shell=True)
      shutil.copy('./%s.frcmod' %(mol.lower()), '../ff/')
    sp.call('antechamber -i '+mol.lower()+'-h.pdb -fi pdb -o '+mol.lower()+'.pdb -fo pdb', shell=True)

    # Create raw complex and clean it
    filenames = ['protein.pdb', '%s.pdb' %mol.lower()]
    with open('./complex-merge.pdb', 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
    with open('complex-merge.pdb') as oldfile, open('complex.pdb', 'w') as newfile:
        for line in oldfile:
            if not 'TER' in line and not 'CONECT' in line:
                newfile.write(line)

    # Align to reference structure using mustang
    sp.call('mustang-3.2.3 -p ./ -i reference.pdb complex.pdb -o aligned -r ON', shell=True)

    # Put in AMBER format and find ligand anchor atoms
    with open('aligned.pdb', 'r') as oldfile, open('aligned-clean.pdb', 'w') as newfile:
        for line in oldfile:
            splitdata = line.split()
            if len(splitdata) > 4:
              if line[21] != 'A':
                newfile.write(line)
    sp.call('pdb4amber -i aligned-clean.pdb -o aligned_amber.pdb -y', shell=True)
    sp.call('vmd -dispdev text -e prep.tcl', shell=True)

    # Check size of anchor file 
    anchor_file = 'anchors.txt'
    if os.stat(anchor_file).st_size == 0:
      os.chdir('../')
      return 'anch1'
    f = open(anchor_file, 'r')
    for line in f:
      splitdata = line.split()
      if len(splitdata) < 3:
        os.rename('./anchors.txt', 'anchors-'+pose+'.txt')
        os.chdir('../')
        return 'anch2'
    os.rename('./anchors.txt', 'anchors-'+pose+'.txt')
    os.chdir('../')

    # Create simulation directory
    if not os.path.exists(pose):
      os.makedirs(pose)
    os.chdir(pose)

    recep_coords = []
    lig_coords = []
    dum_coords = []
    recep_atom = 0
    lig_atom = 0
    dum_atom = 0
    total_atom = 0
    resname_list = []
    resid_list = []
    recep_atomlist = []
    lig_atomlist = []
    dum_atomlist = []

    # Copy a few files
    shutil.copy('../build_files/equil-%s.pdb' %mol.lower(), './')
    shutil.copy('../build_files/%s-noh.pdb' %mol.lower(), './%s.pdb' %mol.lower())
    shutil.copy('../build_files/anchors-'+pose+'.txt', './anchors.txt')

    # Read coordinates for dummy atoms
    for i in range(1, 2):
      shutil.copy('../build_files/dum'+str(i)+'.pdb', './')
      with open('dum'+str(i)+'.pdb') as dum_in:
        lines = (line.rstrip() for line in dum_in)
        lines = list(line for line in lines if line)
        dum_coords.append((float(lines[1][30:38].strip()), float(lines[1][38:46].strip()), float(lines[1][46:54].strip())))
        dum_atomlist.append(lines[1][12:16].strip())
        resname_list.append(lines[1][17:20].strip())
        resid_list.append(float(lines[1][22:26].strip()))
        dum_atom += 1
        total_atom += 1

    # Read coordinates from aligned system
    with open('equil-%s.pdb' %mol.lower()) as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line) # Non-blank lines in a list   

    # Count atoms of receptor and ligand
    for i in range(0, len(lines)):
      if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
              if (lines[i][17:20].strip() != mol) and (lines[i][17:20].strip() != 'DUM'):
                 recep_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                 recep_atomlist.append(lines[i][12:16].strip())
                 resname_list.append(lines[i][17:20].strip())
                 resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                 recep_last = int(lines[i][22:26].strip())
                 recep_atom += 1
                 total_atom += 1
              elif lines[i][17:20].strip() == mol:
                 lig_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                 lig_atomlist.append(lines[i][12:16].strip())
                 resname_list.append(lines[i][17:20].strip())
                 resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                 lig_atom += 1
                 total_atom += 1

    coords = dum_coords + recep_coords + lig_coords
    atom_namelist = dum_atomlist + recep_atomlist + lig_atomlist
    lig_resid = str(recep_last + dum_atom + 1)

    # Read ligand anchors obtained from VMD
    anchor_file = 'anchors.txt'
    f = open(anchor_file, 'r')
    for line in f:
       splitdata = line.split()
       L1 = ":"+lig_resid+"@"+splitdata[0]
       L2 = ":"+lig_resid+"@"+splitdata[1]
       L3 = ":"+lig_resid+"@"+splitdata[2]

    print('Ligand anchors:')
    print(L1)
    print(L2)
    print(L3)
    
    # Write the new pdb file
    build_file = open('build.pdb', 'w')

    # Positions for the dummy atoms 
    for i in range(0, 1):
       build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
       build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))
       build_file.write('%6.2f%6.2f\n'%(0, 0))
       build_file.write('TER\n')

    # Positions of the receptor atoms
    for i in range(dum_atom , dum_atom + recep_atom):
        build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
        build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))

        build_file.write('%6.2f%6.2f\n'%(0, 0))
    build_file.write('TER\n')

    # Positions of the ligand atoms
    for i in range(dum_atom + recep_atom, total_atom):
        build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],mol, float(lig_resid)))
        build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]),float(coords[i][2])))

        build_file.write('%6.2f%6.2f\n'%(0, 0))

    build_file.write('TER\n')
    build_file.write('END\n')
    
    # Write anchors and last protein residue to original pdb file
    with open('equil-%s.pdb' %mol.lower(), 'r') as fin:
        data = fin.read().splitlines(True)
    with open('equil-%s.pdb' %mol.lower(), 'w') as fout:
        fout.write('%-8s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %4s\n' %('REMARK A', P1, P2, P3, L1, L2, L3, first_res, recep_last))
        fout.writelines(data[1:])
    
    # Check for missing residues in receptor structure
    if recep_last != recep_resid_num:
        print('WARNING: Missing residues in the receptor protein sequence. Unless the protein is engineered this is not recommended,') 
        print('a protein modeling tool might be required before running equilibration.') 

    f_in.close()
    build_file.close()

    os.chdir('../')

    return 'all'


def build_rest(hmr, mol, pose, comp, win, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, fwin, l1_x, l1_y, l1_z, l1_range, min_adis, max_adis, sdr_dist):

      
    # Get files or finding new anchors and building some systems

    if not os.path.exists('../build_files'):
      try:
        shutil.copytree('../../../build_files', '../build_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      os.chdir('../build_files')
      # Get last state from equilibrium simulations
      shutil.copy('../../../equil/'+pose+'/md%02d.rst7' %fwin, './')
      for file in glob.glob('../../../equil/%s/full*.prmtop' %pose.lower()):
        shutil.copy(file, './')
      for file in glob.glob('../../../equil/%s/vac*' %pose.lower()):
        shutil.copy(file, './')
      sp.call('cpptraj -p full.prmtop -y md%02d.rst7 -x fe-ini.pdb' %fwin, shell=True)

      # Clean output file
      with open('fe-ini.pdb') as oldfile, open('complex.pdb', 'w') as newfile:
          for line in oldfile:
              if not 'TER' in line and not 'WAT' in line:
                newfile.write(line)


      # Read protein anchors and size from equilibrium
      with open('../../../equil/'+pose+'/equil-%s.pdb' % mol.lower(), 'r') as f:
        data = f.readline().split()
        P1 = data[2].strip()
        P2 = data[3].strip()
        P3 = data[4].strip()
        first_res = data[8].strip()
        recep_last = data[9].strip()

      # Get protein first anchor residue number and protein last residue number from equil simulations
      p1_resid = P1.split('@')[0][1:]
      p1_atom = P1.split('@')[1]
      rec_res = int(recep_last) + 1
      p1_vmd = p1_resid
  

      # Replace names in initial files and VMD scripts
      with open("prep-ini.tcl", "rt") as fin:
        with open("prep.tcl", "wt") as fout:
          for line in fin:
            fout.write(line.replace('MMM', mol).replace('mmm', mol.lower()).replace('NN', p1_atom).replace('P1A', p1_vmd).replace('FIRST','2').replace('LAST',str(rec_res)).replace('STAGE','fe').replace('XDIS','%4.2f' %l1_x).replace('YDIS','%4.2f' %l1_y).replace('ZDIS','%4.2f' %l1_z).replace('RANG','%4.2f' %l1_range).replace('DMAX','%4.2f' %max_adis).replace('DMIN','%4.2f' %min_adis).replace('SDRD','%4.2f' %sdr_dist))

      # Align to reference structure using mustang
      sp.call('mustang-3.2.3 -p ./ -i reference.pdb complex.pdb -o aligned -r ON', shell=True)

      # Put in AMBER format and find ligand anchor atoms
      with open('aligned.pdb', 'r') as oldfile, open('aligned-clean.pdb', 'w') as newfile:
          for line in oldfile:
              splitdata = line.split()
              if len(splitdata) > 4:
                if line[21] != 'A':
                  newfile.write(line)
      sp.call('pdb4amber -i aligned-clean.pdb -o aligned_amber.pdb -y', shell=True)
      sp.call('vmd -dispdev text -e prep.tcl', shell=True)

      # Check size of anchor file 
      anchor_file = 'anchors.txt'
      if os.stat(anchor_file).st_size == 0:
        os.chdir('../')
        return 'anch1'
      f = open(anchor_file, 'r')
      for line in f:
        splitdata = line.split()
        if len(splitdata) < 3:
          os.rename('./anchors.txt', 'anchors-'+pose+'.txt')
          os.chdir('../')
          return 'anch2'
      os.rename('./anchors.txt', 'anchors-'+pose+'.txt')

      # Read ligand anchors obtained from VMD
      lig_resid = str(int(recep_last) + 2)
      anchor_file = 'anchors-'+pose+'.txt'
      f = open(anchor_file, 'r')
      for line in f:
        splitdata = line.split()
        L1 = ":"+lig_resid+"@"+splitdata[0]
        L2 = ":"+lig_resid+"@"+splitdata[1]
        L3 = ":"+lig_resid+"@"+splitdata[2]


      # Write anchors and last protein residue to original pdb file
      with open('fe-%s.pdb' %mol.lower(), 'r') as fin:
          data = fin.read().splitlines(True)
      with open('fe-%s.pdb' %mol.lower(), 'w') as fout:
          fout.write('%-8s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %4s\n' %('REMARK A', P1, P2, P3, L1, L2, L3, first_res, recep_last))
          fout.writelines(data[1:])


      # Get parameters from equilibrium
      if not os.path.exists('../ff'):
        os.makedirs('../ff')
      shutil.copy('../../../equil/ff/%s.mol2' %(mol.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.frcmod' %(mol.lower()), '../ff/')
      shutil.copy('../../../equil/ff/dum.mol2', '../ff/')
      shutil.copy('../../../equil/ff/dum.frcmod', '../ff/')

      os.chdir('../rest/')

    # Copy and replace simulation files for the first window
    if int(win) == 0:
      if os.path.exists('amber_files'):
        shutil.rmtree('./amber_files')
      try:
        shutil.copytree('../../../amber_files', './amber_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      for dname, dirs, files in os.walk('./amber_files'):
        for fname in files:
          fpath = os.path.join(dname, fname)
          with open(fpath) as f:
            s = f.read()
            s = s.replace('_step_', dt).replace('_ntpr_', ntpr).replace('_ntwr_', ntwr).replace('_ntwe_', ntwe).replace('_ntwx_', ntwx).replace('_cutoff_', cut).replace('_gamma_ln_', gamma_ln).replace('_barostat_', barostat).replace('_receptor_ff_', receptor_ff).replace('_ligand_ff_', ligand_ff)
          with open(fpath, "w") as f:
            f.write(s)

      if os.path.exists('run_files'):
        shutil.rmtree('./run_files')
      try:
        shutil.copytree('../../../run_files', './run_files')
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

    if (comp == 'a' or comp == 'l' or comp == 't'):
      # Create window directory
      if not os.path.exists('%s%02d' %(comp, int(win))):
        os.makedirs('%s%02d' %(comp, int(win)))
      os.chdir('%s%02d' %(comp, int(win)))
      # Copy a few files and define new reference state
      if int(win) == 0:
        for file in glob.glob('../../build_files/full*.prmtop'):
          shutil.copy(file, './')
        for file in glob.glob('../../build_files/vac*'):
          shutil.copy(file, './')
        shutil.copy('../../build_files/fe-%s.pdb' %mol.lower(), './')
        shutil.copy('../../build_files/md%02d.rst7' %fwin, './md00.rst7')
        shutil.copy('../../build_files/anchors-'+pose+'.txt', './')
        sp.call('cpptraj -p full.prmtop -y md00.rst7 -x full.rst7 > cpptraj1.log', shell=True)
        shutil.copy('./full.rst7', './full.inpcrd')
        sp.call('cpptraj -p full.prmtop -y md00.rst7 -x full.pdb > cpptraj2.log', shell=True)
      else:
        for file in glob.glob('../%s00/*' %comp):
          shutil.copy(file, './')

    elif comp == 'c':
      # Copy files to c00 to create new box for ligand and copy to the different windows 
      if not os.path.exists('%s%02d' %(comp, int(win))):
        os.makedirs('%s%02d' %(comp, int(win)))
      os.chdir('%s%02d' %(comp, int(win)))
      if int(win) == 0:
        shutil.copy('../../build_files/%s.pdb' %mol.lower(), './')
        shutil.copy('../../build_files/fe-%s.pdb' %mol.lower(), './')
      else:
        for file in glob.glob('../c00/*'):
          shutil.copy(file, './')

    elif comp == 'r':
      # Copy files to r00 to create new box for protein and copy to the different windows 
      if not os.path.exists('%s%02d' %(comp, int(win))):
        os.makedirs('%s%02d' %(comp, int(win)))
      os.chdir('%s%02d' %(comp, int(win)))
      if int(win) == 0:
        shutil.copy('../../build_files/fe-%s.pdb' %mol.lower(), './build-ini.pdb')
        shutil.copy('../../build_files/fe-%s.pdb' %mol.lower(), './')
        shutil.copy('../../build_files/%s.pdb' %mol.lower(), './')
        shutil.copy('../../build_files/full.prmtop', './')
        shutil.copy('../../ff/%s.mol2' %(mol.lower()), './')
        shutil.copy('../../ff/%s.frcmod' %(mol.lower()), './')
        shutil.copy('../../ff/dum.mol2', './')
        shutil.copy('../../ff/dum.frcmod', './')

        # Read coordinates from aligned system
        with open('build-ini.pdb') as f_in:
          lines = (line.rstrip() for line in f_in)
          lines = list(line for line in lines if line) # Non-blank lines in a list   

        # Leave only receptor atoms
        build_file = open('build.pdb', 'w')
        for i in range(0, len(lines)):
          if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
            if (lines[i][17:20].strip() != mol) and (lines[i][17:20].strip() != 'DUM'):
              build_file.write(lines[i]+'\n')
        build_file.close()
      else:
        for file in glob.glob('../r00/*'):
          shutil.copy(file, './')

def build_dec(fwin, hmr, mol, pose, comp, win, water_model, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, sdr_dist, dec_method, l1_x, l1_y, l1_z, l1_range, min_adis, max_adis):


    # Get files or finding new anchors and building some systems

    if (not os.path.exists('../build_files')) or (dec_method == 'sdr' and win == 0):
      if dec_method == 'sdr' and os.path.exists('../build_files'):
        shutil.rmtree('../build_files')    
      try:
        shutil.copytree('../../../build_files', '../build_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      os.chdir('../build_files')
      # Get last state from equilibrium simulations
      shutil.copy('../../../equil/'+pose+'/md%02d.rst7' %fwin, './')
      for file in glob.glob('../../../equil/%s/full*.prmtop' %pose.lower()):
        shutil.copy(file, './')
      for file in glob.glob('../../../equil/%s/vac*' %pose.lower()):
        shutil.copy(file, './')
      sp.call('cpptraj -p full.prmtop -y md%02d.rst7 -x fe-ini.pdb' %fwin, shell=True)

      # Clean output file
      with open('fe-ini.pdb') as oldfile, open('complex.pdb', 'w') as newfile:
          for line in oldfile:
              if not 'TER' in line and not 'WAT' in line:
                newfile.write(line)


      # Read protein anchors and size from equilibrium
      with open('../../../equil/'+pose+'/equil-%s.pdb' % mol.lower(), 'r') as f:
        data = f.readline().split()
        P1 = data[2].strip()
        P2 = data[3].strip()
        P3 = data[4].strip()
        first_res = data[8].strip()
        recep_last = data[9].strip()

      # Get protein first anchor residue number and protein last residue number from equil simulations
      p1_resid = P1.split('@')[0][1:]
      p1_atom = P1.split('@')[1]
      rec_res = int(recep_last) + 1
      p1_vmd = p1_resid
  

      # Replace names in initial files and VMD scripts
      with open("prep-ini.tcl", "rt") as fin:
        with open("prep.tcl", "wt") as fout:
          for line in fin:
            fout.write(line.replace('MMM', mol).replace('mmm', mol.lower()).replace('NN', p1_atom).replace('P1A', p1_vmd).replace('FIRST','2').replace('LAST',str(rec_res)).replace('STAGE','fe').replace('XDIS','%4.2f' %l1_x).replace('YDIS','%4.2f' %l1_y).replace('ZDIS','%4.2f' %l1_z).replace('RANG','%4.2f' %l1_range).replace('DMAX','%4.2f' %max_adis).replace('DMIN','%4.2f' %min_adis).replace('SDRD','%4.2f' %sdr_dist))

      # Align to reference structure using mustang
      sp.call('mustang-3.2.3 -p ./ -i reference.pdb complex.pdb -o aligned -r ON', shell=True)

      # Put in AMBER format and find ligand anchor atoms
      with open('aligned.pdb', 'r') as oldfile, open('aligned-clean.pdb', 'w') as newfile:
          for line in oldfile:
              splitdata = line.split()
              if len(splitdata) > 4:
                if line[21] != 'A':
                  newfile.write(line)
      sp.call('pdb4amber -i aligned-clean.pdb -o aligned_amber.pdb -y', shell=True)
      sp.call('vmd -dispdev text -e prep.tcl', shell=True)

      # Check size of anchor file 
      anchor_file = 'anchors.txt'
      if os.stat(anchor_file).st_size == 0:
        os.chdir('../')
        return 'anch1'
      f = open(anchor_file, 'r')
      for line in f:
        splitdata = line.split()
        if len(splitdata) < 3:
          os.rename('./anchors.txt', 'anchors-'+pose+'.txt')
          os.chdir('../')
          return 'anch2'
      os.rename('./anchors.txt', 'anchors-'+pose+'.txt')

      # Read ligand anchors obtained from VMD
      lig_resid = str(int(recep_last) + 2)
      anchor_file = 'anchors-'+pose+'.txt'
      f = open(anchor_file, 'r')
      for line in f:
        splitdata = line.split()
        L1 = ":"+lig_resid+"@"+splitdata[0]
        L2 = ":"+lig_resid+"@"+splitdata[1]
        L3 = ":"+lig_resid+"@"+splitdata[2]


      # Write anchors and last protein residue to original pdb file
      with open('fe-%s.pdb' %mol.lower(), 'r') as fin:
          data = fin.read().splitlines(True)
      with open('fe-%s.pdb' %mol.lower(), 'w') as fout:
          fout.write('%-8s  %6s  %6s  %6s  %6s  %6s  %6s  %6s  %4s\n' %('REMARK A', P1, P2, P3, L1, L2, L3, first_res, recep_last))
          fout.writelines(data[1:])


      # Get parameters from equilibrium
      if not os.path.exists('../ff'):
        os.makedirs('../ff')
      shutil.copy('../../../equil/ff/%s.mol2' %(mol.lower()), '../ff/')
      shutil.copy('../../../equil/ff/%s.frcmod' %(mol.lower()), '../ff/')
      shutil.copy('../../../equil/ff/dum.mol2', '../ff/')
      shutil.copy('../../../equil/ff/dum.frcmod', '../ff/')

      os.chdir('../'+dec_method+'/')

    # Copy and replace simulation files for the first window
    if int(win) == 0:
      if os.path.exists('amber_files'):
        shutil.rmtree('./amber_files')
      try:
        shutil.copytree('../../../amber_files', './amber_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      for dname, dirs, files in os.walk('./amber_files'):
        for fname in files:
          fpath = os.path.join(dname, fname)
          with open(fpath) as f:
            s = f.read()
            s = s.replace('_step_', dt).replace('_ntpr_', ntpr).replace('_ntwr_', ntwr).replace('_ntwe_', ntwe).replace('_ntwx_', ntwx).replace('_cutoff_', cut).replace('_gamma_ln_', gamma_ln).replace('_barostat_', barostat).replace('_receptor_ff_', receptor_ff).replace('_ligand_ff_', ligand_ff)
          with open(fpath, "w") as f:
            f.write(s)

      if os.path.exists('run_files'):
        shutil.rmtree('./run_files')
      try:
        shutil.copytree('../../../run_files', './run_files')
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

    # Create window directory
    if not os.path.exists('%s%02d' %(comp, int(win))):
      os.makedirs('%s%02d' %(comp, int(win)))
    os.chdir('%s%02d' %(comp, int(win)))
    if int(win) == 0:
      # Copy a few files
      for file in glob.glob('../../build_files/vac_ligand*'):
        shutil.copy(file, './')
      shutil.copy('../../build_files/%s.pdb' %mol.lower(), './')
      shutil.copy('../../build_files/fe-%s.pdb' %mol.lower(), './build-ini.pdb')
      shutil.copy('../../build_files/fe-%s.pdb' %mol.lower(), './')
      shutil.copy('../../build_files/anchors-'+pose+'.txt', './')
      for file in glob.glob('../../ff/%s.*' %mol.lower()):
        shutil.copy(file, './')
      for file in glob.glob('../../ff/dum.*'):
        shutil.copy(file, './')

      dum_coords = []
      recep_coords = []
      lig_coords = []
      dum_atomlist = []
      lig_atomlist = []
      recep_atomlist = []
      dum_atom = 0
      lig_atom = 0
      recep_atom = 0
      total_atom = 0
      resid_lig = 0
      resname_lig = mol
      resname_list = []
      resid_list = []

      # Read coordinates for dummy atoms
      if dec_method == 'sdr':
        for i in range(1, 3):
          shutil.copy('../../build_files/dum'+str(i)+'.pdb', './')
          with open('dum'+str(i)+'.pdb') as dum_in:
            lines = (line.rstrip() for line in dum_in)
            lines = list(line for line in lines if line)
            dum_coords.append((float(lines[1][30:38].strip()), float(lines[1][38:46].strip()), float(lines[1][46:54].strip())))
            dum_atomlist.append(lines[1][12:16].strip())
            resname_list.append(lines[1][17:20].strip())
            resid_list.append(float(lines[1][22:26].strip()))
            dum_atom += 1
            total_atom += 1
      else:
        for i in range(1, 2):
          shutil.copy('../../build_files/dum'+str(i)+'.pdb', './')
          with open('dum'+str(i)+'.pdb') as dum_in:
            lines = (line.rstrip() for line in dum_in)
            lines = list(line for line in lines if line)
            dum_coords.append((float(lines[1][30:38].strip()), float(lines[1][38:46].strip()), float(lines[1][46:54].strip())))
            dum_atomlist.append(lines[1][12:16].strip())
            resname_list.append(lines[1][17:20].strip())
            resid_list.append(float(lines[1][22:26].strip()))
            dum_atom += 1
            total_atom += 1


      # Read coordinates from aligned system
      with open('build-ini.pdb') as f_in:
        lines = (line.rstrip() for line in f_in)
        lines = list(line for line in lines if line) # Non-blank lines in a list   

      # Count atoms of receptor and ligand
      for i in range(0, len(lines)):
        if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
                if (lines[i][17:20].strip() != mol) and (lines[i][17:20].strip() != 'DUM'):
                   recep_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                   recep_atomlist.append(lines[i][12:16].strip())
                   resname_list.append(lines[i][17:20].strip())
                   resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                   recep_last = int(lines[i][22:26].strip())
                   recep_atom += 1
                   total_atom += 1
                elif lines[i][17:20].strip() == mol:
                   lig_coords.append((float(lines[i][30:38].strip()), float(lines[i][38:46].strip()), float(lines[i][46:54].strip())))
                   lig_atomlist.append(lines[i][12:16].strip())
                   resname_list.append(lines[i][17:20].strip())
                   resid_list.append(float(lines[i][22:26].strip()) + dum_atom)
                   lig_atom += 1
                   total_atom += 1


      coords = dum_coords + recep_coords + lig_coords
      atom_namelist = dum_atomlist + recep_atomlist + lig_atomlist
      lig_resid = recep_last + dum_atom + 1

      # Write the new pdb file

      build_file = open('build.pdb', 'w')

      # Positions for the dummy atoms 
      for i in range(0, dum_atom):
         build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
         build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))
         build_file.write('%6.2f%6.2f\n'%(0, 0))
         build_file.write('TER\n')

      # Positions of the receptor atoms
      for i in range(dum_atom , dum_atom + recep_atom):
          build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],resname_list[i], resid_list[i]))
          build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]), float(coords[i][2])))

          build_file.write('%6.2f%6.2f\n'%(0, 0))
      build_file.write('TER\n')

      # Positions of the ligand atoms
      for i in range(dum_atom + recep_atom, total_atom):
          build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, atom_namelist[i],mol, float(lig_resid)))
          build_file.write('%8.3f%8.3f%8.3f'%(float(coords[i][0]), float(coords[i][1]),float(coords[i][2])))

          build_file.write('%6.2f%6.2f\n'%(0, 0))

      build_file.write('TER\n')

      # Extra guests for decoupling

      build_file = open('build.pdb', 'a')
      if (comp == 'e'):
        for i in range(0, lig_atom):
            build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid+1)))
            build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2])))

            build_file.write('%6.2f%6.2f\n'%(0, 0))
        build_file.write('TER\n')
        if dec_method == 'sdr':
          for i in range(0, lig_atom):
              build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid+2)))
              build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2]+sdr_dist)))

              build_file.write('%6.2f%6.2f\n'%(0, 0))
          build_file.write('TER\n')
          for i in range(0, lig_atom):
              build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid+3)))
              build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2]+sdr_dist)))

              build_file.write('%6.2f%6.2f\n'%(0, 0))
        print('Creating new system for decharging...')
      if (comp == 'v' and dec_method == 'sdr'):
        for i in range(0, lig_atom):
            build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid + 1)))
            build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2]+sdr_dist)))

            build_file.write('%6.2f%6.2f\n'%(0, 0))
        print('Creating new system for vdw decoupling...')
      build_file.write('TER\n')
      build_file.write('END\n')
      build_file.close()

      if (comp == 'f' or comp == 'w'):
        # Create system with one or two ligands
        build_file = open('build.pdb', 'w')
        for i in range(0, lig_atom):
            build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid )))
            build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2])))

            build_file.write('%6.2f%6.2f\n'%(0, 0))
        build_file.write('TER\n')
        if comp == 'f':
          for i in range(0, lig_atom):
              build_file.write('%-4s  %5s %-4s %3s  %4.0f    '%('ATOM', i+1, lig_atomlist[i],mol, float(lig_resid + 1)))
              build_file.write('%8.3f%8.3f%8.3f'%(float(lig_coords[i][0]), float(lig_coords[i][1]),float(lig_coords[i][2])))

              build_file.write('%6.2f%6.2f\n'%(0, 0))
          build_file.write('TER\n')
        build_file.write('END\n')
        build_file.close()
        shutil.copy('./build.pdb', './%s.pdb' %mol.lower())
        tleap_vac = open('tleap_vac.in', 'w')
        tleap_vac.write('source leaprc.'+ligand_ff+'\n\n')
        tleap_vac.write('# Load the ligand parameters\n')
        tleap_vac.write('loadamberparams %s.frcmod\n'%(mol.lower()))
        tleap_vac.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), mol.lower()))
        tleap_vac.write('model = loadpdb %s.pdb\n\n' %(mol.lower()))
        tleap_vac.write('check model\n')
        tleap_vac.write('savepdb model vac.pdb\n')
        tleap_vac.write('saveamberparm model vac.prmtop vac.inpcrd\n')
        tleap_vac.write('quit\n\n')
        tleap_vac.close()

        p = sp.call('tleap -s -f tleap_vac.in > tleap_vac.log', shell=True)
        
    else:
      for file in glob.glob('../'+comp+'00/*'):
        shutil.copy(file, './')
    

def create_box(comp, hmr, pose, mol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, buffer_z, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, dec_method):
    
    # Copy and replace simulation files
    if stage != 'fe':
      if os.path.exists('amber_files'):
        shutil.rmtree('./amber_files')
      try:
        shutil.copytree('../amber_files', './amber_files')
      # Directories are the same
      except shutil.Error as e:
        print('Directory not copied. Error: %s' % e)
      # Any error saying that the directory doesn't exist
      except OSError as e:
        print('Directory not copied. Error: %s' % e)
      for dname, dirs, files in os.walk('./amber_files'):
        for fname in files:
          fpath = os.path.join(dname, fname)
          with open(fpath) as f:
            s = f.read()
            s = s.replace('_step_', dt).replace('_ntpr_', ntpr).replace('_ntwr_', ntwr).replace('_ntwe_', ntwe).replace('_ntwx_', ntwx).replace('_cutoff_', cut).replace('_gamma_ln_', gamma_ln).replace('_barostat_', barostat).replace('_receptor_ff_', receptor_ff).replace('_ligand_ff_', ligand_ff)
          with open(fpath, "w") as f:
            f.write(s)
      os.chdir(pose)

    # Copy tleap files that are used for restraint generation and analysis
    shutil.copy('../amber_files/tleap.in.amber16', 'tleap_vac.in')
    shutil.copy('../amber_files/tleap.in.amber16', 'tleap_vac_ligand.in')
    shutil.copy('../amber_files/tleap.in.amber16', 'tleap.in')

    # Copy ligand parameter files
    for file in glob.glob('../ff/*'):
       shutil.copy(file, './')

    # Append tleap file for vacuum
    tleap_vac = open('tleap_vac.in', 'a')
    tleap_vac.write('# Load the ligand parameters\n')        
    tleap_vac.write('loadamberparams %s.frcmod\n'%(mol.lower()))
    tleap_vac.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), mol.lower()))
    tleap_vac.write('model = loadpdb build.pdb\n\n')
    tleap_vac.write('check model\n')
    tleap_vac.write('savepdb model vac.pdb\n')
    tleap_vac.write('saveamberparm model vac.prmtop vac.inpcrd\n')
    tleap_vac.write('quit\n')
    tleap_vac.close()

    # Append tleap file for ligand only
    tleap_vac_ligand = open('tleap_vac_ligand.in', 'a')
    tleap_vac_ligand.write('# Load the ligand parameters\n')        
    tleap_vac_ligand.write('loadamberparams %s.frcmod\n'%(mol.lower()))
    tleap_vac_ligand.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), mol.lower()))
    tleap_vac_ligand.write('model = loadpdb %s.pdb\n\n' %(mol.lower()))
    tleap_vac_ligand.write('check model\n')
    tleap_vac_ligand.write('savepdb model vac_ligand.pdb\n')
    tleap_vac_ligand.write('saveamberparm model vac_ligand.prmtop vac_ligand.inpcrd\n')
    tleap_vac_ligand.write('quit\n')
    tleap_vac_ligand.close()


    # Generate complex in vacuum
    p = sp.call('tleap -s -f tleap_vac.in > tleap_vac.log', shell=True)

    # Generate ligand structure in vacuum
    p = sp.call('tleap -s -f tleap_vac_ligand.in > tleap_vac_ligand.log', shell=True)

    # Find out how many cations/anions are needed for neutralization
    neu_cat = 0
    neu_ani = 0
    f = open('tleap_vac.log', 'r')
    for line in f:
        if "The unperturbed charge of the unit" in line:
            splitline = line.split()
            if float(splitline[6].strip('\'\",.:;#()][')) < 0:
                neu_cat = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
            elif float(splitline[6].strip('\'\",.:;#()][')) > 0:
                neu_ani = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
    f.close()
     
    # Get ligand removed charge when doing LJ calculations
    lig_cat = 0
    lig_ani = 0
    f = open('tleap_vac_ligand.log', 'r')
    for line in f:
        if "The unperturbed charge of the unit" in line:
            splitline = line.split()
            if float(splitline[6].strip('\'\",.:;#()][')) < 0:
                lig_cat = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
            elif float(splitline[6].strip('\'\",.:;#()][')) > 0:
                lig_ani = round(float(re.sub('[+-]', '', splitline[6].strip('\'\"-,.:;#()]['))))
    f.close()
     
    # Adjust ions for LJ and electrostatic Calculations (avoid neutralizing plasma)
    if comp == 'v' and dec_method == 'sdr':
      charge_neut = neu_cat - neu_ani - 2*lig_cat + 2*lig_ani
      neu_cat = 0
      neu_ani = 0
      if charge_neut > 0:
        neu_cat = abs(charge_neut)
      if charge_neut < 0:
        neu_ani = abs(charge_neut)
    if comp == 'e' and dec_method == 'sdr':
      charge_neut = neu_cat - neu_ani - 3*lig_cat + 3*lig_ani
      neu_cat = 0
      neu_ani = 0
      if charge_neut > 0:
        neu_cat = abs(charge_neut)
      if charge_neut < 0:
        neu_ani = abs(charge_neut)
    
    # Define volume density for different water models
    if water_model == 'TIP3P':
       water_box = water_model.upper()+'BOX'
       ratio = 0.0576
    elif water_model == 'SPCE': 
       water_box = 'SPCBOX'
       ratio = 0.0576
    elif water_model == 'TIP4PEW': 
       water_box = water_model.upper()+'BOX'
       ratio = 0.0573

    # Fixed number of water molecules
    if num_waters != 0:

      # Create the first box guess to get the initial number of waters and cross sectional area
      buff = 50.0  
      scripts.write_tleap(mol, water_model, water_box, buff, buffer_x, buffer_y)
      num_added = scripts.check_tleap()
      cross_area = scripts.cross_sectional_area()

      # First iteration to estimate box volume and number of ions
      res_diff = num_added - num_waters 
      buff_diff = res_diff/(ratio*cross_area)
      buff -= buff_diff 
      print(buff)
      if buff < 0:
        print ('Not enough water molecules to fill the system in the z direction, please increase the number of water molecules')
        sys.exit(1)
      # Get box volume and number of added ions
      scripts.write_tleap(mol, water_model, water_box, buff, buffer_x, buffer_y)
      box_volume = scripts.box_volume()
      print(box_volume)
      num_cations = round(0.85*ion_def[2]*6.02e23*box_volume*1e-27) # 0.85 factor to account for some shrinking of the box during equilibration
      print(num_cations)

      # Update target number of residues according to the ion definitions 
      if (neut == 'no'):
        target_num = int(num_waters - neu_cat + neu_ani + 2*int(num_cations)) 
      elif (neut == 'yes'):
        target_num = int(num_waters + neu_cat + neu_ani)
    
      # Define a few parameters for solvation iteration
      buff = 50.0  
      count = 0
      max_count = 10
      rem_limit = 16
      factor = 1
      ind = 0.90   
      buff_diff = 1.0  

      # Iterate to get the correct number of waters
      while num_added != target_num:
          count += 1
          if count > max_count:
          # Try different parameters
               rem_limit += 4
               if ind > 0.5:
                 ind = ind - 0.02
               else:
                 ind = 0.90
               factor = 1
               max_count = max_count + 10
          tleap_remove = None
          # Manually remove waters if inside removal limit
          if num_added > target_num and (num_added - target_num) < rem_limit:
              difference = num_added - target_num
              tleap_remove = [target_num + 1 + i for i in range(difference)]
              scripts.write_tleap(mol, water_model, water_box, buff, buffer_x, buffer_y, tleap_remove)
              scripts.check_tleap()
              break
          # Set new buffer size based on chosen water density
          res_diff = num_added - target_num - (rem_limit/2)
          buff_diff = res_diff/(ratio*cross_area)
          buff -= (buff_diff * factor)
          if buff < 0:
             print ('Not enough water molecules to fill the system in the z direction, please increase the number of water molecules')
             sys.exit(1)
          # Set relaxation factor  
          factor = ind * factor
          # Get number of waters
          scripts.write_tleap(mol, water_model, water_box, buff, buffer_x, buffer_y)
          num_added = scripts.check_tleap()       
    # Fixed z buffer 
    elif buffer_z != 0:
      buff = buffer_z
      tleap_remove = None
      # Get box volume and number of added ions
      scripts.write_tleap(mol, water_model, water_box, buff, buffer_x, buffer_y)
      box_volume = scripts.box_volume()
      print(box_volume)
      num_cations = round(0.85*ion_def[2]*6.02e23*box_volume*1e-27) # 0.85 factor to account for some shrinking of the box during equilibration
      print(num_cations)
 
    # Number of cations and anions   
    num_cat = num_cations
    num_ani = num_cations - neu_cat + neu_ani
    
    # Write the final tleap file with the correct system size and removed water molecules
    shutil.copy('tleap.in', 'tleap_solvate.in')
    tleap_solvate = open('tleap_solvate.in', 'a')
    tleap_solvate.write('# Load the ligand parameters\n')        
    tleap_solvate.write('loadamberparams %s.frcmod\n'%(mol.lower()))
    tleap_solvate.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), mol.lower()))
    tleap_solvate.write('model = loadpdb build.pdb\n\n')
    tleap_solvate.write('# Load the water and jc ion parameters\n')        
    tleap_solvate.write('source leaprc.water.%s\n'%(water_model.lower()))
    tleap_solvate.write('loadamberparams frcmod.ionsjc_%s\n\n'%(water_model.lower()))
    tleap_solvate.write('# Create water box with chosen model\n')
    tleap_solvate.write('solvatebox model ' + water_box + ' {'+ str(buffer_x) +' '+ str(buffer_y) +' '+ str(buff) +'}\n\n')
    if tleap_remove is not None:
        tleap_solvate.write('# Remove a few waters manually\n')
        for water in tleap_remove:
            tleap_solvate.write('remove model model.%s\n' % water)
        tleap_solvate.write('\n')
    # Ionize/neutralize system
    if (neut == 'no'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[0], num_cat))
        tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[1], num_ani))
    elif (neut == 'yes'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        if neu_cat != 0:
          tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[0], neu_cat))
        if neu_ani != 0:
          tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[1], neu_ani))
    tleap_solvate.write('\n')
    tleap_solvate.write('desc model\n')
    tleap_solvate.write('savepdb model full.pdb\n')
    tleap_solvate.write('saveamberparm model full.prmtop full.inpcrd\n')
    tleap_solvate.write('quit')
    tleap_solvate.close()
    p = sp.call('tleap -s -f tleap_solvate.in > tleap_solvate.log', shell=True)
    
    f = open('tleap_solvate.log', 'r')
    for line in f:
        if "Could not open file" in line:
           print ('WARNING!!!')
           print (line)
           sys.exit(1)
        if "WARNING: The unperturbed charge of the unit:" in line:
           print (line)
           print ('The system is not neutralized properly after solvation')
        if "addIonsRand: Argument #2 is type String must be of type: [unit]" in line:
           print('Aborted.The ion types specified in the input file could be wrong.')              
           print('Please check the tleap_solvate.log file, and the ion types specified in the input file.\n')
           sys.exit(1)
    f.close()


    # Apply hydrogen mass repartitioning
    print('Applying mass repartitioning...')
    shutil.copy('../amber_files/parmed-hmr.in', './')
    sp.call('parmed -O -n -i parmed-hmr.in > parmed-hmr.log', shell=True)

    if stage != 'fe':
      os.chdir('../')
   
def ligand_box(mol, lig_buffer, water_model, neut, ion_def, comp, ligand_ff):
    # Define volume density for different water models
    if water_model == 'TIP3P':
       water_box = water_model.upper()+'BOX'
    elif water_model == 'SPCE': 
       water_box = 'SPCBOX'
    elif water_model == 'TIP4PEW': 
       water_box = water_model.upper()+'BOX'

    # Copy ligand parameter files
    for file in glob.glob('../../ff/%s.*' %mol.lower()):
        shutil.copy(file, './')

    # Write and run preliminary tleap file
    tleap_solvate = open('tmp_tleap.in', 'w')
    tleap_solvate.write('source leaprc.'+ligand_ff+'\n\n')        
    tleap_solvate.write('# Load the ligand parameters\n')        
    tleap_solvate.write('loadamberparams %s.frcmod\n'%(mol.lower()))
    tleap_solvate.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), mol.lower()))
    tleap_solvate.write('model = loadpdb %s.pdb\n\n' %(mol.lower()))
    tleap_solvate.write('# Load the water and jc ion parameters\n')        
    tleap_solvate.write('source leaprc.water.%s\n'%(water_model.lower()))
    tleap_solvate.write('loadamberparams frcmod.ionsjc_%s\n\n'%(water_model.lower()))
    tleap_solvate.write('check model\n')
    tleap_solvate.write('savepdb model vac.pdb\n')
    tleap_solvate.write('saveamberparm model vac.prmtop vac.inpcrd\n\n')
    tleap_solvate.write('# Create water box with chosen model\n')
    tleap_solvate.write('solvatebox model ' + water_box + ' '+str(lig_buffer)+'\n\n')
    tleap_solvate.write('quit\n')
    tleap_solvate.close()

    # Get box volume and number of added ions
    box_volume = scripts.box_volume()
    print(box_volume)
    num_cations = round(0.85*ion_def[2]*6.02e23*box_volume*1e-27) # 0.85 factor to account for some shrinking of the box during equilibration
    print(num_cations)


     # Write and run tleap file
    tleap_solvate = open('tleap_solvate.in', 'a')
    tleap_solvate.write('source leaprc.'+ligand_ff+'\n\n')        
    tleap_solvate.write('# Load the ligand parameters\n')        
    tleap_solvate.write('loadamberparams %s.frcmod\n'%(mol.lower()))
    tleap_solvate.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), mol.lower()))
    tleap_solvate.write('model = loadpdb %s.pdb\n\n' %(mol.lower()))
    tleap_solvate.write('# Load the water and jc ion parameters\n')        
    tleap_solvate.write('source leaprc.water.%s\n'%(water_model.lower()))
    tleap_solvate.write('loadamberparams frcmod.ionsjc_%s\n\n'%(water_model.lower()))
    tleap_solvate.write('check model\n')
    tleap_solvate.write('savepdb model vac.pdb\n')
    tleap_solvate.write('saveamberparm model vac.prmtop vac.inpcrd\n\n')
    tleap_solvate.write('# Create water box with chosen model\n')
    tleap_solvate.write('solvatebox model ' + water_box + ' '+str(lig_buffer)+'\n\n')
    if (neut == 'no'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        tleap_solvate.write('addionsrand model %s %d\n' % (ion_def[0], num_cations))
        tleap_solvate.write('addionsrand model %s 0\n' % (ion_def[1]))
    elif (neut == 'yes'):
        tleap_solvate.write('# Add ions for neutralization/ionization\n')
        tleap_solvate.write('addionsrand model %s 0\n' % (ion_def[0]))
        tleap_solvate.write('addionsrand model %s 0\n' % (ion_def[1]))
    tleap_solvate.write('\n')
    tleap_solvate.write('desc model\n')
    tleap_solvate.write('savepdb model full.pdb\n')
    tleap_solvate.write('saveamberparm model full.prmtop full.inpcrd\n')
    tleap_solvate.write('quit\n')
    tleap_solvate.close()
    p = sp.call('tleap -s -f tleap_solvate.in > tleap_solvate.log', shell=True)
 
    # Apply hydrogen mass repartitioning
    print('Applying mass repartitioning...')
    shutil.copy('../amber_files/parmed-hmr.in', './')
    sp.call('parmed -O -n -i parmed-hmr.in > parmed-hmr.log', shell=True)

    # Copy a few files for consistency
    if (comp != 'f' and comp != 'w'):
      shutil.copy('./vac.pdb','./vac_ligand.pdb')
      shutil.copy('./vac.prmtop','./vac_ligand.prmtop')
