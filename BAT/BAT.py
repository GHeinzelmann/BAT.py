#!/usr/bin/env python2
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
from lib import build 
from lib import scripts 
from lib import setup
from lib import analysis

ion_def = []
poses_list = []
poses_def = []
release_eq = []
translate_apr = []
attach_rest = []
lambdas = []  
weights = []  
components = []  
aa1_poses = []  
aa2_poses = []  

# Predefine some variables to avoid errors

a_steps1 = 0
a_steps2 = 0
l_steps1 = 0
l_steps2 = 0
t_steps1 = 0
t_steps2 = 0
c_steps1 = 0
c_steps2 = 0
r_steps1 = 0
r_steps2 = 0
e_steps1 = 0
e_steps2 = 0
v_steps1 = 0
v_steps2 = 0
f_steps1 = 0
f_steps2 = 0
w_steps1 = 0
w_steps2 = 0
u_steps1 = 0
u_steps2 = 0
sdr_dist = 0

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

        if lines[i][0] == 'pull_ligand':
            if lines[i][1].lower() == 'yes':
                pull_ligand = 'yes'
            elif lines[i][1].lower() == 'no':
                pull_ligand = 'no'
            else:
                print('Wrong input! Please use yes or no to indicate whether to pull out the ligand or not.')
                sys.exit(1)

        elif lines[i][0] == 'temperature':
            temperature = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'eq_steps1':
            eq_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'eq_steps2':
            eq_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'prep_steps1':
            prep_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'prep_steps2':
            prep_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'a_steps1':
            a_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'a_steps2':
            a_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l_steps1':
            l_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l_steps2':
            l_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 't_steps1':
            t_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 't_steps2':
            t_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'u_steps1':
            u_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'u_steps2':
            u_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'c_steps1':
            c_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'c_steps2':
            c_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'r_steps1':
            r_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'r_steps2':
            r_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'e_steps1':
            e_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'e_steps2':
            e_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'v_steps1':
            v_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'v_steps2':
            v_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'f_steps1':
            f_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'f_steps2':
            f_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'w_steps1':
            w_steps1 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'w_steps2':
            w_steps2 = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'pull_spacing':
            pull_spacing = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'poses_list':
            newline = lines[i][1].strip('\'\"-,.:;#()][').split(',')
            for j in range(0, len(newline)):
                poses_list.append(scripts.check_input('int', newline[j], input_file, lines[i][0]))
        elif lines[i][0] == 'calc_type':
            calc_type = lines[i][1].lower()
        elif lines[i][0] == 'celpp_receptor':
            celp_st = lines[i][1]
        elif lines[i][0] == 'p1':
            H1 = lines[i][1]
        elif lines[i][0] == 'p2':
            H2 = lines[i][1]
        elif lines[i][0] == 'p3':
            H3 = lines[i][1]
        elif lines[i][0] == 'ligand_name':
            mol = lines[i][1]
        elif lines[i][0] == 'fe_type':
            if lines[i][1].lower() == 'rest':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'dd':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'sdr':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'pmf':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'pmf-rest':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'sdr-rest':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'dd-rest':
                fe_type = lines[i][1].lower()
            elif lines[i][1].lower() == 'custom':
                fe_type = lines[i][1].lower()
            else:
                print('Free energy type not recognized, please choose rest (restraints only), dd (double decoupling), sdr (simultaneous decoupling-recoupling), pmf (umbrella sampling), dd-rest, sdr-rest, pmf-rest, or custom')
                sys.exit(1)
        elif lines[i][0] == 'dec_int':
            if lines[i][1].lower() == 'mbar':
                dec_int = lines[i][1].lower()
            elif lines[i][1].lower() == 'ti':
                dec_int = lines[i][1].lower()
            else:
                print('Decoupling integration method not recognized, please choose ti or mbar')
                sys.exit(1)
        elif lines[i][0] == 'dec_method':
            if lines[i][1].lower() == 'dd':
                dec_method = lines[i][1].lower()
            elif lines[i][1].lower() == 'sdr':
                dec_method = lines[i][1].lower()
            else:
                print('Decoupling method not recognized, please choose dd or sdr')
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
        elif lines[i][0] == 'water_model':
            if lines[i][1].lower() == 'tip3p':
                water_model = lines[i][1].upper()
            elif lines[i][1].lower() == 'tip4pew':
                water_model = lines[i][1].upper()
            elif lines[i][1].lower() == 'spce':
                water_model = lines[i][1].upper()
            else:
                print('Water model not supported. Please choose TIP3P, TIP4PEW or SPCE')
                sys.exit(1)
        elif lines[i][0] == 'num_waters':
            num_waters = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'neutralize_only':
            if lines[i][1].lower() == 'yes':
                neut = 'yes'
            elif lines[i][1].lower() == 'no':
                neut = 'no'
            else:
                print('Wrong input! Please choose neutralization only or add extra ions')
                sys.exit(1)
        elif lines[i][0] == 'cation':
            cation = lines[i][1]
        elif lines[i][0] == 'anion':
            anion = lines[i][1]
        elif lines[i][0] == 'num_cations':
            num_cations = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'num_cat_ligbox':
            num_cat_ligbox = scripts.check_input('int', lines[i][1], input_file, lines[i][0])
        elif lines[i][0] == 'buffer_x':
            buffer_x = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'buffer_y':
            buffer_y = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_buffer':
            lig_buffer = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_distance_force':
            rec_distance_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_angle_force':
            rec_angle_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_dihcf_force':
            rec_dihcf_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_discf_force':
            rec_discf_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_distance_force':
            lig_distance_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_angle_force':
            lig_angle_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_dihcf_force':
            lig_dihcf_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'lig_discf_force':
            lig_discf_force = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'sdr_dist':
            sdr_dist = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l1_x':
            l1_x = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l1_y':
            l1_y = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l1_z':
            l1_z = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l1_zm':
            l1_zm = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'l1_range':
            l1_range = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'min_adis':
            min_adis = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'max_adis':
            max_adis = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'rec_bb':
            if lines[i][1].lower() == 'yes':
                rec_bb = 'yes'
            elif lines[i][1].lower() == 'no':
                rec_bb = 'no'
            else:
                print('Wrong input! Please use yes or no to indicate whether protein backbone restraints'
                      'will be used.')
                sys.exit(1)
        elif lines[i][0] == 'bb_start':
            bb_start = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'bb_end':
            bb_end = scripts.check_input('int', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'bb_equil':
            if lines[i][1].lower() == 'yes':
                bb_equil = lines[i][1].lower()
            else:
                bb_equil = 'no'
        elif lines[i][0] == 'release_eq':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                release_eq.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'translate_apr':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                translate_apr.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'attach_rest':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                attach_rest.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'lambdas':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                lambdas.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'weights':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                weights.append(scripts.check_input('float', strip_line[j], input_file, lines[i][0]))
        elif lines[i][0] == 'components':
            strip_line = lines[i][1].strip('\'\"-,.:;#()][').split()
            for j in range(0, len(strip_line)):
                components.append(strip_line[j])
        elif lines[i][0] == 'ntpr':
            ntpr = lines[i][1]
        elif lines[i][0] == 'ntwr':
            ntwr = lines[i][1]
        elif lines[i][0] == 'ntwe':
            ntwe = lines[i][1]
        elif lines[i][0] == 'ntwx':
            ntwx = lines[i][1]
        elif lines[i][0] == 'cut':
            cut = lines[i][1]
        elif lines[i][0] == 'gamma_ln':
            gamma_ln = lines[i][1]
        elif lines[i][0] == 'barostat':
           barostat = lines[i][1]
        elif lines[i][0] == 'receptor_ff':
            receptor_ff = lines[i][1]
        elif lines[i][0] == 'ligand_ff':
            if lines[i][1].lower() == 'gaff':
                ligand_ff = 'gaff'
            elif lines[i][1].lower() == 'gaff2':
                ligand_ff = 'gaff2'
            else:
                print('Wrong input! Available options for ligand force-field are gaff and gaff2')
                sys.exit(1)
        elif lines[i][0] == 'ligand_ph':
            ligand_ph = scripts.check_input('float', lines[i][1], input_file, lines[i][0]) 
        elif lines[i][0] == 'dt':
            dt = lines[i][1]



# Number of simulations, 1 equilibrium and 1 production
apr_sim = 2

# Define free energy components
if fe_type == 'rest':
  components = ['c', 'a', 'l', 't', 'r'] 
  dec_method = 'dd'
elif fe_type == 'sdr':
  components = ['e', 'v'] 
  dec_method = 'sdr'
elif fe_type == 'dd':
  components = ['e', 'v', 'f', 'w'] 
  dec_method = 'dd'
elif fe_type == 'pmf':
  components = ['u'] 
  dec_method = 'dd'
elif fe_type == 'pmf-rest':
  components = ['c', 'a', 'l', 't', 'r', 'u'] 
  dec_method = 'dd'
elif fe_type == 'sdr-rest':
  components = ['c', 'a', 'l', 't', 'r', 'e', 'v'] 
  dec_method = 'sdr'
elif fe_type == 'dd-rest':
  components = ['c', 'a', 'l', 't', 'r', 'e', 'v', 'f', 'w'] 
  dec_method = 'dd'

# Pull ligand out or not
if pull_ligand == 'no':
  translate_apr = [ 0.00 ]
  pull_spacing = 1.0
  prep_steps2 = 0  

# Do not apply protein backbone restraints
if rec_bb == 'no':
  bb_start =  1
  bb_end   =  0
  bb_equil = 'no'

# Create poses definitions
if calc_type == 'dock':
  for i in range(0, len(poses_list)):
    poses_def.append('pose'+str(poses_list[i]))
elif calc_type == 'crystal':
  poses_def = [celp_st]

# Total distance
apr_distance = translate_apr[-1]
rng = 0

# Create restraint definitions
rest = [rec_distance_force, rec_angle_force, rec_dihcf_force, rec_discf_force, lig_distance_force, lig_angle_force, lig_dihcf_force, lig_discf_force]

# Create ion definitions
ion_def = [cation, anion, num_cations]
ion_lig = [cation, anion, num_cat_ligbox]

# Define number of steps for all stages
dic_steps1 = {}
dic_steps2 = {}
dic_steps1['a'] = a_steps1
dic_steps2['a'] = a_steps2
dic_steps1['l'] = l_steps1
dic_steps2['l'] = l_steps2
dic_steps1['t'] = t_steps1
dic_steps2['t'] = t_steps2
dic_steps1['c'] = c_steps1
dic_steps2['c'] = c_steps2
dic_steps1['r'] = r_steps1
dic_steps2['r'] = r_steps2
dic_steps1['v'] = v_steps1
dic_steps2['v'] = v_steps2
dic_steps1['e'] = e_steps1
dic_steps2['e'] = e_steps2
dic_steps1['w'] = w_steps1
dic_steps2['w'] = w_steps2
dic_steps1['f'] = f_steps1
dic_steps2['f'] = f_steps2

if stage == 'equil':
  comp = 'q'
  win = 0
  trans_dist = 0
  # Create equilibrium systems for all poses listed in the input file
  for i in range(0, len(poses_def)):
    rng = len(release_eq) - 1
    pose = poses_def[i]
    if not os.path.exists('./all-poses/'+pose+'.pdb'):
      continue
    print('Setting up '+str(poses_def[i]))
    # Get number of simulations
    num_sim = len(release_eq)
    # Create aligned initial complex
    anch = build.build_equil(pose, celp_st, mol, H1, H2, H3, calc_type, l1_x, l1_y, l1_z, l1_zm, l1_range, min_adis, max_adis, ligand_ff, ligand_ph)
    if anch == 'anch1':
      aa1_poses.append(pose)
      os.chdir('../')
      continue
    if anch == 'anch2':
      aa2_poses.append(pose)
      os.chdir('../')
      continue
    # Solvate system with ions
    print('Creating box...')
    build.create_box(comp, hmr, pose, mol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, dec_method)
    # Apply restraints and prepare simulation files
    print('Equil release weights:')
    for i in range(0, len(release_eq)):
      weight = release_eq[i]
      print('%s' %str(weight))
      setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
      shutil.copy('./'+pose+'/disang.rest', './'+pose+'/disang%02d.rest' %int(i))
    shutil.copy('./'+pose+'/disang%02d.rest' %int(0), './'+pose+'/disang.rest')
    setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, eq_steps1, eq_steps2, rng)
    os.chdir('../')
  if len(aa1_poses) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand first anchor L1 for', aa1_poses)
    print ('The ligand is most likely not in the defined binding site in these systems.')
  if len(aa2_poses) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand L2 or L3 anchors for', aa2_poses)
    print ('Try reducing the min_adis parameter in the input file.')

elif stage == 'prep':
  win = 0
  weight = 100.0
  comp = 's'
  # Prepare systems after equilibration for poses listed in the input file
  for i in range(0, len(poses_def)):
    pose = poses_def[i]
    if not os.path.exists('./equil/'+pose):
      continue
    print('Setting up '+str(poses_def[i]))
    # Get number of simulations
    num_sim = int(apr_distance/pull_spacing)+1
    rng = num_sim - 1
    # Create aligned initial complex
    fwin = len(release_eq) - 1
    anch = build.build_prep(pose, mol, fwin, l1_x, l1_y, l1_z, l1_zm, l1_range, min_adis, max_adis)
    if anch == 'anch1':
      aa1_poses.append(pose)
      os.chdir('../')
      continue
    if anch == 'anch2':
      aa2_poses.append(pose)
      os.chdir('../')
      continue
    # Solvate system with ions
    print('Creating box...')
    build.create_box(comp, hmr, pose, mol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, dec_method)
    # Apply restraints and prepare simulation files
    print('Creating preparation steps...')
    for i in range(0, num_sim):
      trans_dist = float(i*pull_spacing)
      setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
      shutil.copy('./'+pose+'/disang.rest', './'+pose+'/disang%03d.rest' %int(i))
    shutil.copy('./'+pose+'/disang%03d.rest' %int(0), './'+pose+'/disang.rest')
    setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, prep_steps1, prep_steps2, rng)
    os.chdir('../')
  if len(aa1_poses) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand first anchor L1 for', aa1_poses)
    print ('The ligand most likely left the binding site during equilibration.')
  if len(aa2_poses) != 0:
    print('\n')
    print ('WARNING: Could not find the ligand L2 or L3 anchors for', aa2_poses)
    print ('Try reducing the min_adis parameter in the input file.')
elif stage == 'fe':
  # Create systems for all poses after preparation
  num_sim = apr_sim
  # Create and move to free energy directory
  if not os.path.exists('fe'):
    os.makedirs('fe')
  os.chdir('fe')
  for i in range(0, len(poses_def)):
    pose = poses_def[i]
    if not os.path.exists('../prep/'+pose):
      continue
    print('Setting up '+str(poses_def[i]))
    # Create and move to pose directory
    if not os.path.exists(pose):
      os.makedirs(pose)
    os.chdir(pose)
    # Generate folder and restraints for all components and windows
    for j in range(0, len(components)):
      comp = components[j]
      # Translation (umbrella)
      if (comp == 'u'):
        if not os.path.exists('pmf'):
          os.makedirs('pmf')
        os.chdir('pmf')
        weight = 100.0
        for k in range(0, len(translate_apr)):
          trans_dist = translate_apr[k]
          win = k
          print('window: %s%02d distance: %s' %(comp, int(win), str(trans_dist)))
          build.build_apr(hmr, mol, pose, comp, win, trans_dist, pull_spacing, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt)
          setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, u_steps1, u_steps2, rng)
        os.chdir('../')  
      # Ligand conformational release in a small box
      elif (comp == 'c'):          
        if not os.path.exists('rest'):
          os.makedirs('rest')
        os.chdir('rest')
        trans_dist = 0
        for k in range(0, len(attach_rest)):
          weight = attach_rest[k]
          win = k
          if int(win) == 0:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_apr(hmr, mol, pose, comp, win, trans_dist, pull_spacing, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt)
            print('Creating box for ligand only...')
            build.ligand_box(mol, lig_buffer, water_model, neut, ion_lig, comp, ligand_ff)
            setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
            setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, c_steps1, c_steps2, rng)
          else:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_apr(hmr, mol, pose, comp, win, trans_dist, pull_spacing, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt)
            setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
            setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, c_steps1, c_steps2, rng)
        os.chdir('../')  
      # Receptor conformational release in a separate box
      elif (comp == 'r'):          
        if not os.path.exists('rest'):
          os.makedirs('rest')
        os.chdir('rest')
        trans_dist = 0
        for k in range(0, len(attach_rest)):
          weight = attach_rest[k]
          win = k
          if int(win) == 0:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_apr(hmr, mol, pose, comp, win, trans_dist, pull_spacing, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt)
            print('Creating box for apo protein...')
            build.create_box(comp, hmr, pose, mol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, dec_method)
            setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
            setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, r_steps1, r_steps2, rng)
          else:
            print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
            build.build_apr(hmr, mol, pose, comp, win, trans_dist, pull_spacing, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt)
            setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
            setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, r_steps1, r_steps2, rng)
        os.chdir('../')  
      # Simultaneous/double decoupling
      elif (comp == 'v' or comp == 'e'):          
        steps1 = dic_steps1[comp]
        steps2 = dic_steps2[comp]
        if not os.path.exists(dec_method):
          os.makedirs(dec_method)
        os.chdir(dec_method)
        trans_dist = 0
        for k in range(0, len(lambdas)):
          weight = lambdas[k]
          win = k
          print('window: %s%02d lambda: %s' %(comp, int(win), str(weight)))
          if int(win) == 0:
            build.build_dec(hmr, mol, pose, comp, win, water_model, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, sdr_dist, dec_method)
            build.create_box(comp, hmr, pose, mol, num_waters, water_model, ion_def, neut, buffer_x, buffer_y, stage, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, dec_method)
            setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
            setup.dec_files(temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, weight, lambdas, dec_method)
          else:
            build.build_dec(hmr, mol, pose, comp, win, water_model, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, sdr_dist, dec_method)
            setup.dec_files(temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, weight, lambdas, dec_method)
        os.chdir('../')  
      # Bulk systems for dd
      elif (comp == 'f' or comp == 'w'):
        steps1 = dic_steps1[comp]
        steps2 = dic_steps2[comp]
        if not os.path.exists('dd'):
          os.makedirs('dd')
        os.chdir('dd')
        trans_dist = 0
        for k in range(0, len(lambdas)):
          weight = lambdas[k]
          win = k
          if int(win) == 0:
            print('window: %s%02d lambda: %s' %(comp, int(win), str(weight)))
            build.build_dec(hmr, mol, pose, comp, win, water_model, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, sdr_dist, dec_method)
            print('Creating box for ligand decharging in bulk...')
            build.ligand_box(mol, lig_buffer, water_model, neut, ion_lig, comp, ligand_ff)
            setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
            setup.dec_files(temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, weight, lambdas, dec_method)
          else:
            print('window: %s%02d lambda: %s' %(comp, int(win), str(weight)))
            build.build_dec(hmr, mol, pose, comp, win, water_model, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt, sdr_dist, dec_method)
            setup.dec_files(temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, weight, lambdas, dec_method)
        os.chdir('../')
      # Attachments in the bound system
      else:          
        if not os.path.exists('rest'):
          os.makedirs('rest')
        os.chdir('rest')
        trans_dist = 0
        for k in range(0, len(attach_rest)):
          weight = attach_rest[k]
          win = k
          print('window: %s%02d weight: %s' %(comp, int(win), str(weight)))
          build.build_apr(hmr, mol, pose, comp, win, trans_dist, pull_spacing, ntpr, ntwr, ntwe, ntwx, cut, gamma_ln, barostat, receptor_ff, ligand_ff, dt)
          setup.restraints(pose, rest, bb_start, bb_end, weight, stage, mol, trans_dist, comp, bb_equil, sdr_dist, dec_method)
          steps1 = dic_steps1[comp]
          steps2 = dic_steps2[comp]
          setup.sim_files(hmr, temperature, mol, num_sim, pose, comp, win, stage, steps1, steps2, rng)
        os.chdir('../')  
    os.chdir('../')
elif stage == 'analysis':
  # Free energies MBAR/TI and analytical calculations
  for i in range(0, len(poses_def)):
    pose = poses_def[i]
    analysis.fe_values(blocks, components, temperature, pose, attach_rest, translate_apr, lambdas, weights, dec_int, dec_method, rest)
    os.chdir('../../')
