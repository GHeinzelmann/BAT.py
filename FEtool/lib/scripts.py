#!/usr/bin/env python2
import datetime as dt
import glob as glob
import os as os
import re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys

def help_message():
    print('Use the flags -i and -s for the input file and current stage of the calculations')
    print('Example: python APRprotein.py -i input.in -s equil')

def write_tleap(mol, water_model, water_box, buff, buffer_x, buffer_y, tleap_remove=None):
    shutil.copy('tleap.in', 'tmp_tleap.in')
    tmp_file = open('tmp_tleap.in', 'a')
    tmp_file.write('# Load the ligand parameters\n')        
    tmp_file.write('loadamberparams %s.frcmod\n'%(mol.lower()))
    tmp_file.write('%s = loadmol2 %s.mol2\n\n'%(mol.upper(), mol.lower()))
    tmp_file.write('model = loadpdb build.pdb\n\n')
    tmp_file.write('# Load the water and jc ion parameters\n')        
    tmp_file.write('source leaprc.water.%s\n'%(water_model.lower()))
    tmp_file.write('loadamberparams frcmod.ionsjc_%s\n\n'%(water_model.lower()))
    tmp_file.write('solvatebox model ' + water_box + ' {'+ str(buffer_x) +' '+ str(buffer_y) +' '+ str(buff) + '}\n')
    if tleap_remove is not None:
        for water in tleap_remove:
            tmp_file.write('remove model model.%s\n' % water)
        tmp_file.write('\n')
    tmp_file.write('quit')
    tmp_file.close()


def check_tleap():
    p = sp.call('tleap -s -f tmp_tleap.in > tmp.log', shell=True)
    # Get how many residues were added to the system
    num_added = None
    f = open('tmp.log', 'r')
    for line in f:
        if "Added" in line and "residues" in line:
            num_added = line[line.find('Added') + len('Added'):].split()
            num_added = int(num_added[0])
    f.close()
    if num_added is None:
  	print('Problem solvating system, check your input file and the tleap log file\n')
        sys.exit(1)
    return num_added


def cross_sectional_area():
    p = sp.call('tleap -s -f tmp_tleap.in > tmp.log', shell=True)
    # Get the total box size in the three axes and the xy cross sectional area
    num_added = None
    f = open('tmp.log', 'r')
    for line in f:
        if "Total vdw box size" in line:
            splitdata = line.split()
            x_axis = float(splitdata[4])
            y_axis = float(splitdata[5])
    f.close()
    cross_area=float(x_axis * y_axis)
    return cross_area

def check_input(param_type, param_value, filename, param_name):
    if not param_value:
    # If the parameter value was not defined	
         if param_type == 'string':
             return 'None'
         elif param_type == 'list':
             return []
         elif param_type == 'float':
             return 0.0
         elif param_type == 'int':
             return 0
    # If the value is provided in the input file, check for correct format
    if param_type == 'float':
        try:
            float(param_value)
        except ValueError:
            print(param_value)
            print('Use a floating point value for %s in the input file' % param_name)
            sys.exit()

        return float(param_value)

    elif param_type == 'int':
        try:
            int(param_value)
        except ValueError:
            print(param_value)
            print('Use an integer value for %s in the input file' % param_name)
            sys.exit()

        if int(param_value) < 0:
            print(param_value)
            print('Use a non-negative value for %s in the input file' % param_name)
            sys.exit()
        else:
            return int(param_value)

    elif (param_type == 'list') or (param_type == 'string'):
        return param_value

def num_to_mask(pdb_file):
    atm_num = []
    atm_num.append(0)
    # Get the relation between atom mask and atom index
    with open(pdb_file) as f_in:
      lines = (line.rstrip() for line in f_in)
      lines = list(line for line in lines if line)
      for i in range(0, len(lines)):
	if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
	  atomname = lines[i][12:16].strip() 
	  resid = lines[i][22:26].strip() 
	  number = lines[i][6:11].strip()  
	  atm_num.append(':'+resid+'@'+atomname)
    return atm_num
