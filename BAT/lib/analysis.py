#!/usr/bin/env python2
import glob as glob
import os as os
import re as re
import shutil as shutil
import signal as signal
import subprocess as sp
import sys as sys
import math
import numpy as np 
from lib.pymbar import MBAR # multistate Bennett acceptance ratio
from lib.pymbar import timeseries # timeseries analysis
from pathlib import Path

def fe_openmm(components, temperature, pose, dec_method, rest, attach_rest, lambdas, dic_itera1, dic_itera2, itera_steps, dt, dlambda, dec_int, weights, blocks):

    # Total simulation time
    total_time = 0
    for i in components:
      if i == 'a' or i == 'l' or i == 't' or i == 'c' or i == 'r' or i == 'm' or i == 'n':
        total_time = total_time + (dic_itera1[i]+dic_itera2[i])*itera_steps*len(attach_rest)*float(dt)/1000
      else:
        total_time = total_time + (dic_itera1[i]+dic_itera2[i])*itera_steps*len(lambdas)*float(dt)/1000
#    print(total_time)


    # Set initial values to zero
    fe_a = fe_bd = fe_t = fe_m = fe_n = fe_v = fe_e = fe_c = fe_r = fe_l = fe_f = fe_w = fe_vs = fe_es = 0
    fb_a = fb_bd = fb_t = fb_m = fb_n = fb_v = fb_e = fb_c = fb_r = fb_l = fb_f = fb_w = fb_es = fb_vs = 0
    sd_a = sd_bd = sd_t = sd_m = sd_n = sd_v = sd_e = sd_c = sd_r = sd_l = sd_f = sd_w = sd_vs = sd_es = 0

    # Get free energies for the whole run
    os.chdir('fe')
    os.chdir(pose)
    # Create Results folder
    if not os.path.exists('Results'):
      os.makedirs('Results')
    # Copy complex pdb structure
    shutil.copy('./build_files/complex.pdb','./Results/')
    for i in range(0, len(components)):
      comp = components[i]
      if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
        os.chdir('rest')
        os.chdir('%s-comp' %(comp))
        if (comp == 't' or comp == 'm'):
          # Calculate analytical release for dd and sdr
          with open('disang.rest', "r") as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if '#Lig_TR' in line)
            splitdata = lines[0].split()
            r0 = float(splitdata[6].strip(','))
            splitdata = lines[1].split()
            a1_0  = float(splitdata[6].strip(','))
            splitdata = lines[2].split()
            t1_0  = float(splitdata[6].strip(','))
            splitdata = lines[3].split()
            a2_0  = float(splitdata[6].strip(','))
            splitdata = lines[4].split()
            t2_0  = float(splitdata[6].strip(','))
            splitdata = lines[5].split()
            t3_0  = float(splitdata[6].strip(','))
            k_r = rest[2]
            k_a = rest[3]
            fe_bd = fe_int_op(r0, a1_0, t1_0, a2_0, t2_0, t3_0, k_r, k_a, temperature)
        out_file=Path('./output.dat')
        if out_file.exists(): 
          with open(out_file, "r") as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if line) # Non-blank lines in a list   
            for k in range(0, len(lines)):
              splitdata = lines[k].split()
              if (splitdata[0].strip() == 'Relative' and splitdata[6].strip() == 'whole'):
                if comp == 'c':
                  fe_c = -1.00*float(splitdata[9])
                elif comp == 'a':
                  fe_a = float(splitdata[9])
                elif comp == 't':
                  fe_t = float(splitdata[9])
                elif comp == 'n':
                  fe_n = -1.00*float(splitdata[9])
                elif comp == 'm':
                  fe_m = float(splitdata[9])
                elif comp == 'l':
                  fe_l = float(splitdata[9])
                elif comp == 'r':
                  fe_r = -1.00*float(splitdata[9])
      elif comp == 'e' or comp == 'v' or comp == 'f' or comp == 'w':
        os.chdir(dec_method)
        os.chdir('%s-comp' %(comp))
        out_file=Path('./output.dat')
        if dec_int == 'mbar':
          if out_file.exists(): 
            with open(out_file, "r") as f_in:
              lines = (line.rstrip() for line in f_in)
              lines = list(line for line in lines if line) # Non-blank lines in a list   
              for k in range(0, len(lines)):
                splitdata = lines[k].split()
                if (splitdata[0].strip() == 'Relative' and splitdata[6].strip() == 'whole'):
                  if comp == 'e' and dec_method == 'dd':
                    fe_e = -1.00*float(splitdata[9])
                  elif comp == 'v' and dec_method == 'dd':
                    fe_v = -1.00*float(splitdata[9])
                  elif comp == 'e' and dec_method == 'sdr':
                    fe_es = -1.00*float(splitdata[9])
                  elif comp == 'v' and dec_method == 'sdr':
                    fe_vs = -1.00*float(splitdata[9])
                  elif comp == 'f':
                    fe_f = float(splitdata[9])
                  elif comp == 'w':
                    fe_w = float(splitdata[9])
        elif dec_int == 'ti':
          ### Determine Number of windows
          K = 0
          filename = './'+comp+'%02.0f/output.dat' % K
          while os.path.isfile(filename):
            K = K+1
            filename = './'+comp+'%02.0f/output.dat' % K
          deltagop = 0
          deltagop_er = 0
          for k in range(K):
            filename = './'+comp+'%02.0f/output.dat' % k
            wind_d = 0
            wind_er = 0
            if  os.path.exists(filename):
              with open(filename, "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line) # Non-blank lines in a list   
                for j in range(0, len(lines)):
                  splitdata = lines[j].split()
                  if (splitdata[0].strip() == 'Relative' and splitdata[6].strip() == 'whole'):
                    wind_d = float(splitdata[9])
                    wind_er = float(splitdata[12])
            deltagop = deltagop + wind_d*weights[k]   
            deltagop_er = deltagop_er + wind_er*weights[k]   
          if comp == 'e' and dec_method == 'dd':
            fe_e = -1.00*float(deltagop/dlambda)
          elif comp == 'v' and dec_method == 'dd':
            fe_v = -1.00*float(deltagop/dlambda)
          elif comp == 'e' and dec_method == 'sdr':
            fe_es = -1.00*float(deltagop/dlambda)
          elif comp == 'v' and dec_method == 'sdr':
            fe_vs = -1.00*float(deltagop/dlambda)
          elif comp == 'f':
            fe_f = float(deltagop/dlambda)
          elif comp == 'w':
            fe_w = float(deltagop/dlambda)
      os.chdir('../../') 

    os.chdir('../../') 

    # Get MBAR free energy averages for the blocks
    os.chdir('fe')
    os.chdir(pose)
    blstd_a = []
    blstd_l = []
    blstd_t = []
    blstd_c = []
    blstd_r = []
    blstd_e = []
    blstd_v = []
    blstd_es = []
    blstd_vs = []
    blstd_f = []
    blstd_w = []
    blstd_m = []
    blstd_n = []
    for k in range(0, blocks):
      # Reset free energy values
      fb_a = fb_bd = fb_t = fb_m = fb_n = fb_v = fb_e = fb_c = fb_r = fb_l = fb_f = fb_w = fb_es = fb_vs = 0
      for i in range(0, len(components)):
        comp = components[i]
        if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
          os.chdir('rest')
          with open('./'+comp+'-comp/output.dat', "r") as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if 'Relative' in line and 'block' in line)
            splitdata = lines[k].split()
            if comp == 'c':
              fb_c = -1.00*float(splitdata[8])
              blstd_c.append(fb_c)
            elif comp == 'a':
              fb_a = float(splitdata[8])
              blstd_a.append(fb_a)
            elif comp == 't':
              fb_t = float(splitdata[8])
              blstd_t.append(fb_t)
            elif comp == 'n':
              fb_n = -1.00*float(splitdata[8])
              blstd_n.append(fb_n)
            elif comp == 'm':
              fb_m = float(splitdata[8])
              blstd_m.append(fb_m)
            elif comp == 'l':
              fb_l = float(splitdata[8])
              blstd_l.append(fb_l)
            elif comp == 'r':
              fb_r = -1.00*float(splitdata[8])
              blstd_r.append(fb_r)
          os.chdir('../')
        elif comp == 'v' or comp == 'e' or comp == 'f' or comp == 'w':
          os.chdir(dec_method)
          if dec_int == 'mbar': 
            with open('./'+comp+'-comp/output.dat', "r") as f_in:
              lines = (line.rstrip() for line in f_in)
              lines = list(line for line in lines if 'Relative' in line and 'block' in line)
              splitdata = lines[k].split()
              if comp == 'e' and dec_method == 'dd':
                fb_e = -1.00*float(splitdata[8])
                blstd_e.append(fb_e)
              if comp == 'e' and dec_method == 'sdr':
                fb_es = -1.00*float(splitdata[8])
                blstd_es.append(fb_es)
              elif comp == 'v' and dec_method == 'dd':
                fb_v = -1.00*float(splitdata[8])
                blstd_v.append(fb_v)
              elif comp == 'v' and dec_method == 'sdr':
                fb_vs = -1.00*float(splitdata[8])
                blstd_vs.append(fb_vs)
              if comp == 'f':
                fb_f = float(splitdata[8])
                blstd_f.append(fb_f)
              elif comp == 'w':
                fb_w = float(splitdata[8])
                blstd_w.append(fb_w)
          elif dec_int == 'ti': 
              ### Determine Number of windows
              K = 0
              filename = './'+comp+'-comp/'+comp+'%02.0f/output.dat' % K
              while os.path.isfile(filename):
                K = K+1
                filename = './'+comp+'-comp/'+comp+'%02.0f/output.dat' % K
              deltagop = 0
              for j in range(K):
                filename = './'+comp+'-comp/'+comp+'%02.0f/output.dat' % j
                wind_d = 0
                if  os.path.exists(filename):
                  with open(filename, "r") as f_in:
                    lines = (line.rstrip() for line in f_in)
                    lines = list(line for line in lines if 'Relative' in line and 'block' in line)
                    splitdata = lines[k].split()
                    wind_d = float(splitdata[8])
                    deltagop = deltagop + wind_d*weights[j]   
              if comp == 'e' and dec_method == 'dd':
                fb_e = -1.00*float(deltagop/dlambda)
                blstd_e.append(fb_e)
              elif comp == 'v' and dec_method == 'dd':
                fb_v = -1.00*float(deltagop/dlambda)
                blstd_v.append(fb_v)
              elif comp == 'e' and dec_method == 'sdr':
                fb_es = -1.00*float(deltagop/dlambda)
                blstd_es.append(fb_es)
              elif comp == 'v' and dec_method == 'sdr':
                fb_vs = -1.00*float(deltagop/dlambda)
                blstd_vs.append(fb_vs)
              elif comp == 'f':
                fb_f = float(deltagop/dlambda)
                blstd_f.append(fb_f)
              elif comp == 'w':
                fb_w = float(deltagop/dlambda)
                blstd_w.append(fb_w)
          os.chdir('../')

      # mevc modification
      if fb_m != 0 and fb_n == 0 and fb_c != 0:
        fb_n = fb_c 
        fb_c = 0 
        

      fb_bd = fe_bd
      blck_sdr = -1*(fb_a + fb_l + fb_t + fb_es + fb_vs + fb_bd + fb_c + fb_r)
      blck_dd = -1*(fb_a + fb_l + fb_t + fb_e + fb_v + fb_w + fb_f + fb_bd + fb_c + fb_r)
      blckm_dd = -1*(fb_m + fb_e + fb_v + fb_w + fb_f + fb_bd + fb_n)
      blckm_sdr = -1*(fb_m + fb_es + fb_vs + fb_bd + fb_n)

      # Write results for the blocks
      resfile = open('./Results/Res-b%02d.dat' %(k+1), 'w')
      if dec_method == 'dd':
        if fb_t != 0 or fb_c != 0 or fb_r != 0 or fb_a != 0 or fb_l != 0:
          resfile.write('\n----------------------------------------------\n')
          resfile.write('All components DD method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach protein CF;', fb_a))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand CF;', fb_l))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand TR;', fb_t))
          resfile.write('%-20s %8.2f\n' % ('Site Elect ('+dec_int.upper()+');', fb_e))
          resfile.write('%-20s %8.2f\n' % ('Site LJ ('+dec_int.upper()+');', fb_v))
          resfile.write('%-20s %8.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fb_w))
          resfile.write('%-20s %8.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fb_f))
          resfile.write('%-20s %8.2f\n' % ('Release ligand TR;',fb_bd))
          resfile.write('%-20s %8.2f\n' % ('Release ligand CF;', fb_c))
          resfile.write('%-20s %8.2f\n\n' % ('Release protein CF;', fb_r))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blck_dd))
        # Merged results
        if fb_m != 0 or fb_n != 0:
          fb_rel = fb_bd + fb_n
          resfile.write('\n----------------------------------------------\n')
          resfile.write('Merged components DD method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach all;', fb_m))
          resfile.write('%-20s %8.2f\n' % ('Site Elect ('+dec_int.upper()+');', fb_e))
          resfile.write('%-20s %8.2f\n' % ('Site LJ ('+dec_int.upper()+');', fb_v))
          resfile.write('%-20s %8.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fb_w))
          resfile.write('%-20s %8.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fb_f))
          resfile.write('%-20s %8.2f\n\n' % ('Release all;', fb_rel))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blckm_dd))
      if dec_method == 'sdr':
        if fb_t != 0 or fb_c != 0 or fb_r != 0 or fb_a != 0 or fb_l != 0:
          resfile.write('\n----------------------------------------------\n')
          resfile.write('All components SDR method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach protein CF;', fb_a))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand CF;', fb_l))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand TR;', fb_t))
          resfile.write('%-20s %8.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fb_es))
          resfile.write('%-20s %8.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fb_vs))
          resfile.write('%-20s %8.2f\n' % ('Release ligand TR;', fb_bd))
          resfile.write('%-20s %8.2f\n' % ('Release ligand CF;', fb_c))
          resfile.write('%-20s %8.2f\n\n' % ('Release protein CF;', fb_r))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blck_sdr))
        # Merged results
        if fb_m != 0 or fb_n != 0:
          fb_rel = fb_bd + fb_n
          resfile.write('\n----------------------------------------------\n')
          resfile.write('Merged components SDR method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach all;', fb_m))
          resfile.write('%-20s %8.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fb_es))
          resfile.write('%-20s %8.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fb_vs))
          resfile.write('%-20s %8.2f\n\n' % ('Release all;', fb_rel))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blckm_sdr))
      resfile.write('\n----------------------------------------------\n\n')
      resfile.write('Energies in kcal/mol\n')
      resfile.close()

    # Get sigmas of ther blocks

    for i in range(0, len(components)):
      comp = components[i]
      if comp == 'a':
        sd_a = np.std(blstd_a)
      if comp == 'l':
        sd_l = np.std(blstd_l)
      if comp == 't':
        sd_t = np.std(blstd_t)
      if comp == 'c':
        sd_c = np.std(blstd_c)
      if comp == 'r':
        sd_r = np.std(blstd_r)
      if comp == 'm':
        sd_m = np.std(blstd_m)
      if comp == 'n':
        sd_n = np.std(blstd_n)
      if comp == 'e' and dec_method == 'dd':
        sd_e = np.std(blstd_e)
      if comp == 'v' and dec_method == 'dd':
        sd_v = np.std(blstd_v)
      if comp == 'e' and dec_method == 'sdr':
        sd_es = np.std(blstd_es)
      if comp == 'v' and dec_method == 'sdr':
        sd_vs = np.std(blstd_vs)
      if comp == 'f':
        sd_f = np.std(blstd_f)
      if comp == 'w':
        sd_w = np.std(blstd_w)

    # mevc modification
    if fe_m != 0 and fe_n == 0 and fe_c != 0:
      fe_n = fe_c 
      sd_n = sd_c
      fe_c = 0 
      sd_c = 0 

    # Write final results
    total_dd = -1*(fe_a + fe_l + fe_t + fe_e + fe_v + fe_w + fe_f + fe_bd + fe_c + fe_r)
    merged_dd = -1*(fe_m + fe_e + fe_v + fe_w + fe_f + fe_bd + fe_n)
    total_sdr = -1*(fe_a + fe_l + fe_t + fe_es + fe_vs + fe_bd + fe_c + fe_r)
    merged_sdr = -1*(fe_m + fe_es + fe_vs + fe_bd + fe_n)
    sd_dd = math.sqrt(sd_a**2 + sd_l**2 + sd_t**2 + sd_e**2 + sd_v**2 + sd_w**2 + sd_f**2 + sd_bd**2 + sd_c**2 + sd_r**2)
    sd_merg_dd = math.sqrt(sd_m**2 + sd_e**2 + sd_v**2 + sd_w**2 + sd_f**2 + sd_bd**2 + sd_n**2)
    sd_sdr = math.sqrt(sd_a**2 + sd_l**2 + sd_t**2 + sd_es**2 + sd_vs**2 + sd_bd**2 + sd_c**2 + sd_r**2)
    sd_merg_sdr = math.sqrt(sd_m**2 + sd_es**2 + sd_vs**2 + sd_bd**2 + sd_n**2)

    resfile = open('./Results/Results.dat', 'w')
    if dec_method == 'dd':
      if fe_t != 0 or fe_c != 0 or fe_r != 0 or fe_a != 0 or fe_l != 0:
        resfile.write('\n----------------------------------------------\n')
        resfile.write('All components DD method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach protein CF;', fe_a, sd_a))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand CF;', fe_l, sd_l))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand TR;', fe_t, sd_t))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site Elect ('+dec_int.upper()+');', fe_e, sd_e))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site LJ ('+dec_int.upper()+');', fe_v, sd_v))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fe_w, sd_w))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fe_f, sd_f))
        resfile.write('%-20s %8.2f;    \n' % ('Release ligand TR;',fe_bd))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Release ligand CF;', fe_c, sd_c))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release protein CF;', fe_r, sd_r))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', total_dd, sd_dd))
      # Merged results
      if fe_m != 0 or fe_n != 0:
        fe_rel = fe_bd + fe_n
        resfile.write('\n----------------------------------------------\n')
        resfile.write('Merged components DD method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach all;', fe_m, sd_m))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site Elect ('+dec_int.upper()+');', fe_e, sd_e))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site LJ ('+dec_int.upper()+');', fe_v, sd_v))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fe_w, sd_w))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fe_f, sd_f))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release all;', fe_rel, sd_n))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', merged_dd, sd_merg_dd))
    if dec_method == 'sdr':
      if fe_t != 0 or fe_c != 0 or fe_r != 0 or fe_a != 0 or fe_l != 0:
        resfile.write('\n----------------------------------------------\n')
        resfile.write('All components SDR method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach protein CF;', fe_a, sd_a))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand CF;', fe_l, sd_l))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand TR;', fe_t, sd_t))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fe_es, sd_es))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fe_vs, sd_vs))
        resfile.write('%-20s %8.2f;    \n' % ('Release ligand TR;',fe_bd))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Release ligand CF;', fe_c, sd_c))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release protein CF;', fe_r, sd_r))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', total_sdr, sd_sdr))
      # Merged results
      if fe_m != 0 or fe_n != 0:
        fe_rel = fe_bd + fe_n
        resfile.write('\n----------------------------------------------\n')
        resfile.write('Merged components SDR method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach all;', fe_m, sd_m))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fe_es, sd_es))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fe_vs, sd_vs))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release all;', fe_rel, sd_n))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', merged_sdr, sd_merg_sdr))
    resfile.write('\n----------------------------------------------\n\n')
    resfile.write('Energies in kcal/mol\n\n')
    resfile.write('Total simulation time (based on input file): %6.1f nanoseconds\n\n' % total_time)
    resfile.write('Please cite:\n\n')
    resfile.write('G. Heinzelmann and M. K. Gilson (2021). “Automation of absolute protein-ligand binding free energy calculations for docking refinement and compound evaluation”. Scientific Reports, 11, 1116.\n\n')
    resfile.write('D. J. Huggins (2022) "Comparing the Performance of Different AMBER Protein Forcefields, Partial Charge Assignments, and Water Models for Absolute Binding Free Energy Calculations." Journal of Chemical Theory and Computation, 18, 2616.\n\n')
    resfile.close()


def fe_values(blocks, components, temperature, pose, attach_rest, lambdas, weights, dec_int, dec_method, rest, dic_steps1, dic_steps2, dt):

    # Total simulation time
    total_time = 0
    for i in components:
      if i == 'a' or i == 'l' or i == 't' or i == 'c' or i == 'r' or i == 'm' or i == 'n':
        total_time = total_time + (dic_steps1[i]+dic_steps2[i])*len(attach_rest)*float(dt)/1000
      else:
        total_time = total_time + (dic_steps1[i]+dic_steps2[i])*len(lambdas)*float(dt)/1000
#    print(total_time)
         

    # Set initial values to zero
    fe_a = fe_bd = fe_t = fe_m = fe_n = fe_v = fe_e = fe_c = fe_r = fe_l = fe_f = fe_w = fe_vs = fe_es = 0
    fb_a = fb_bd = fb_t = fb_m = fb_n = fb_v = fb_e = fb_c = fb_r = fb_l = fb_f = fb_w = fb_es = fb_vs = 0
    sd_a = sd_bd = sd_t = sd_m = sd_n = sd_v = sd_e = sd_c = sd_r = sd_l = sd_f = sd_w = sd_vs = sd_es = 0

    # Acquire simulation data
    os.chdir('fe')
    os.chdir(pose)
    for i in range(0, len(components)):
      comp = components[i]
      if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
        os.chdir('rest')
        for j in range(0, len(attach_rest)):
          data = []
          win = j
          os.chdir('%s%02d' %(comp, int(win)))
          if (comp == 't' or comp == 'm') and win == 0:
            # Calculate analytical release for dd and sdr
            with open('disang.rest', "r") as f_in:
              lines = (line.rstrip() for line in f_in)
              lines = list(line for line in lines if '#Lig_TR' in line)
              splitdata = lines[0].split()
              r0 = float(splitdata[6].strip(','))
              splitdata = lines[1].split()
              a1_0  = float(splitdata[6].strip(','))
              splitdata = lines[2].split()
              t1_0  = float(splitdata[6].strip(','))
              splitdata = lines[3].split()
              a2_0  = float(splitdata[6].strip(','))
              splitdata = lines[4].split()
              t2_0  = float(splitdata[6].strip(','))
              splitdata = lines[5].split()
              t3_0  = float(splitdata[6].strip(','))
              k_r = rest[2]
              k_a = rest[3]
              fe_bd = fe_int(r0, a1_0, t1_0, a2_0, t2_0, t3_0, k_r, k_a, temperature)
          # Get restraint trajectory file
          sp.call('cpptraj -i restraints.in > restraints.log 2>&1', shell=True)
          # Separate in blocks
          with open("restraints.dat", "r") as fin:
            for line in fin:
              if not '#' in line:
                data.append(line)
          for k in range(0, blocks):
            fout = open('rest%02d.dat' % (k+1), "w")
            for t in range(k*int(round(len(data)//blocks)), (k+1)*int(round(len(data)//blocks))):
              fout.write(data[t])
            fout.close()
          os.chdir('../')
      elif comp == 'e' or comp == 'v' or comp == 'f' or comp == 'w':
        os.chdir(dec_method)
        if dec_int == 'ti':
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
        elif dec_int == 'mbar':
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
      os.chdir('../')

    os.chdir('../../')

    # Get free energies for the whole run    
    for i in range(0, len(components)):
      comp = components[i]
      if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
        rest_file = 'restraints.dat'
        mode = 'all'
        fe_mbar(comp, pose, mode, rest_file, temperature)
        mode = 'sub'
        fe_mbar(comp, pose, mode, rest_file, temperature)
      else:
        if dec_int == 'ti':
          rest_file = 'dvdl.dat'
          mode = 'all'
          fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)
        elif dec_int == 'mbar':
          rest_file = 'energies.dat'
          mode = 'all'
          fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)
          mode = 'sub'
          fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

    # Get free energies for the blocks
    for i in range(0, len(components)):
      for k in range(0, blocks):
        comp = components[i]
        if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
          rest_file = 'rest%02d.dat' % (k+1)
          mode = 'b%02d' % (k+1)
          fe_mbar(comp, pose, mode, rest_file, temperature)
        else:
          if dec_int == 'ti':
            rest_file = 'dvdl%02d.dat' % (k+1) 
            mode = 'b%02d' % (k+1) 
            fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)
          elif dec_int == 'mbar':
            rest_file = 'ener%02d.dat' % (k+1) 
            mode = 'b%02d' % (k+1) 
            fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature)

    sys.stdout = sys.__stdout__

    # Calculate final results
    os.chdir('fe')
    os.chdir(pose)
    # Get MBAR free energy averages
    for i in range(0, len(components)):
      comp = components[i]
      if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
        os.chdir('rest')
        with open('./data/mbar-'+comp+'-all.dat', "r") as f_in:
          lines = (line.rstrip() for line in f_in)
          lines = list(line for line in lines if line)
          data = lines[-1]
          splitdata = data.split()
          if comp == 'c':
            fe_c = -1.00*float(splitdata[1])
          elif comp == 'a':
            fe_a = float(splitdata[1])
          elif comp == 't':
            fe_t = float(splitdata[1])
          elif comp == 'n':
            fe_n = -1.00*float(splitdata[1])
          elif comp == 'm':
            fe_m = float(splitdata[1])
          elif comp == 'l':
            fe_l = float(splitdata[1])
          elif comp == 'r':
            fe_r = -1.00*float(splitdata[1])
        os.chdir('../')
      elif comp == 'v' or comp == 'e' or comp == 'f' or comp == 'w':
        os.chdir(dec_method)
        with open('./data/'+dec_int+'-'+comp+'-all.dat', "r") as f_in:
          lines = (line.rstrip() for line in f_in)
          lines = list(line for line in lines if line)
          data = lines[-1]
          splitdata = data.split()
          if comp == 'e' and dec_method == 'dd':
            fe_e = float(splitdata[1])
          if comp == 'e' and dec_method == 'sdr':
            fe_es = float(splitdata[1])
          elif comp == 'v' and dec_method == 'dd':
            fe_v = float(splitdata[1])
          elif comp == 'v' and dec_method == 'sdr':
            fe_vs = float(splitdata[1])
          if comp == 'w':
            fe_w = -1.00*float(splitdata[1])
          elif comp == 'f':
            fe_f = -1.00*float(splitdata[1])
        os.chdir('../')


    # Get errors
    for i in range(0, len(components)):
      comp = components[i]
      if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
        os.chdir('rest')
        b_data = [] 
        for k in range(0, blocks):
          with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if line)
            data = lines[-1]
            splitdata = data.split()
            b_data.append(float(splitdata[1]))
          if comp == 'c':
            sd_c = np.std(b_data)
          elif comp == 'a':
            sd_a = np.std(b_data)
          elif comp == 't':
            sd_t = np.std(b_data)
          elif comp == 'm':
            sd_m = np.std(b_data)
          elif comp == 'n':
            sd_n = np.std(b_data)
          elif comp == 'l':
            sd_l = np.std(b_data)
          elif comp == 'r':
            sd_r = np.std(b_data)
        os.chdir('../')
      elif comp == 'e' or comp == 'v' or comp == 'f' or comp == 'w': 
        os.chdir(dec_method)
        if dec_int == 'mbar':
          if comp == 'e' and dec_method == 'dd':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_e = np.std(b_data)
          if comp == 'e' and dec_method == 'sdr':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_es = np.std(b_data)
          if comp == 'v' and dec_method == 'sdr':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_vs = np.std(b_data)
          if comp == 'v' and dec_method == 'dd':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_v = np.std(b_data)
          if comp == 'f':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_f = np.std(b_data)
          if comp == 'w':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_w = np.std(b_data)
        elif dec_int == 'ti':
          if comp == 'e' and dec_method == 'dd':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_e = np.std(b_data)
          if comp == 'e' and dec_method == 'sdr':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_es = np.std(b_data)
          if comp == 'v' and dec_method == 'dd':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_v = np.std(b_data)
          if comp == 'v' and dec_method == 'sdr':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_vs = np.std(b_data)
          if comp == 'f':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_f = np.std(b_data)
          if comp == 'w':
            b_data = [] 
            for k in range(0, blocks):
              with open('./data/ti-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
                lines = (line.rstrip() for line in f_in)
                lines = list(line for line in lines if line)
                data = lines[-1]
                splitdata = data.split()
                b_data.append(float(splitdata[1]))
            sd_w = np.std(b_data)
        os.chdir('../')

    # Create Results folder
    if not os.path.exists('Results'):
      os.makedirs('Results')

    # Copy complex pdb structure
      shutil.copy('./build_files/complex.pdb','./Results/')

    # Get MBAR free energy averages for the blocks
    for k in range(0, blocks):
      # Reset free energy values
      fb_a = fb_bd = fb_t = fb_m = fb_n = fb_v = fb_e = fb_c = fb_r = fb_l = fb_f = fb_w = fb_es = fb_vs = 0
      for i in range(0, len(components)):
        comp = components[i]
        if comp == 'a' or comp == 'l' or comp == 't' or comp == 'c' or comp == 'r' or comp == 'm' or comp == 'n':
          os.chdir('rest')
          with open('./data/mbar-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if line)
            data = lines[-1]
            splitdata = data.split()
            if comp == 'c':
              fb_c = -1.00*float(splitdata[1])
            elif comp == 'a':
              fb_a = float(splitdata[1])
            elif comp == 't':
              fb_t = float(splitdata[1])
            elif comp == 'n':
              fb_n = -1.00*float(splitdata[1])
            elif comp == 'm':
              fb_m = float(splitdata[1])
            elif comp == 'l':
              fb_l = float(splitdata[1])
            elif comp == 'r':
              fb_r = -1.00*float(splitdata[1])
          os.chdir('../')
        elif comp == 'v' or comp == 'e' or comp == 'f' or comp == 'w':
          os.chdir(dec_method)
          with open('./data/'+dec_int+'-'+comp+'-b%02d.dat' %(k+1), "r") as f_in:
            lines = (line.rstrip() for line in f_in)
            lines = list(line for line in lines if line)
            data = lines[-1]
            splitdata = data.split()
            if comp == 'e' and dec_method == 'dd':
              fb_e = float(splitdata[1])
            if comp == 'e' and dec_method == 'sdr':
              fb_es = float(splitdata[1])
            elif comp == 'v' and dec_method == 'dd':
              fb_v = float(splitdata[1])
            elif comp == 'v' and dec_method == 'sdr':
              fb_vs = float(splitdata[1])
            if comp == 'f':
              fb_f = -1.00*float(splitdata[1])
            elif comp == 'w':
              fb_w = -1.00*float(splitdata[1])
          os.chdir('../')

      # mevc modification
      if fb_m != 0 and fb_n == 0 and fb_c != 0:
        fb_n = fb_c
        fb_c = 0

      fb_bd = fe_bd
      blck_sdr = -1*(fb_a + fb_l + fb_t + fb_es + fb_vs + fb_bd + fb_c + fb_r)
      blck_dd = -1*(fb_a + fb_l + fb_t + fb_e + fb_v + fb_w + fb_f + fb_bd + fb_c + fb_r)
      blckm_dd = -1*(fb_m + fb_e + fb_v + fb_w + fb_f + fb_bd + fb_n)
      blckm_sdr = -1*(fb_m + fb_es + fb_vs + fb_bd + fb_n)

      # Write results for the blocks
      resfile = open('./Results/Res-b%02d.dat' %(k+1), 'w')
      if dec_method == 'dd' and os.path.exists('./dd/data/'):
        if fb_t != 0 or fb_c != 0 or fb_r != 0 or fb_a != 0 or fb_l != 0:
          resfile.write('\n----------------------------------------------\n')
          resfile.write('All components DD method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach protein CF;', fb_a))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand CF;', fb_l))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand TR;', fb_t))
          resfile.write('%-20s %8.2f\n' % ('Site Elect ('+dec_int.upper()+');', fb_e))
          resfile.write('%-20s %8.2f\n' % ('Site LJ ('+dec_int.upper()+');', fb_v))
          resfile.write('%-20s %8.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fb_w))
          resfile.write('%-20s %8.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fb_f))
          resfile.write('%-20s %8.2f\n' % ('Release ligand TR;',fb_bd))
          resfile.write('%-20s %8.2f\n' % ('Release ligand CF;', fb_c))
          resfile.write('%-20s %8.2f\n\n' % ('Release protein CF;', fb_r))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blck_dd))
        # Merged results
        if fb_m != 0 or fb_n != 0:
          fb_rel = fb_bd + fb_n
          resfile.write('\n----------------------------------------------\n')
          resfile.write('Merged components DD method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach all;', fb_m))
          resfile.write('%-20s %8.2f\n' % ('Site Elect ('+dec_int.upper()+');', fb_e))
          resfile.write('%-20s %8.2f\n' % ('Site LJ ('+dec_int.upper()+');', fb_v))
          resfile.write('%-20s %8.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fb_w))
          resfile.write('%-20s %8.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fb_f))
          resfile.write('%-20s %8.2f\n\n' % ('Release all;', fb_rel))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blckm_dd))
      if dec_method == 'sdr' and os.path.exists('./sdr/data/'):
        if fb_t != 0 or fb_c != 0 or fb_r != 0 or fb_a != 0 or fb_l != 0:
          resfile.write('\n----------------------------------------------\n')
          resfile.write('All components SDR method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach protein CF;', fb_a))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand CF;', fb_l))
          resfile.write('%-20s %8.2f\n' % ('Attach ligand TR;', fb_t))
          resfile.write('%-20s %8.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fb_es))
          resfile.write('%-20s %8.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fb_vs))
          resfile.write('%-20s %8.2f\n' % ('Release ligand TR;', fb_bd))
          resfile.write('%-20s %8.2f\n' % ('Release ligand CF;', fb_c))
          resfile.write('%-20s %8.2f\n\n' % ('Release protein CF;', fb_r))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blck_sdr))
        # Merged results
        if fb_m != 0 or fb_n != 0:
          fb_rel = fb_bd + fb_n
          resfile.write('\n----------------------------------------------\n')
          resfile.write('Merged components SDR method')
          resfile.write('\n----------------------------------------------\n\n')
          resfile.write('%-20s %-10s\n\n' % ('Component', 'Free Energy'))
          resfile.write('%-20s %8.2f\n' % ('Attach all;', fb_m))
          resfile.write('%-20s %8.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fb_es))
          resfile.write('%-20s %8.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fb_vs))
          resfile.write('%-20s %8.2f\n\n' % ('Release all;', fb_rel))
          resfile.write('%-20s %8.2f\n' % ('Binding free energy;', blckm_sdr))
      resfile.write('\n----------------------------------------------\n\n')
      resfile.write('Energies in kcal/mol\n')
      resfile.close()

    # mevc modification
    if fe_m != 0 and fe_n == 0 and fe_c != 0:
      fe_n = fe_c
      sd_n = sd_c
      fe_c = 0
      sd_c = 0

    # Write final results
    total_dd = -1*(fe_a + fe_l + fe_t + fe_e + fe_v + fe_w + fe_f + fe_bd + fe_c + fe_r)
    merged_dd = -1*(fe_m + fe_e + fe_v + fe_w + fe_f + fe_bd + fe_n)
    total_sdr = -1*(fe_a + fe_l + fe_t + fe_es + fe_vs + fe_bd + fe_c + fe_r)
    merged_sdr = -1*(fe_m + fe_es + fe_vs + fe_bd + fe_n)
    sd_dd = math.sqrt(sd_a**2 + sd_l**2 + sd_t**2 + sd_e**2 + sd_v**2 + sd_w**2 + sd_f**2 + sd_bd**2 + sd_c**2 + sd_r**2)
    sd_merg_dd = math.sqrt(sd_m**2 + sd_e**2 + sd_v**2 + sd_w**2 + sd_f**2 + sd_bd**2 + sd_n**2)
    sd_sdr = math.sqrt(sd_a**2 + sd_l**2 + sd_t**2 + sd_es**2 + sd_vs**2 + sd_bd**2 + sd_c**2 + sd_r**2)
    sd_merg_sdr = math.sqrt(sd_m**2 + sd_es**2 + sd_vs**2 + sd_bd**2 + sd_n**2)

    resfile = open('./Results/Results.dat', 'w')
    if dec_method == 'dd' and os.path.exists('./dd/data/'):
      if fe_t != 0 or fe_c != 0 or fe_r != 0 or fe_a != 0 or fe_l != 0:
        resfile.write('\n----------------------------------------------\n')
        resfile.write('All components DD method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach protein CF;', fe_a, sd_a))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand CF;', fe_l, sd_l))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand TR;', fe_t, sd_t))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site Elect ('+dec_int.upper()+');', fe_e, sd_e))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site LJ ('+dec_int.upper()+');', fe_v, sd_v))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fe_w, sd_w))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fe_f, sd_f))
        resfile.write('%-20s %8.2f;    \n' % ('Release ligand TR;',fe_bd))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Release ligand CF;', fe_c, sd_c))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release protein CF;', fe_r, sd_r))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', total_dd, sd_dd))
      # Merged results
      if fe_m != 0 or fe_n != 0:
        fe_rel = fe_bd + fe_n
        resfile.write('\n----------------------------------------------\n')
        resfile.write('Merged components DD method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach all;', fe_m, sd_m))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site Elect ('+dec_int.upper()+');', fe_e, sd_e))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Site LJ ('+dec_int.upper()+');', fe_v, sd_v))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk LJ ('+dec_int.upper()+');', fe_w, sd_w))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Bulk Elect ('+dec_int.upper()+');', fe_f, sd_f))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release all;', fe_rel, sd_n))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', merged_dd, sd_merg_dd))
    if dec_method == 'sdr' and os.path.exists('./sdr/data/'):
      if fe_t != 0 or fe_c != 0 or fe_r != 0 or fe_a != 0 or fe_l != 0:
        resfile.write('\n----------------------------------------------\n')
        resfile.write('All components SDR method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach protein CF;', fe_a, sd_a))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand CF;', fe_l, sd_l))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach ligand TR;', fe_t, sd_t))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fe_es, sd_es))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fe_vs, sd_vs))
        resfile.write('%-20s %8.2f;    \n' % ('Release ligand TR;',fe_bd))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Release ligand CF;', fe_c, sd_c))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release protein CF;', fe_r, sd_r))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', total_sdr, sd_sdr))
      # Merged results
      if fe_m != 0 or fe_n != 0:
        fe_rel = fe_bd + fe_n
        resfile.write('\n----------------------------------------------\n')
        resfile.write('Merged components SDR method')
        resfile.write('\n----------------------------------------------\n\n')
        resfile.write('%-20s %-10s %-4s\n\n' % ('Component', 'Free Energy;', 'Sigma'))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Attach all;', fe_m, sd_m))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Electrostatic ('+dec_int.upper()+');', fe_es, sd_es))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Lennard-Jones ('+dec_int.upper()+');', fe_vs, sd_vs))
        resfile.write('%-20s %8.2f;    %3.2f\n\n' % ('Release all;', fe_rel, sd_n))
        resfile.write('%-20s %8.2f;    %3.2f\n' % ('Binding free energy;', merged_sdr, sd_merg_sdr))
    resfile.write('\n----------------------------------------------\n\n')
    resfile.write('Energies in kcal/mol\n\n')
    resfile.write('Total simulation time (based on input file): %6.1f nanoseconds\n\n' % total_time)
    resfile.write('Please cite:\n\n')
    resfile.write('G. Heinzelmann and M. K. Gilson (2021). “Automation of absolute protein-ligand binding free energy calculations for docking refinement and compound evaluation”. Scientific Reports, 11, 1116.\n\n')
    resfile.write('D. J. Huggins (2022) "Comparing the Performance of Different AMBER Protein Forcefields, Partial Charge Assignments, and Water Models for Absolute Binding Free Energy Calculations." Journal of Chemical Theory and Computation, 18, 2616.\n\n')
    resfile.close()


def fe_mbar(comp, pose, mode, rest_file, temperature):

    kB = 1.381e-23 * 6.022e23 / (4.184 * 1000.0) # Boltzmann constant in kJ/mol/K
    beta = 1/(kB * temperature) # beta
    N_max = 20000 # Max frames for any simulation window, you should check this if you did some long runs


    ### Change to pose directory
    os.chdir('fe')
    os.chdir(pose)
    os.chdir('rest')
    if not os.path.exists('data'):
      os.makedirs('data')

    # Define log file
    sys.stdout = open('./data/mbar-'+comp+'-'+mode+'.log', 'w')

    ### Determine Number of windows
    K = 0
    filename = './'+comp+'%02.0f/%s' % (K, rest_file)
    while os.path.isfile(filename):
      K = K+1
      filename = './'+comp+'%02.0f/%s' % (K, rest_file)

    ## Determine Number of restraints
    infile = open('./'+comp+'00/disang.rest', 'r')
    disang = infile.readlines()
    infile.close()
    R = 0
    if (comp == 't'):
      for line in disang:
        cols = line.split()
        if len(cols) != 0 and (cols[-1] == "#Lig_TR"):
          R += 1
    elif (comp == 'l' or comp == 'c'):
      for line in disang:
        cols = line.split()
        if len(cols) != 0 and (cols[-1] == "#Lig_C" or cols[-1] == "#Lig_D"):
          R += 1
    elif (comp == 'a' or comp == 'r'):
      for line in disang:
        cols = line.split()
        if len(cols) != 0 and (cols[-1] == "#Rec_C" or cols[-1] == "#Rec_D"):
          R += 1
    elif (comp == 'm' or comp == 'n'):
      for line in disang:
        cols = line.split()
        if len(cols) != 0 and (cols[-1] == "#Rec_C" or cols[-1] == "#Rec_D" or cols[-1] == "#Lig_TR" or cols[-1] == "#Lig_C" or cols[-1] == "#Lig_D"):
          R += 1

    print  ("K= %5.0f  R= %5.0f" % ( K, R ))

    ### Calculate Statistical Inefficiency (g)
    def calcg(data):
      sum = 0
      randnum = ("%05.0f" % (int(100000*np.random.random())))
      datafn = '/dev/shm/series.'+randnum+'.dat'
      acffn = '/dev/shm/acf.'+randnum+'.dat'
      cppfn = '/dev/shm/pt-acf.'+randnum+'.in'
      np.savetxt(datafn,data)
      cpptin = open(cppfn, 'w')
      cpptin.write("readdata "+datafn+" name "+randnum+"\nautocorr "+randnum+" out "+acffn+" noheader\n")
      cpptin.close()

      FNULL = open(os.devnull, 'w')
      sp.call(['cpptraj','-i',cppfn], stdout=FNULL, stderr=sp.STDOUT)

      with open(acffn, 'r') as acf:
        for line in acf:
          col = line.split()
          t = float(col[0]) - 1.0
      T = t

      with open(acffn, 'r') as acf:
        for line in acf:
          col = line.split()
          t = float(col[0]) - 1.0
          v = float(col[1])
          if t == 0:
            continue
          if v < 0.0:
            break
          sum += ( 1 - (t/T) )*(v)

      sp.call(['rm',datafn,acffn,cppfn])

      return 1+(2*sum)

    ### Allocate storage for simulation data
    N = np.zeros([K], np.int32)                       # N_k[k] is the number of snapshots to be used from umbrella simulation k
    Neff = np.zeros([K], np.int32)
    Nind = np.zeros([K], np.int32)
    rty = ['d']*R                                     # restraint type (distance or angle)
    rfc = np.zeros([K,R], np.float64)                 # restraint force constant
    req = np.zeros([K,R], np.float64)                 # restraint target value
    val = np.zeros([N_max,K,R], np.float64)           # value of the restrained variable at each frame n
    g = np.zeros([K], np.float64)
    u=np.zeros([N_max], np.float64)

    ### Read the simulation data
    for k in range(K):
      # Read Equilibrium Value and Force Constant
      filename = './'+comp+'%02.0f/disang.rest' % k
      infile = open(filename, 'r')
      disang = infile.readlines()
      infile.close()
      r = 0
      for line in disang:
        cols = line.split()
        if (comp == 't'):
          if len(cols) != 0 and (cols[-1] == "#Lig_TR"):
            natms = len(cols[2].split(','))-1
            req[k,r] = float(cols[6].replace(",", ""))
            if natms == 2:
              rty[r] = 'd'
              rfc[k,r] = float(cols[12].replace(",", ""))
            elif natms == 3:
              rty[r] = 'a'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            elif natms == 4:
              rty[r] = 't'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            else:
              sys.exit("not sure about restraint type!")
            r += 1
        elif (comp == 'l' or comp == 'c'):
          if len(cols) != 0 and (cols[-1] == "#Lig_C" or cols[-1] == "#Lig_D"):
            natms = len(cols[2].split(','))-1
            req[k,r] = float(cols[6].replace(",", ""))
            if natms == 2:
              rty[r] = 'd'
              rfc[k,r] = float(cols[12].replace(",", ""))
            elif natms == 3:
              rty[r] = 'a'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            elif natms == 4:
              rty[r] = 't'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            else:
              sys.exit("not sure about restraint type!")
            r += 1
        elif (comp == 'a' or comp == 'r'):
          if len(cols) != 0 and (cols[-1] == "#Rec_C" or cols[-1] == "#Rec_D"):
            natms = len(cols[2].split(','))-1
            req[k,r] = float(cols[6].replace(",", ""))
            if natms == 2:
              rty[r] = 'd'
              rfc[k,r] = float(cols[12].replace(",", ""))
            elif natms == 3:
              rty[r] = 'a'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            elif natms == 4:
              rty[r] = 't'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            else:
              sys.exit("not sure about restraint type!")
            r += 1
        elif (comp == 'm' or comp == 'n'):
          if len(cols) != 0 and (cols[-1] == "#Rec_C" or cols[-1] == "#Rec_D" or cols[-1] == "#Lig_TR" or cols[-1] == "#Lig_C" or cols[-1] == "#Lig_D"):
            natms = len(cols[2].split(','))-1
            req[k,r] = float(cols[6].replace(",", ""))
            if natms == 2:
              rty[r] = 'd'
              rfc[k,r] = float(cols[12].replace(",", ""))
            elif natms == 3:
              rty[r] = 'a'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            elif natms == 4:
              rty[r] = 't'
              rfc[k,r] = float(cols[12].replace(",", ""))*(np.pi/180.0)*(np.pi/180.0)  ### Convert to degrees
            else:
              sys.exit("not sure about restraint type!")
            r += 1

      # Read in Values for restrained variables for each simulation
      filename = './'+comp+'%02.0f/%s' % (k, rest_file)
      infile = open(filename, 'r')
      restdat = infile.readlines()     # slice off first 20 lines  readlines()[20:]
      infile.close()
      # Parse Data
      n = 0
      for line in restdat:
        if line[0] != '#' and line[0] != '@' and n < N_max:
          cols = line.split()
          for r in range(R):
            if rty[r] == 't': # Do phase corrections
              tmp = float(cols[r+1])
              if tmp < req[k,r]-180.0:
                val[n,k,r] = tmp + 360
              elif tmp > req[k,r]+180.0:
                val[n,k,r] = tmp - 360
              else:
                val[n,k,r] = tmp
            else:
              val[n,k,r] = float(cols[r+1])
          n += 1
       
      N[k] = n

      # Calculate Reduced Potential 
      if comp != 'u': ### Attach/Release Restraints
        if rfc[k,0] == 0:
          tmp=np.ones([R],np.float64)*0.001  ########## CHECK THIS!! might interfere on protein attach
          u[0:N[k]] = np.sum(beta*tmp[0:R]*((val[0:N[k],k,0:R]-req[k,0:R])**2), axis=1)
        else:
          u[0:N[k]] = np.sum(beta*rfc[k,0:R]*((val[0:N[k],k,0:R]-req[k,0:R])**2), axis=1)
      else: ### Umbrella/Translation
        u[0:N[k]] = (beta*rfc[k,0]*((val[0:N[k],k,0]-req[k,0])**2))


      if mode == 'sub':
        g[k] = calcg(u[0:N[k]])
        subs = timeseries.subsampleCorrelatedData(np.zeros([N[k]]),g=g[k])
        Nind[k] = len(subs)
        Neff[k] = Nind[k]
      else:
        g[k] = 1.00
        Neff[k] = N[k]

      print  ("Processed Window %5.0f.  N= %12.0f.  g= %10.3f   Neff= %12.0f" % ( k, N[k], g[k], Neff[k] ))

    Upot = np.zeros([K,K,np.max(Neff)], np.float64)

    # Calculate Restraint Energy
    for k in range(K):
      if mode == 'sub': #subsampling
        subs = timeseries.subsampleCorrelatedData(np.zeros([N[k]]),g=g[k])
        for l in range(K):
          if comp != 'u': # Attach Restraints
            Upot[k,l,0:Neff[k]] = np.sum(beta*rfc[l,0:R]*((val[subs[0:Neff[k]],k,0:R]-req[l,0:R])**2), axis=1)
          else: # Umbrella/Translation
            Upot[k,l,0:Neff[k]] = (beta*rfc[l,0]*((val[subs[0:Neff[k]],k,0]-req[l,0])**2))
      else:
        Neff[k] = N[k]
        for l in range(K): # all samples
          if comp != 'u': # Attach Restraints
            Upot[k,l,0:Neff[k]] = np.sum(beta*rfc[l,0:R]*((val[0:Neff[k],k,0:R]-req[l,0:R])**2), axis=1)
          else: # Umbrella/Translation
            Upot[k,l,0:Neff[k]] = (beta*rfc[l,0]*((val[0:Neff[k],k,0]-req[l,0])**2))

    val=[]

    print  ("Running MBAR... ") 
    mbar = MBAR(Upot, Neff)

    print  ("Calculate Free Energy Differences Between States")
    [Deltaf, dDeltaf] = mbar.getFreeEnergyDifferences()

    min = np.argmin(Deltaf[0])

    # Write to file
    print  ("Free Energy Differences (in units of kcal/mol)")
    print  ("%9s %8s %8s %12s %12s" % ('bin', 'f', 'df', 'deq', 'dfc'))
    datfile = open('./data/mbar-'+comp+'-'+mode+'.dat', 'w')
    for k in range(K):
      if comp != 'u': # Attach/release
        print ("%10.5f %10.5f %10.5f %12.7f %12.7f" % ( rfc[k,0]/rfc[-1,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] ))
        datfile.write ( "%10.5f %10.5f %10.5f %12.7f %12.7f\n" % ( rfc[k,0]/rfc[-1,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] ) )
      else: # Umbrella/Translation
        print ("%10.5f %10.5f %10.5f %12.7f %12.7f" % ( req[k,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] ))
        datfile.write ( "%10.5f %10.5f %10.5f %12.7f %12.7f\n" % ( req[k,0], Deltaf[0,k]/beta, dDeltaf[0,k]/beta, req[k,0], rfc[k,0] ) )
    datfile.close()
    print ("\n\n")
    
    os.chdir('../../../')


def fe_int(r1_0, a1_0, t1_0, a2_0, t2_0, t3_0, k_r, k_a, temperature):

    R = 1.987204118e-3 # kcal/mol-K, a.k.a. boltzman constant
    beta = 1/(temperature*R)
    r1lb,r1ub,r1st = [0.0,100.0,0.0001]
    a1lb,a1ub,a1st = [0.0,np.pi,0.00005] 
    t1lb,t1ub,t1st = [-np.pi,np.pi,0.00005] 
    a2lb,a2ub,a2st = [0.0,np.pi,0.00005] 
    t2lb,t2ub,t2st = [-np.pi,np.pi,0.00005] 
    t3lb,t3ub,t3st = [-np.pi,np.pi,0.00005] 

    def dih_per(lb,ub,st,t_0):
      drange = np.arange(lb,ub,st) 
      delta = (drange-np.radians(t_0))   
      for i in range(0, len(delta)):
        if delta[i] >= np.pi:
          delta[i] = delta[i]-(2*np.pi)
        if delta[i] <= -np.pi:
          delta[i] = delta[i]+(2*np.pi)
      return delta
     
    def f_r1(val):
      return (val**2)*np.exp(-beta*k_r*(val-r1_0)**2)
    def f_a1(val):
      return np.sin(val)*np.exp(-beta*k_a*(val-np.radians(a1_0))**2)
    def f_a2(val):
      return np.sin(val)*np.exp(-beta*k_a*(val-np.radians(a2_0))**2)
    def f_t1(delta):
      return np.exp(-beta*k_a*(delta)**2)
    def f_t2(delta):
      return np.exp(-beta*k_a*(delta)**2)
    def f_t3(delta):
      return np.exp(-beta*k_a*(delta)**2)
      

    ### Integrate translation and rotation
    r1_int,a1_int,t1_int,a2_int,t2_int,t3_int = [0.0,0.0,0.0,0.0,0.0,0.0]
    intrange = np.arange(r1lb,r1ub,r1st)
    r1_int = np.trapz(f_r1(intrange),intrange)
    intrange = np.arange(a1lb,a1ub,a1st)
    a1_int = np.trapz(f_a1(intrange),intrange)
    intrange = dih_per(t1lb,t1ub,t1st,t1_0)
    t1_int = np.trapz(f_t1(intrange),intrange)
    intrange = np.arange(a2lb,a2ub,a2st)
    a2_int = np.trapz(f_a2(intrange),intrange)
    intrange = dih_per(t2lb,t2ub,t2st,t2_0)
    t2_int = np.trapz(f_t2(intrange),intrange)
    intrange = dih_per(t3lb,t3ub,t3st,t3_0)
    t3_int = np.trapz(f_t3(intrange),intrange)
    return R*temperature*np.log((1/(8.0*np.pi*np.pi))*(1.0/1660.0)*r1_int*a1_int*t1_int*a2_int*t2_int*t3_int)

def fe_int_op(r1_0, a1_0, t1_0, a2_0, t2_0, t3_0, k_r, k_a, temperature):

    R = 1.987204118e-3 # kcal/mol-K, a.k.a. boltzman constant
    beta = 1/(temperature*R)

    # Numerical integration limits and spacing (trapezoid)
    r1lb,r1ub,r1st = [0.0,100.0,0.00001]
    a1lb,a1ub,a1st = [0.0,np.pi,0.000001]
    t1lb,t1ub,t1st = [-np.pi,np.pi,0.000001]
    a2lb,a2ub,a2st = [0.0,np.pi,0.000001]

    # Potential energy expressions
    def f_r1(val):
      return (val**2)*np.exp(-beta*k_r*(val-r1_0)**2)
    def f_a1(val):
      return np.sin(val)*np.exp(-beta*k_a*(val-np.radians(a1_0))**2)
    def f_a2(val):
      return np.sin(val)*np.exp(-beta*k_a*(val-np.radians(a2_0))**2)
    def f_t1(val):
      return np.exp(-beta*k_a*(1+np.cos(val-np.radians(t1_0)-np.pi))*2)

     
    # Integrate translation and rotation
    r1_int,a1_int,t1_int,a2_int = [0.0,0.0,0.0,0.0]
    intrange = np.arange(r1lb,r1ub,r1st)
    r1_int = np.trapz(f_r1(intrange),intrange)
    intrange = np.arange(a1lb,a1ub,a1st)
    a1_int = np.trapz(f_a1(intrange),intrange)
    intrange = np.arange(t1lb,t1ub,t1st)
    t1_int = np.trapz(f_t1(intrange),intrange)
    intrange = np.arange(a2lb,a2ub,a2st)
    a2_int = np.trapz(f_a2(intrange),intrange)

    # Output total TR release free energy 
    return R*temperature*np.log((1/(8.0*np.pi*np.pi))*(1.0/1660.0)*r1_int*a1_int*t1_int*a2_int*t1_int*t1_int)

def fe_dd(comp, pose, mode, lambdas, weights, dec_int, dec_method, rest_file, temperature):

    kB = 1.381e-23 * 6.022e23 / (4.184 * 1000.0) # Boltzmann constant in kJ/mol/K
    beta = 1/(kB * temperature) # beta
    N_max = 20000 # Max frames for any simulation window, you should check this if you did some long runs


    os.chdir('fe')
    os.chdir(pose)
    os.chdir(dec_method)
    if not os.path.exists('data'):
      os.makedirs('data')

    # Define log file
    sys.stdout = open('./data/'+dec_int+'-'+comp+'-'+mode+'.dat', 'w')

    ### Determine Number of windows
    K = 0
    filename = './'+comp+'%02.0f/%s' % (K, rest_file)
    while os.path.isfile(filename):
      K = K+1
      filename = './'+comp+'%02.0f/%s' % (K, rest_file)

    if dec_int == 'ti':
      deltag = 0
      dvdl = []
      for k in range(K):
        data = []
        # Read in Values for restrained variables for each simulation
        filename = './'+comp+'%02.0f/%s' % (k, rest_file)
        infile = open(filename, 'r')
        restdat = infile.readlines()     # slice off first 20 lines  readlines()[20:]
        infile.close()
        # Parse Data
        for line in restdat:
          data.append(float(line.split()[1]))
        dvdl.append(float(sum(data)/len(data)))

      for i in range(0, len(dvdl)):
        print ('%-10s%6.5f,  %-8s%9.5f' % ('lambda =', float(lambdas[i]), 'dvdl =', float(dvdl[i])))

      for i in range(K):
        deltag = deltag + dvdl[i]*weights[i]

      print ('\n%-8s %9.5f' % ('deltaG  ', float(deltag)))
    elif dec_int == 'mbar':

      ### Allocate storage for simulation data
      N = np.zeros([K], np.int32)                       # N_k[k] is the number of snapshots to be used from umbrella simulation k
      Neff = np.zeros([K], np.int32)
      Nind = np.zeros([K], np.int32)
      val = np.zeros([N_max,K,K], np.float64)           # value of the restrained variable at each frame n
      g = np.zeros([K], np.float64)
      u=np.zeros([N_max], np.float64)

      ### Calculate Statistical Inefficiency (g)
      def calcg(data):
        sum = 0
        randnum = ("%05.0f" % (int(100000*np.random.random())))
        datafn = '/dev/shm/series.'+randnum+'.dat'
        acffn = '/dev/shm/acf.'+randnum+'.dat'
        cppfn = '/dev/shm/pt-acf.'+randnum+'.in'
        np.savetxt(datafn,data)
        cpptin = open(cppfn, 'w')
        cpptin.write("readdata "+datafn+" name "+randnum+"\nautocorr "+randnum+" out "+acffn+" noheader\n")
        cpptin.close()

        FNULL = open(os.devnull, 'w')
        sp.call(['cpptraj','-i',cppfn], stdout=FNULL, stderr=sp.STDOUT)

        with open(acffn, 'r') as acf:
          for line in acf:
            col = line.split()
            t = float(col[0]) - 1.0
        T = t

        with open(acffn, 'r') as acf:
          for line in acf:
            col = line.split()
            t = float(col[0]) - 1.0
            v = float(col[1])
            if t == 0:
              continue
            if v < 0.0:
              break
            sum += ( 1 - (t/T) )*(v)

        sp.call(['rm',datafn,acffn,cppfn])

        return 1+(2*sum)

      for k in range(K):
        # Read in Values for restrained variables for each simulation
        filename = './'+comp+'%02.0f/%s' % (k, rest_file)
        infile = open(filename, 'r')
        restdat = infile.readlines()     # slice off first 20 lines  readlines()[20:]
        infile.close()
        # Parse Data
        n = 0
        lambdas = []
        for line in restdat:         
          cols = line.split()
          if len(cols) >= 1:
            lambdas.append(float(cols[1])) 
          if len(cols) == 0:
            break
        for line in restdat:         
          cols = line.split()
          if len(cols) >= 1:
            if '**' not in cols[2]:
              lamb = float(cols[1].strip())
              val[n,k,lambdas.index(lamb)] = cols[2]
          if len(cols) == 0:
            n += 1
        N[k] = n
        
        # Calculate reduced potential
        u[0:N[k]] = beta*(val[0:N[k],k,k])

        # Subsample or not
        if mode == 'sub':
          g[k] = calcg(u[0:N[k]])
          subs = timeseries.subsampleCorrelatedData(np.zeros([N[k]]),g=g[k])
          Nind[k] = len(subs)
          Neff[k] = Nind[k]
        else:
          g[k] = 1.00
          Neff[k] = N[k]

        print  ("Processed Window %5.0f.  N= %12.0f.  g= %10.3f   Neff= %12.0f" % ( k, N[k], g[k], Neff[k] ))


      # Calculate decoupling energy
      Upot = np.zeros([K,K,np.max(Neff)], np.float64)
      for k in range(K):
        for l in range(K): 
          Upot[k,l,0:Neff[k]] = beta*(val[0:Neff[k],k,l])

      val = []

      print  ("\nRunning MBAR... ")
      mbar = MBAR(Upot, Neff)

      print  ("Calculate Free Energy Differences Between States")
      [Deltaf, dDeltaf] = mbar.getFreeEnergyDifferences()

      min = np.argmin(Deltaf[0])

      # Write to file
      print  ("\nFree Energy Differences (in units of kcal/mol)")
      print  ("%9s %8s %8s" % ('lambda', 'f', 'df'))
      for k in range(K):
        print ("%10.5f %10.5f %10.5f" % ( lambdas[k], Deltaf[0,k]/beta, dDeltaf[0,k]/beta))
      print ("\n\n")

    os.chdir('../../../')


