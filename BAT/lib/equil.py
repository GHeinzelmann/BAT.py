#!/usr/bin/env python
# coding: utf-8
import math
import openmmtools
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
from sys import stdout
import numpy as np
from importlib import reload
from openmmtools.multistate import MultiStateReporter, MultiStateSampler, ReplicaExchangeSampler, ParallelTemperingSampler, SAMSSampler, ReplicaExchangeAnalyzer
from openmmtools.states import GlobalParameterState
import os

kB = 0.008314472471220214 * unit_definitions.kilojoules_per_mole/unit_definitions.kelvin
temp = TMPRT*unit_definitions.kelvin
kT = kB * temp
kcal = 4.1868 * unit_definitions.kilojoules_per_mole
kTtokcal = kT/kcal * unit_definitions.kilocalories_per_mole
runflag='run'
stge='STG'
steps=TOTST
out_freq=int(steps/10)

#State with a lambda controlling the restraints
class RestraintComposableState(GlobalParameterState):
	lambda_restraints=GlobalParameterState.GlobalParameter('lambda_restraints', standard_value=1.0)

args = sys.argv[1:]
if len(args) > 1:
	raise ValueError('Only take one argument runflag')

if len(args) == 0:
	runflag='run'
else:
	runflag=args[0]
	allowedrunflags = ['run', 'extend', 'recover']
	if runflag not in allowedrunflags:
		raise ValueError('Please select runflag from {}'.format(allowedrunflags))

# Generate system from AMBER
prmtop = AmberPrmtopFile('PRMFL')
inpcrd = AmberInpcrdFile('full.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=CTF*nanometer, constraints=HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)

################# Read AMBER restraint file and set ligand decoupling atoms #############

atoms_list = []
eq_values = []
fcn_values = []
eq_tr = []
eq_pc = []
eq_lc = []
fcn_tr = []
fcn_pc = []
fcn_lc = []
fcn_com = []

# Protein conf restraints
# Open AMBER restraint file for reading 
with open('disang'+stge+'.rest', "r") as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if '#Rec_C' in line or '#Rec_D' in line)
    for i in range(0, len(lines)):
      splitdata = lines[i].split()
      atoms = splitdata[2].split(',')[0:-1]
      eql = float(splitdata[6].strip(','))
      fcn = float(splitdata[12].strip(','))
      atoms_list.append(atoms)
      eq_values.append(eql)
      fcn_values.append(fcn)
for i in range(0, len(atoms_list)):
    # Adjust atom numbering 
    for j in range(0, len(atoms_list[i])):
      atoms_list[i][j]=int(atoms_list[i][j])-1 
    if len(atoms_list[i]) == 2:
      eq_pc.append(eq_values[i])
      fcn_pc.append(fcn_values[i]*2)
    elif len(atoms_list[i]) == 3:
      eq_pc.append(eq_values[i]*np.pi/180)
      fcn_pc.append(fcn_values[i]*2)
    elif len(atoms_list[i]) == 4:
      eq_pc.append((eq_values[i]*np.pi/180)-np.pi)
      fcn_pc.append(fcn_values[i]*2)
atoms_pc=list(atoms_list)

atoms_list = []
eq_values = []
fcn_values = []

# Ligand TR restraints
# Open AMBER restraint file for reading 
with open('disang'+stge+'.rest', "r") as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if '#Lig_TR' in line)
    for i in range(0, len(lines)):
      splitdata = lines[i].split()
      atoms = splitdata[2].split(',')[0:-1]
      eql = float(splitdata[6].strip(','))
      fcn = float(splitdata[12].strip(','))
      atoms_list.append(atoms)
      eq_values.append(eql)
      fcn_values.append(fcn)
for i in range(0, len(atoms_list)):
    # Adjust atom numbering 
    for j in range(0, len(atoms_list[i])):
      atoms_list[i][j]=int(atoms_list[i][j])-1 
    if len(atoms_list[i]) == 2:
      eq_tr.append(eq_values[i])
      fcn_tr.append(fcn_values[i]*2)
    elif len(atoms_list[i]) == 3:
      eq_tr.append(eq_values[i]*np.pi/180)
      fcn_tr.append(fcn_values[i]*2)
    elif len(atoms_list[i]) == 4:
      eq_tr.append((eq_values[i]*np.pi/180)-np.pi)
      fcn_tr.append(fcn_values[i]*2)
atoms_tr=list(atoms_list)

atoms_list = []
eq_values = []
fcn_values = []

# Ligand conf restraints
# Open AMBER restraint file for reading 
with open('disang'+stge+'.rest', "r") as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if '#Lig_C' in line or '#Lig_D' in line)
    for i in range(0, len(lines)):
      splitdata = lines[i].split()
      atoms = splitdata[2].split(',')[0:-1]
      eql = float(splitdata[6].strip(','))
      fcn = float(splitdata[12].strip(','))
      atoms_list.append(atoms)
      eq_values.append(eql)
      fcn_values.append(fcn)
for i in range(0, len(atoms_list)):
    # Adjust atom numbering 
    for j in range(0, len(atoms_list[i])):
      atoms_list[i][j]=int(atoms_list[i][j])-1 
    if len(atoms_list[i]) == 2:
      eq_lc.append(eq_values[i])
      fcn_lc.append(fcn_values[i]*2)
    elif len(atoms_list[i]) == 3:
      eq_lc.append(eq_values[i]*np.pi/180)
      fcn_lc.append(fcn_values[i]*2)
    elif len(atoms_list[i]) == 4:
      eq_lc.append((eq_values[i]*np.pi/180)-np.pi)
      fcn_lc.append(fcn_values[i]*2)
atoms_lc=list(atoms_list)

atoms_list = []
eq_values = []
fcn_values = []

# Center of mass restraints

# Read atoms from pdb file
with open('full.pdb') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

dum_atom = 0

# Count dummy atoms
for i in range(0, len(lines)):
  if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
    if lines[i][17:20].strip() == 'DUM':
      dum_atom += 1

# Read dummy coordinates 
with open('full.inpcrd') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

splitdata = lines[2].split()
if dum_atom >= 1:    
  xr = splitdata[0]
  yr = splitdata[1]
  zr = splitdata[2]
if dum_atom > 1:    
  xb = splitdata[3]
  yb = splitdata[4]
  zb = splitdata[5]

# Open AMBER collective variable file for reading 
with open('cv.in', "r") as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line)
    for i in range(0, len(lines)):
      if 'cv_ni' in lines[i]:
        splitdata = lines[i].split()
        atoms = splitdata[5].split(',')[2:-1]
        atoms_list.append(atoms)
      if 'anchor_strength' in lines[i]:
        splitdata = lines[i].split()
        fcn = float(splitdata[2].strip(','))
        fcn_values.append(fcn)

for i in range(0, len(atoms_list)):
    # Adjust atom numbering 
    for j in range(0, len(atoms_list[i])):
      atoms_list[i][j]=int(atoms_list[i][j])-1 
    fcn_com.append(fcn_values[i])

com_atoms=list(atoms_list)

atoms_list = []
eq_values = []

# Get ligand atoms 

# Read atoms list from pdb file
with open('full.pdb') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

# Create ligand list
for i in range(0, len(lines)):
    if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
      if lines[i][17:20].strip() == 'LIG':
        atoms_list.append(lines[i][6:11].strip())

for i in range(0, len(atoms_list)):
    # Adjust atom numbering 
    atoms_list[i]=int(atoms_list[i])-1 

dec_atoms=list(atoms_list)

print('')
print('Ligand atoms')
print('')
print(dec_atoms)
print('')

#########################################################################################


#Pressure control
system.addForce(MonteCarloBarostat(1*unit_definitions.atmospheres, TMPRT*unit_definitions.kelvin, 25))

restraint_state = RestraintComposableState(lambda_restraints=1.0)

#Apply harmonic restraints

# Set bond restraint properties 
bond_energy_function = "lambda_restraints*(K/2)*(r-r0)^2;"
harmonicforce=CustomBondForce(bond_energy_function)
harmonicforce.addPerBondParameter('r0')
harmonicforce.addPerBondParameter('K')
harmonicforce.addGlobalParameter('lambda_restraints', 1.0)
new_force_index=system.getNumForces()
harmonicforce.setForceGroup(new_force_index)

# Set angle restraint properties 
angle_energy_function = "lambda_restraints*0.5*k*(theta-theta0)^2"
angleforce=CustomAngleForce(angle_energy_function)
angleforce.addPerAngleParameter('theta0')
angleforce.addPerAngleParameter('k')
angleforce.addGlobalParameter('lambda_restraints', 1.0)
angleforce.setForceGroup(new_force_index)

# Set torsion restraint properties 
torsion_energy_function = "lambda_restraints*k*(1+cos(n*theta-theta0))"
torsionforce=CustomTorsionForce(torsion_energy_function)
torsionforce.addPerTorsionParameter('n')
torsionforce.addPerTorsionParameter('theta0')
torsionforce.addPerTorsionParameter('k')
torsionforce.addGlobalParameter('lambda_restraints', 1.0)
torsionforce.setForceGroup(new_force_index)

# Protein conf restraints
print('Protein conformational restraints:')
print('')
for i in range(0, len(atoms_pc)):
    if len(atoms_pc[i]) == 2:
      bondlength=eq_pc[i]*unit_definitions.angstroms
      bondforce=fcn_pc[i]*unit_definitions.kilocalorie_per_mole/unit_definitions.angstroms**2
      harmonicforce.addBond(atoms_pc[i][0],atoms_pc[i][1], [bondlength,bondforce])
      print('harmonicforce.addBond('+str(atoms_pc[i][0])+','+str(atoms_pc[i][1])+', ['+str(bondlength)+','+str(bondforce)+']')
    elif len(atoms_pc[i]) == 4:
      torsionconst=fcn_pc[i]*unit_definitions.kilocalorie_per_mole
      torsionforce.addTorsion(atoms_pc[i][0],atoms_pc[i][1],atoms_pc[i][2],atoms_pc[i][3], [1, eq_pc[i]*unit_definitions.radian,torsionconst])
      print('torsionforce.addTorsion('+str(atoms_pc[i][0])+','+str(atoms_pc[i][1])+','+str(atoms_pc[i][2])+','+str(atoms_pc[i][3])+', [1, '+str(eq_pc[i])+','+str(torsionconst)+'])')
print('')

# TR restraints
print('Ligand TR restraints:')
print('')
for i in range(0, len(atoms_tr)):
    if len(atoms_tr[i]) == 2:
      bondlength=eq_tr[i]*unit_definitions.angstroms
      bondforce=fcn_tr[i]*unit_definitions.kilocalorie_per_mole/unit_definitions.angstroms**2
      harmonicforce.addBond(atoms_tr[i][0],atoms_tr[i][1], [bondlength,bondforce])
      print('harmonicforce.addBond('+str(atoms_tr[i][0])+','+str(atoms_tr[i][1])+', ['+str(bondlength)+','+str(bondforce)+']')
    elif len(atoms_tr[i]) == 3:
      restrainedangle=eq_tr[i]*unit_definitions.radians
      angleconst=fcn_tr[i]*unit_definitions.kilocalorie_per_mole/unit_definitions.radian**2
      angleforce.addAngle(atoms_tr[i][0],atoms_tr[i][1],atoms_tr[i][2], [restrainedangle,angleconst])
      print('angleforce.addAngle('+str(atoms_tr[i][0])+','+str(atoms_tr[i][1])+','+str(atoms_tr[i][2])+', ['+str(restrainedangle)+','+str(angleconst)+']')
    elif len(atoms_tr[i]) == 4:
      torsionconst=fcn_tr[i]*unit_definitions.kilocalorie_per_mole
      torsionforce.addTorsion(atoms_tr[i][0],atoms_tr[i][1],atoms_tr[i][2],atoms_tr[i][3], [1, eq_tr[i]*unit_definitions.radian,torsionconst])
      print('torsionforce.addTorsion('+str(atoms_tr[i][0])+','+str(atoms_tr[i][1])+','+str(atoms_tr[i][2])+','+str(atoms_tr[i][3])+', [1, '+str(eq_tr[i])+','+str(torsionconst)+'])')
print('')

# Ligand conf restraints
print('Ligand conformational restraints:')
print('')
for i in range(0, len(atoms_lc)):
    if len(atoms_lc[i]) == 4:
      torsionconst=fcn_lc[i]*unit_definitions.kilocalorie_per_mole
      torsionforce.addTorsion(atoms_lc[i][0],atoms_lc[i][1],atoms_lc[i][2],atoms_lc[i][3], [1, eq_lc[i]*unit_definitions.radian,torsionconst])
      print('torsionforce.addTorsion('+str(atoms_lc[i][0])+','+str(atoms_lc[i][1])+','+str(atoms_lc[i][2])+','+str(atoms_lc[i][3])+', [1, '+str(eq_lc[i])+','+str(torsionconst)+'])')
print('')


# Center of Mass (COM) restraints for receptor and the bulk ligand 

bondGroups = []
bondParameters = []

if dum_atom >= 1:
  # Configure COM expression
  comforce=CustomCentroidBondForce(1, "0.5*k*((x1-x0)^2+(y1-y0)^2+(z1-z0)^2)")
  comforce.addPerBondParameter('k')
  comforce.addPerBondParameter('x0')
  comforce.addPerBondParameter('y0')
  comforce.addPerBondParameter('z0')
  comforce.setForceGroup(new_force_index)

  # Add receptor COM restraints
  print('Receptor COM restraints:')
  print('')
  comforce.addGroup(com_atoms[0])
  bondGroups.append(0)
  bondParameters.append(float(fcn_com[0])*unit_definitions.kilocalorie_per_mole/unit_definitions.angstroms**2)
  bondParameters.append(float(xr)*unit_definitions.angstroms)
  bondParameters.append(float(yr)*unit_definitions.angstroms)
  bondParameters.append(float(zr)*unit_definitions.angstroms)
  comforce.addBond(bondGroups, bondParameters)
  print('comforce.addBond('+str(com_atoms[0])+', '+str(bondParameters[0])+', '+str(bondParameters[1])+', '+str(bondParameters[2])+', '+str(bondParameters[3])+')')
  print('')

  bondGroups = []
  bondParameters = []
  dummy_group = [ 0 ]
  dum_fcn = 20.0

  # Restrain COM dummy to use as reference later
  print('Dummy COM restraints:')
  print('')
  comforce.addGroup(dummy_group)
  bondGroups.append(1)
  bondParameters.append(float(dum_fcn)*unit_definitions.kilocalorie_per_mole/unit_definitions.angstroms**2)
  bondParameters.append(float(xr)*unit_definitions.angstroms)
  bondParameters.append(float(yr)*unit_definitions.angstroms)
  bondParameters.append(float(zr)*unit_definitions.angstroms)
  comforce.addBond(bondGroups, bondParameters)
  print('comforce.addBond('+str(dummy_group)+', '+str(bondParameters[0])+', '+str(bondParameters[1])+', '+str(bondParameters[2])+', '+str(bondParameters[3])+')')
  print('')
 
system.addForce(harmonicforce) # after
system.addForce(angleforce) # after
system.addForce(torsionforce) # after 
system.addForce(comforce) # after 

#setup integrator

integrator=LangevinIntegrator(TMPRT*unit_definitions.kelvin, GAMMA_LN/unit_definitions.picoseconds, TSTP*unit_definitions.femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator)

if stge != '00':
  prv = int(stge) - 1
  with open('md%02d.chk' %prv, 'rb') as f:
    simulation.context.loadCheckpoint(f.read())
else:
  simulation.context.setPositions(inpcrd.positions)
  if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  # Minimize when starting
  simulation.minimizeEnergy()
  simulation.context.setVelocitiesToTemperature(TMPRT*unit_definitions.kelvin)

# Run simulation

simulation.reporters.append(DCDReporter('md'+stge+'.dcd', out_freq))

simulation.reporters.append(CheckpointReporter('md'+stge+'.chk', out_freq))

simulation.reporters.append(StateDataReporter(stdout, out_freq, step=True,

        potentialEnergy=True, temperature=True, progress=True, remainingTime=True,

        speed=True, totalSteps=steps, separator='      '))

simulation.step(steps)

