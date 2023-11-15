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
comp='CMPN'

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

#Dynamics
timeStep=TSTP*unit_definitions.femtoseconds
stepsPerIteration=SPITR
productionIterations=PRIT
equilibrationIterations=EQIT
iterationsPerCheckpoint=ITCH
extendIterations=1000


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
with open('disang.rest', "r") as f_in:
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
with open('disang.rest', "r") as f_in:
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
with open('disang.rest', "r") as f_in:
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

if dum_atom >= 1:    
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

# Get ligand atoms for decoupling

# Read atoms list from pdb file
with open('full.pdb') as f_in:
    lines = (line.rstrip() for line in f_in)
    lines = list(line for line in lines if line) # Non-blank lines in a list   

# Create decoupling list
for i in range(0, len(lines)):
    if (lines[i][0:6].strip() == 'ATOM') or (lines[i][0:6].strip() == 'HETATM'):
      if lines[i][17:20].strip() == 'LIG':
        atoms_list.append(lines[i][6:11].strip())

for i in range(0, len(atoms_list)):
    # Adjust atom numbering 
    atoms_list[i]=int(atoms_list[i])-1 

dec_atoms=list(atoms_list)

print('')
print('Ligand decoupling atoms')
print('')
print(dec_atoms)
print('')

#########################################################################################


#Pressure control
system.addForce(MonteCarloBarostat(1*unit_definitions.atmospheres, TMPRT*unit_definitions.kelvin, 25))

restraint_state = RestraintComposableState(lambda_restraints=1.0)
ligand_a_atoms = list(dec_atoms)
ligand_a_bonds = []
ligand_a_angles =[]
ligand_a_torsions =[]

if comp == 'r':

  # No ligand restraints
  atoms_lc = []
  atoms_tr = []

if comp == 'c':

  # Only ligand conf restraints
  atoms_pc = []
  atoms_tr = []

#Fixed harmonic restraints

if comp == 'l' or comp == 't':

  # Set bond restraint properties 
  bond_energy_function = "(K/2)*(r-r0)^2;"
  harmonicforce=CustomBondForce(bond_energy_function)
  harmonicforce.addPerBondParameter('r0')
  harmonicforce.addPerBondParameter('K')
  new_force_index=system.getNumForces()
  harmonicforce.setForceGroup(new_force_index)

  # Set angle restraint properties 
  angle_energy_function = "0.5*k*(theta-theta0)^2"
  angleforce=CustomAngleForce(angle_energy_function)
  angleforce.addPerAngleParameter('theta0')
  angleforce.addPerAngleParameter('k')
  angleforce.setForceGroup(new_force_index)

  # Set torsion restraint properties 
  torsion_energy_function = "k*(1+cos(n*theta-theta0))"
  torsionforce=CustomTorsionForce(torsion_energy_function)
  torsionforce.addPerTorsionParameter('n')
  torsionforce.addPerTorsionParameter('theta0')
  torsionforce.addPerTorsionParameter('k')
  torsionforce.setForceGroup(new_force_index)


  # Fixed Protein conf restraints
  print('Fixed protein conformational restraints:')
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

  if comp == 't':

    # Ligand conf restraints
    print('Fixed ligand conformational restraints:')
    print('')
    for i in range(0, len(atoms_lc)):
      if len(atoms_lc[i]) == 4:
        torsionconst=fcn_lc[i]*unit_definitions.kilocalorie_per_mole
        torsionforce.addTorsion(atoms_lc[i][0],atoms_lc[i][1],atoms_lc[i][2],atoms_lc[i][3], [1, eq_lc[i]*unit_definitions.radian,torsionconst])
        print('torsionforce.addTorsion('+str(atoms_lc[i][0])+','+str(atoms_lc[i][1])+','+str(atoms_lc[i][2])+','+str(atoms_lc[i][3])+', [1, '+str(eq_lc[i])+','+str(torsionconst)+'])')
    print('')

  system.addForce(harmonicforce) # after
  system.addForce(angleforce) # after
  system.addForce(torsionforce) # after 

#Lambda-dependent harmonic restraints

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

if comp == 'a' or comp == 'r' or comp == 'm' or comp == 'n':

  # Lambda dependent Protein conf restraints
  print('Lambda dependent protein conformational restraints:')
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

if comp == 'l' or comp == 'c' or comp == 'm' or comp == 'n':

  # Lambda dependent Ligand conf restraints
  print('Lambda dependent ligand conformational restraints:')
  print('')
  for i in range(0, len(atoms_lc)):
    if len(atoms_lc[i]) == 4:
      torsionconst=fcn_lc[i]*unit_definitions.kilocalorie_per_mole
      torsionforce.addTorsion(atoms_lc[i][0],atoms_lc[i][1],atoms_lc[i][2],atoms_lc[i][3], [1, eq_lc[i]*unit_definitions.radian,torsionconst])
      print('torsionforce.addTorsion('+str(atoms_lc[i][0])+','+str(atoms_lc[i][1])+','+str(atoms_lc[i][2])+','+str(atoms_lc[i][3])+', [1, '+str(eq_lc[i])+','+str(torsionconst)+'])')
  print('')

if comp == 't' or comp == 'm': 

  # TR restraints
  print('Lambda dependent ligand TR restraints:')
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

if dum_atom > 1:

  bondGroups = []
  bondParameters = []

  # Add bulk ligand COM restraints
  print('Bulk ligand COM restraints:')
  print('')
  comforce.addGroup(com_atoms[1])
  bondGroups.append(1)
  bondParameters.append(float(fcn_com[1])*unit_definitions.kilocalorie_per_mole/unit_definitions.angstroms**2)
  bondParameters.append(float(xb)*unit_definitions.angstroms)
  bondParameters.append(float(yb)*unit_definitions.angstroms)
  bondParameters.append(float(zb)*unit_definitions.angstroms)
  comforce.addBond(bondGroups, bondParameters)
  print('comforce.addBond('+str(com_atoms[1])+', '+str(bondParameters[0])+', '+str(bondParameters[1])+', '+str(bondParameters[2])+', '+str(bondParameters[3])+')')
  print('')
 
system.addForce(harmonicforce) # after
system.addForce(angleforce) # after
system.addForce(torsionforce) # after 
if dum_atom >= 1:
  system.addForce(comforce) # after 

#Get date from forcefield

num_a_atoms=len(ligand_a_atoms)
num_a_bonds=len(ligand_a_bonds)
num_a_angles=len(ligand_a_angles)
num_a_torsions=len(ligand_a_torsions)

#Setup alchemical system
reload(openmmtools.alchemy)
factory = openmmtools.alchemy.AbsoluteAlchemicalFactory(consistent_exceptions=False, split_alchemical_forces = True, alchemical_pme_treatment = 'exact') # RIZZI CHECK
reference_system = system

#OpenMM crashes if one adds nonexistent alchemical angles or torsions
if(num_a_bonds == 0):
	if(num_a_angles == 0):
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms,alchemical_torsions = ligand_a_torsions, name='A')
	else:
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_angles = ligand_a_angles, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_angles = ligand_a_angles, alchemical_torsions = ligand_a_torsions, name='A')
else:
	if(num_a_angles == 0):
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_torsions = ligand_a_torsions, name='A')
	else:
		if(num_a_torsions == 0):
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_angles = ligand_a_angles, name='A')
		else:
			alchemical_region_A = openmmtools.alchemy.AlchemicalRegion(alchemical_atoms = ligand_a_atoms, alchemical_bonds=ligand_a_bonds, alchemical_angles = ligand_a_angles, alchemical_torsions = ligand_a_torsions, name='A')

alchemical_system_in = factory.create_alchemical_system(reference_system, alchemical_regions = [alchemical_region_A])

alchemical_state_A = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'A')
reload(openmmtools.alchemy)
TS = openmmtools.states.ThermodynamicState(alchemical_system_in, temperature=TMPRT*unit_definitions.kelvin, pressure=1*unit_definitions.bar)
composable_states = [alchemical_state_A, restraint_state]
compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
reload(openmmtools.alchemy)
integrator=LangevinIntegrator(TMPRT*unit_definitions.kelvin, GAMMA_LN/unit_definitions.picoseconds, TSTP*unit_definitions.femtoseconds)
context = compound_state.create_context(integrator)
alchemical_system_in=context.getSystem()

#Use offsets to interpolate
alchemical_state_A = openmmtools.alchemy.AlchemicalState.from_system(alchemical_system_in, parameters_name_suffix = 'A')
reload(openmmtools.alchemy)
TS = openmmtools.states.ThermodynamicState(alchemical_system_in, temperature=TMPRT*unit_definitions.kelvin, pressure=1*unit_definitions.bar)
composable_states = [alchemical_state_A, restraint_state]
compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
reload(openmmtools.alchemy)

#DEBUG info
sys = compound_state.get_system()
file = open('DEBUG_restraint.xml','w')
file.write(XmlSerializer.serialize(sys))
file.close()

#Setup lambda values
lambdas = LAMBDAS

nstates=len(lambdas)
print("There will be ", nstates, " states in total")
print("stepsPerIteration:", stepsPerIteration, " productionIterations: ", productionIterations, "equilibrationIterations: ", equilibrationIterations)
print("Timestep: ", timeStep)
box_vec = alchemical_system_in.getDefaultPeriodicBoxVectors()
print("Box vectors:", box_vec)

#Sanity check
print("")
print("Lambdas matrix")
print("Lrestraints")
for j in range(len(lambdas)):
    print("%-15.7f" % lambdas[j])

sampler_states = list()
thermodynamic_states = list()

for k in range(nstates):
    compound_state = openmmtools.states.CompoundThermodynamicState(thermodynamic_state=TS, composable_states=composable_states)
    if comp == 'r':
      compound_state.lambda_sterics_A=0.0
      compound_state.lambda_electrostatics_A=0.0
    else:
      compound_state.lambda_sterics_A=1.0
      compound_state.lambda_electrostatics_A=1.0
    compound_state.lambda_restraints=lambdas[k]
    sys = compound_state.get_system()
    sampler_states.append(openmmtools.states.SamplerState(positions=inpcrd.positions, box_vectors=box_vec))
    thermodynamic_states.append(compound_state)

print("Integrator: LangevinSplittingDynamicsMove")
print("Sampler: ReplicaExchangeSampler")

lsd_move = openmmtools.mcmc.LangevinSplittingDynamicsMove(timestep=timeStep, collision_rate=GAMMA_LN/unit_definitions.picoseconds, n_steps=stepsPerIteration)
print('Minimizing......')
for k in range(nstates):
      sampler_state = sampler_states[k]
      thermodynamic_state = thermodynamic_states[k]
      integrator=LangevinIntegrator(TMPRT*unit_definitions.kelvin, GAMMA_LN/unit_definitions.picoseconds, TSTP*unit_definitions.femtoseconds)
      context = thermodynamic_state.create_context(integrator)
      system = context.getSystem()
      for force in system.getForces(): # RIZZI CHECK
              if isinstance(force, CustomBondForce):
                      force.updateParametersInContext(context)
              elif isinstance(force, HarmonicBondForce):
                      force.updateParametersInContext(context)
              elif isinstance(force, HarmonicAngleForce):
                      force.updateParametersInContext(context)
              elif isinstance(force, PeriodicTorsionForce):
                      force.updateParametersInContext(context)
              elif isinstance(force, CustomAngleForce):
                      force.updateParametersInContext(context)
              elif isinstance(force, NonbondedForce):
                      force.updateParametersInContext(context)
              elif isinstance(force, CustomNonbondedForce):
                      force.updateParametersInContext(context)
              elif isinstance(force, CustomTorsionForce):
                      force.updateParametersInContext(context)
      sampler_state.apply_to_context(context)
      initial_energy = thermodynamic_state.reduced_potential(context)
      print("Sampler state {}: initial energy {:8.3f}kT".format(k, initial_energy))
      LocalEnergyMinimizer.minimize(context)
      sampler_state.update_from_context(context)
      final_energy = thermodynamic_state.reduced_potential(context)
      print("Sampler state {}: final energy {:8.3f}kT".format(k, final_energy))
      del context
print('Minimized......')

if runflag == 'run':
      repex_simulation = ReplicaExchangeSampler(mcmc_moves=lsd_move, number_of_iterations=productionIterations)

tmp_dir = './trajectory/'
storage = os.path.join(tmp_dir, 'restraint.nc')
reporter = MultiStateReporter(storage, checkpoint_interval=iterationsPerCheckpoint)
if runflag != 'run':
    repex_simulation = ReplicaExchangeSampler.from_storage(reporter)
else:
    repex_simulation.create(thermodynamic_states, sampler_states, reporter)
    print('Equilibrating......')
    repex_simulation.equilibrate(equilibrationIterations)
    print('Simulating......')
if runflag == 'recover' or runflag == 'run':
        repex_simulation.run()
elif runflag == 'extend':
        repex_simulation.extend(extendIterations)

#will add all iterations even if coming from a previous restart
all_iters = repex_simulation.iteration
print('All iterations = {}'.format(all_iters))

final_en = 0
blocks = BLCKS
print('Blocks to analyze = {}'.format(blocks))
for i in range(0, blocks):
    init_n=int(round(i*(all_iters/blocks)))
    max_n=int(round((i+1)*(all_iters/blocks)))
    print(init_n)
    print(max_n)
    analyzer = ReplicaExchangeAnalyzer(reporter, n_equilibration_iterations=init_n)
    analyzer.max_n_iterations = max_n
    Delta_f_ij, dDelta_f_ij = analyzer.get_free_energy()
    print("Relative free energy change for block {0} = {1} +- {2}".format(i+1, Delta_f_ij[0, nstates - 1]*kTtokcal, dDelta_f_ij[0, nstates - 1]*kTtokcal))
    block_en = Delta_f_ij[0, nstates - 1]*kT/kcal
    print(block_en)
    final_en = final_en + block_en/blocks

print(final_en)

analyzer = ReplicaExchangeAnalyzer(reporter)
analyzer.max_n_iterations = all_iters
Delta_f_ij, dDelta_f_ij = analyzer.get_free_energy()
print("Relative free energy change for the whole {0} = {1} +- {2}".format('decoupling/recoupling', Delta_f_ij[0, nstates - 1]*kTtokcal, dDelta_f_ij[0, nstates - 1]*kTtokcal))

[matrix,eigenvalues,ineff]=analyzer.generate_mixing_statistics()
print("Mixing Stats")
print(matrix)
print(eigenvalues)
print(ineff)
