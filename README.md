*Note: The master branch is for the 2.0 version of BAT. For the version associated with Ref. [1], which includes the APR method, download the package from the BATv1.0 branch.* 

--------------------------------------------------------

*New features on BAT 2.0:*

 
*- Relative restraints between the receptor and the ligand (figure below), without the need for three fixed dummy atoms.*

*- Center of mass restraints on the receptor, and the bulk ligand when the SDR method is applied, so their internal degrees of freedom are not affected.*

*- Only two stages, equilibrium and free energy simulations. The preparation stage is no longer required, since the APR method is not available in the 2.0 version.*

*- Possible choice between a fixed number of waters, or fixed solvation buffers in the three cartesian axes.*

*- Automatic determination of the number of ions based on the chosen salt concentration.*

*- Simpler procedure to add new receptors.*

-------------------------------------------------------

*See also: GHOAT.py, a fully automated tool for guest-host ABFE calculations using SDR with pmemd.cuda:* 

https://github.com/GHeinzelmann/GHOAT.py 

*A tutorial and a detailed user guide are available, as well as necessary parameter and input files for several hosts.*


# BAT.py v2.0

The Binding Affinity Tool (BAT.py) is a python tool for fully automated absolute binding free energy calculations. Its workflow encompasses the creation of the bound complex, generation of parameters using Antechamber, preparation of the simulation files, and post-processing to retrieve the binding free energy. By using the _pmemd.cuda_ software from AMBER20, it is able to perform several calculations at a reduced computational cost using graphics processing units (GPUs).


The 2.0 version of BAT.py can perform binding free energy calculations by two alchemical routes in the presence of restraints, either with double decoupling (DD) procedure or with the simultaneous decoupling recoupling (SDR) method. For the use of the APR method in addition to DD and SDR, download the 1.0 version of the code at the BATv1.0 branch. BAT.py is compatible with the simulation package AMBER20, also requiring a few installed programs to work properly, which are listed in the next section. 

![](doc/figure.png)

# Getting started

To use BAT.py, download the files from this repository, which already contain an example for ligand binding to the second bromodomain of the BRD4 protein - BRD4(2). In order to perform all the steps from BAT.py, the following programs must be installed and in your path:

VMD (Visual Molecular Dynamics) [2] - https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

Openbabel 2.4.1 [3] - https://github.com/openbabel/openbabel/releases/tag/openbabel-2-4-1 &dagger;

MUSTANG v3.2.3 (MUltiple (protein) STructural AligNment alGorithm) [4] - http://lcb.infotech.monash.edu.au/mustang/

AmberTools20 or later [5] - http://ambermd.org/AmberTools.php &Dagger;

_pmemd.cuda_ software from AMBER20 [5] - http://ambermd.org/GetAmber.php

&dagger; Had problems with wrong protonation using Openbabel 3, so keeping the 2.4.1 version for now, might change in the future. 

&Dagger; If pdb4amber from Ambertools does not work, add the following line to your amber.sh file: export PYTHONPATH=$PYTHONPATH:$AMBERHOME/lib/python3.8/site-packages/pdb4amber-1.7.dev0-py3.8.egg/

The folder ./BAT/all-poses contains an example of system input files, with a docked receptor from the 5uez crystal structure (LMCSS-5uf0\_5uez\_docked.pdb), as well as 9 docked poses for the ligand with the 5uf0 crystal structure (pose0.pdb to pose8.pdb). The docking files were generated and converted to .pdb using Autodock Vina and AutodockTools, following a protocol adapted from the CELPP challenge tutorial (https://docs.google.com/document/d/1iJcPUktbdrRftAA8cuVa32Ri1TPr2XvZVqTccDja2OM). Inside the ./all-poses folder there is also the original crystal structure file for 5uf0. Below we show an example of using these files to calculate the standard binding free energies of the top 5 docked poses and the crystal structure, with all the necessary steps in the calculation. 

# Running a sample calculation

The simulations and analysis from this example will be performed inside the ./BAT folder. The simulations are divided in two steps, equilibration (folder ./equil) and free energy calculation (folder ./fe). The input file with all the needed BAT.py parameters for double decoupling is called input-dd.in, with the meaning of each explained in more detail in the user guide, located inside the ./doc folder. For our sample calculation, we will use the values already provided in the input files included in this distribution. Briefly, the poses\_list parameter sets up the calculation for the first 5 poses from Autodock Vina, all in the ./all-poses folder. The input files can be modified to perform the calculations in the 5uf0 crystal structure, by changing the calc\_type option to "crystal", the celpp\_receptor option to "5uf0", and the ligand\_name option to "89J", which is the ligand residue name in the 5uf0 pdb structure. 

![](doc/workflow.png)

## Equilibration

The equilibration step starts from the docked complex or the crystal structure, gradually releasing restraints applied on the ligand and then performing a final simulation with an unrestrained ligand. The necessary simulation parameters for the ligand are also generated in this stage, using the General Amber Force Field versions 1 and 2 (GAFF or GAFF2) [6], and the AM1-BCC charge model [7,8]. To run this step, inside the program main folder type:

python BAT.py -i input-dd.in -s equil

BAT.py is compatible with python 3.8 versions. If you have another version, or you find that this command gives an error, you can use the python version included in the Ambertools20 distribution:

$AMBERHOME/miniconda/bin/python BAT.py -i input-dd.in -s equil

This command will create an ./equil folder, with one folder inside for each of the docked poses (pose0, pose1, etc.). In order to run the simulations for each pose, you can use the run-local.bash script (to run them locally), or the PBS-run script, which is designed to run in a queue system such as TORQUE. Both of these files might have to be adjusted, depending on your computer or server configuration, which can be done in the templates located in the ./BAT/run\_files folder. The number of simulations and the applied restraints will depend on the _release_eq_ array defined in the input file. 


## Free energy calculation 

### Simulations

The free energy stage starts from the equilibrated system, rebuilding the latter, as well as redefining the anchor atoms and the restraints for use in the free energy calculation. In this example we will use the double decoupling method (DD) with retraints to obtain the binding free energies. For charged ligands, one should use the simultaneous decoupling recoupling (SDR) method instead, as explained in the user guide. Again in the program main folder, type:

python BAT.py -i input-dd.in -s fe

For each pose or crystal structure, a folder will be created inside ./fe, and inside there will be two folders, ./restraints and ./dd. The restraints folder contains all the simulations needed for the application/removal of restraints. The ./dd folder contains the coupling/decoupling of the ligand electrostatic/LJ interactions, both in the binding site and in bulk. A script called run-all-dd.bash, inside the ./run_files folder, can be used to run these simulatons quickly using the PBS scripts provided. A similar script can be written to do the same, using your particular running protocol. 

### Analysis

Once all of the simulations are concluded, it is time to process the output files and obtain the binding free energies. Here a few parameters concerning the analysis can be set in the input file, such as using TI or MBAR [9] for double decoupling, number of blocks for block data analysis, and the Gaussian weights if TI is used for double decoupling. Inside the main folder type:

python BAT.py -i input-dd.in -s analysis

You should see a ./Results directory inside each ./fe/pose folder, containing all the components and the final calculated binding free energy, located in the Results.dat file. This folder also contains the results for each of the chosen data blocks, which is useful to check for convergence and fluctuations, and is also used to calculate the uncertainties. This fully automated procedure can be readily applied for any other ligand that binds to the second BRD4 bromodomain, and with minimal adjustments it can be extended to several other proteins.

## Using the SDR method

The SDR method is suitable for ligands with net charge, since it keeps the two free energy topologies with the same charge during the transformations. To apply SDR, a few parameters have to be changed or added, which is shown in the input-sdr.in file included in this example. All that is needed is to perform the free energy and analysis steps again with this input file. The DD and SDR methods should produce consistent results if the reference state is the same from the equilibrium stage.  

# Extending it to other systems

## Additional ligands to BRD4(2)

The sample system shown here uses a particular ligand that binds to the second bromodomain of the BRD4 protein - BRD4(2). The system alignment, parameter generation and assignment of the ligand anchor atoms is done automatically, so these same calculations can be extended to any other ligand that binds to this receptor. The only thing needed is the files in the ./all-poses folder to be changed, including the docked receptor and poses pdb files, as well as the crystal structure if desired.     

## Additional receptors

To include a new receptor system, some additional input data is needed. They include a reference.pdb file to align the system using MUSTANG, three chosen protein anchors, and possibly a few variables for ligand anchor atom search. These can be found inside the ./systems-library folder for three other bromodomains (CREBBP, BRD4(1) and BAZ2B), the T4 Lysozyme, and the MCL-1 protein. Other systems will be added with time, as the program is further tested and validated.    

# More information

A paper explaining the whole BAT.py theoretical aspects and calculation procedure is available in Ref [1]. For more information you can contact the author directly:

Germano Heinzelmann <br/>
Departamento de Física, Universidade Federal de Santa Catarina <br/>
Florianópolis - SC  88040-970 Brasil <br/>
email: germanohei@gmail.com <br/>

# Acknowledgments

Germano Heinzelmann thanks FAPESC and CNPq for the research grants, and Prof. Michael Gilson for the support on developing the code.

# References


1. G. Heinzelmann and M. K. Gilson (2021). “Automation of absolute protein-ligand binding free energy calculations for docking refinement and compound evaluation”. Scientific Reports, 11, 1116.

2. W. Humphrey, A. Dalke and K. Schulten. (1996)  "VMD - Visual Molecular Dynamics", Journal of Molecular Graphics, 14, 33-38.

3. N. M. O'Boyle, M. Banck, C. A. James, C. Morley, T. Vandermeersch, and G. R. HutchisonEmail. (2011) "Open Babel: An open chemical toolbox." Journal of Cheminformatics, 3, 33.

4. A. S. Konagurthu, J. Whisstock, P. J. Stuckey, and A. M. Lesk. (2006) “MUSTANG: A multiple structural alignment algorithm”. Proteins, 64, 559-574.

5. D.A. Case, K. Belfon, I.Y. Ben-Shalom, S.R. Brozell, D.S. Cerutti, T.E. Cheatham, III, V.W.D. Cruzeiro, T.A. Darden, R.E. Duke, G. Giambasu, M.K. Gilson, H. Gohlke, A.W. Goetz, R. Harris, S. Izadi, S.A. Izmailov, K. Kasavajhala, A. Kovalenko, R. Krasny, T. Kurtzman, T.S. Lee, S. LeGrand, P. Li, C. Lin, J. Liu, T. Luchko, R. Luo, V. Man, K.M. Merz, Y. Miao, O. Mikhailovskii, G. Monard, H. Nguyen, A. Onufriev, F.Pan, S. Pantano, R. Qi, D.R. Roe, A. Roitberg, C. Sagui, S. Schott-Verdugo, J. Shen, C. Simmerling, N.R.Skrynnikov, J. Smith, J. Swails, R.C. Walker, J. Wang, L. Wilson, R.M. Wolf, X. Wu, Y. Xiong, Y. Xue, D.M. York and P.A. Kollman (2020), AMBER 2020, University of California, San Francisco.

6. J. Wang, R.M. Wolf, J.W. Caldwell, and P. A. Kollman, D. A. Case (2004) "Development and testing of a general AMBER force field". Journal of Computational Chemistry, 25, 1157-1174.

7. A. Jakalian, B. L. Bush, D. B. Jack, and C.I. Bayly (2000) "Fast, efficient generation of high‐quality atomic charges. AM1‐BCC model: I. Method". Journal of Computational Chemistry, 21, 132-146.

8. A. Jakalian, D. B. Jack, and C.I. Bayly (2002) "Fast, efficient generation of high‐quality atomic charges. AM1‐BCC model: II. Parameterization and validation". Journal of Computational Chemistry, 16, 1623-1641.

9. M. R. Shirts and J. Chodera (2008) “Statistically optimal analysis of samples from multiple equilibrium states.” Journal of  Chemical Physics, 129, 129105.





