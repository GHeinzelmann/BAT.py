#!/usr/bin/env python

__author__ = 'germano.heinzelmann@ufsc.br'


from d3r.celppade.custom_dock import Dock

class autodockvina(Dock):
    """Abstract class defining methods for a custom docking solution
    for CELPP
    """
    Dock.SCI_PREPPED_LIG_SUFFIX = '_prepared.pdbqt'
    Dock.SCI_PREPPED_PROT_SUFFIX = '_prepared.pdbqt'


    def ligand_technical_prep(self, sci_prepped_lig, targ_info_dict = {}):
        """
        'Technical preparation' is the step immediate preceding
        docking. During this step, you may perform any file
        conversions or processing that are specific to your docking
        program. Implementation of this function is optional.
        :param sci_prepped_lig: Scientifically prepared ligand file
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: A list of result files to be copied into the
        subsequent docking folder. The base implementation merely
        returns the input string in a list (ie. [sci_prepped_lig]) 
        """
        return super(autodockvina,
                     self).ligand_technical_prep(sci_prepped_lig,
                                                 targ_info_dict = targ_info_dict)

    def receptor_technical_prep(self, 
                                sci_prepped_receptor, 
                                pocket_center, 
                                targ_info_dict = {}):
        """
        'Technical preparation' is the step immediately preceding
        docking. During this step, you may perform any file
        conversions or processing that are specific to your docking
        program. Implementation of this function is optional.
        :param sci_prepped_receptor: Scientifically prepared receptor file
        :param pocket_center: list of floats [x,y,z] of predicted pocket center
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: A list of result files to be copied into the
        subsequent docking folder. This implementation merely
        returns the input string in a list (ie [sci_prepped_receptor])
        """
        
        #return [sci_prepped_receptor, pocket_center]
        return super(autodockvina,
                     self).receptor_technical_prep(sci_prepped_receptor, 
                                                   pocket_center,
                                                   targ_info_dict=targ_info_dict)




    def dock(self, 
             tech_prepped_lig_list, 
             tech_prepped_receptor_list, 
             output_receptor_pdb, 
             output_lig_mol, 
             targ_info_dict={}):
        """
        This function is the only one which the contestant MUST
        implement.  The dock() step runs the actual docking
        algorithm. Its first two arguments are the return values from
        the technical preparation functions for the ligand and
        receptor. These arguments are lists of file names (strings),
        which can be assumed to be in the current directory. 
        If prepare_ligand() and ligand_technical_prep() are not
        implemented by the contestant, tech_prepped_lig_list will
        contain a single string which names a SMILES file in the
        current directory.
        If receptor_scientific_prep() and receptor_technical_prep() are not
        implemented by the contestant, tech_prepped_receptor_list will
        contain a single string which names a PDB file in the current
        directory.
        The outputs from this step must be two files - a pdb with the
        filename specified in the output_receptor_pdb argument, and a
        mol with the filename specified in the output_ligand_mol
        argument.
        :param tech_prepped_lig_list: The list of file names resturned by ligand_technical_prep. These have been copied into the current directory.
        :param tech_prepped_receptor_list: The list of file names resturned by receptor_technical_prep. These have been copied into the current directory.
        :param output_receptor_pdb: The final receptor (after docking) must be converted to pdb format and have exactly this file name.
        :param output_lig mol: The final ligand (after docking) must be converted to mol format and have exactly this file name.
        :param targ_info_dict: A dictionary of information about this target and the candidates chosen for docking.
        :returns: True if docking is successful, False otherwise. Unless overwritten, this implementation always returns False
        """

        print targ_info_dict
        receptor_pdbqt = tech_prepped_receptor_list[0]
        ligand_pdbqt = tech_prepped_lig_list[0]
        pocket_center = targ_info_dict['pocket_center']
        
        vina_command = ('vina --receptor ' + receptor_pdbqt + '  --ligand  ' +
                         ligand_pdbqt + ' --center_x ' + str(pocket_center[0]) +
                        ' --center_y ' + str(pocket_center[1]) + 
                        ' --center_z ' + str(pocket_center[2]) + 
                        ' --size_x 10 --size_y 10 --size_z 10 --seed 999 ' +
                        ' 1> vina.stdout 2> vina.stderr')
        print "Running: " + vina_command
        os.system(vina_command)

        out_dock_file = ligand_pdbqt.replace('.pdbqt','_out.pdbqt')
        
        

        poses = open(out_dock_file).read()
        x = poses.split('ENDMDL')
        for i in range(0, len(x)-1):
           x[i] = os.linesep.join([s for s in x[i].splitlines() if s])
           po_file = open('pose%s.pdbqt' %i, 'w')
           po_file.write (x[i])
           ovar = 'echo "\nENDMDL" >> pose%s.pdbqt' %i
           ovar2 = '. $MGL_ROOT/bin/mglenv.sh; python $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f pose%s.pdbqt -o pose%s.pdb' %(i, i)     
           po_file.close()
           os.system(ovar)
           os.system(ovar2)

        os.system("sed -e '/ENDMDL/,$d' " + out_dock_file + " > top_pose.pdbqt")
        os.system("echo ENDMDL >> top_pose.pdbqt")
        
        os.system('. $MGL_ROOT/bin/mglenv.sh; python $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f top_pose.pdbqt -o top_pose.pdb') 
        os.system("babel -ipdb top_pose.pdb -omol " + output_lig_mol)

        os.system('. $MGL_ROOT/bin/mglenv.sh; python $MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24/pdbqt_to_pdb.py -f ' + receptor_pdbqt + ' -o ' + output_receptor_pdb) 

        





if ("__main__") == (__name__):
    import os
    import logging
    import shutil
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-l", "--ligsciprepdir", metavar="PATH", help = "PATH where we can find the scientific ligand prep output")
    parser.add_argument("-p", "--protsciprepdir", metavar="PATH", help = "PATH where we can find the scientific protein prep output")
    parser.add_argument("-o", "--outdir", metavar = "PATH", help = "PATH where we will put the docking output")
    # Leave option for custom logging config here
    logger = logging.getLogger()
    logging.basicConfig( format  = '%(asctime)s: %(message)s', datefmt = '%m/%d/%y %I:%M:%S', filename = 'final.log', filemode = 'w', level   = logging.INFO )
    opt = parser.parse_args()
    lig_sci_prep_dir = opt.ligsciprepdir
    prot_sci_prep_dir = opt.protsciprepdir
    dock_dir = opt.outdir
    #running under this dir
    abs_running_dir = os.getcwd()
    log_file_path = os.path.join(abs_running_dir, 'final.log')
    log_file_dest = os.path.join(os.path.abspath(dock_dir), 'final.log')
    docker = autodockvina()
    docker.run_dock(prot_sci_prep_dir,
                    lig_sci_prep_dir,
                    dock_dir)
    #move the final log file to the result dir
    shutil.move(log_file_path, log_file_dest)
