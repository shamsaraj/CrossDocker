#!/usr/bin/env python
#
# 
#
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_ligand4.py,v 1.5.4.1 2009/04/15 17:41:57 rhuey Exp $
#
import os 

from MolKit import Read

from AutoDockTools.MoleculePreparation import AD4LigandPreparation




import sys


def usage():
    "Print helpful, accurate usage statement to stdout."
    print "Usage: prepare_ligand4.py -l filename"
    print
    print "    Description of command..."
    print "         -l     ligand_filename (.pdb or .mol2 or .pdbq format)"
    print "    Optional parameters:"
    print "        [-v]    verbose output"
    print "        [-o pdbqt_filename] (default output filename is ligand_filename_stem + .pdbqt)"
    print "        [-d]    dictionary to write types list and number of active torsions "

    print "        [-A]    type(s) of repairs to make:\n\t\t bonds_hydrogens, bonds, hydrogens (default is to do no repairs)"
    print "        [-C]    do not add charges (default is to add gasteiger charges)"
    print "        [-p]    preserve input charges on atom type, eg -p Zn"
    print "               (default is not to preserve charges on any specific atom type)"
    print "        [-U]    cleanup type:\n\t\t nphs_lps, nphs, lps, '' (default is 'nphs_lps') "
    print "        [-B]    type(s) of bonds to allow to rotate "
    print "               (default sets 'backbone' rotatable and 'amide' + 'guanidinium' non-rotatable)"
    print "        [-R]    index for root"
    print "        [-F]    check for and use largest non-bonded fragment (default is not to do this)"
    print "        [-M]    interactive (default is automatic output)"
    print "        [-I]    string of bonds to inactivate composed of "
    print "                   of zero-based atom indices eg 5_13_2_10  "
    print "                   will inactivate atoms[5]-atoms[13] bond "
    print "                               and atoms[2]-atoms[10] bond "
    print "                      (default is not to inactivate any specific bonds)"
    print "        [-Z]    inactivate all active torsions     "
    print "                      (default is leave all rotatable active except amide and guanidinium)"
    print "        [-g]    attach all nonbonded fragments "
    print "                      (default is not to do this)"


# process command arguments
#try:
    #opt_list, args = getopt.getopt(sys.argv[1:], 'l:vo:d:A:Cp:U:B:R:MFI:Zgh')
#except getopt.GetoptError, msg:
    #print 'prepare_ligand4.py: %s' %msg
    #usage()
    #sys.exit(2)
def PL (LF):
    # initialize required parameters
    #-l: ligand
    ligand_filename =  LF
    # optional parameters
    verbose = None
    add_bonds = False
    #-A: repairs to make: add bonds and/or hydrogens
    repairs = ""
    #-C  default: add gasteiger charges
    charges_to_add = 'gasteiger'
    #-p preserve charges on specific atom types
    preserve_charge_types=''
    #-U: cleanup by merging nphs_lps, nphs, lps
    cleanup  = "nphs_lps"
    #-B named rotatable bond type(s) to allow to rotate
    #allowed_bonds = ""
    allowed_bonds = "backbone"
    #-r  root
    root = 'auto'
    #-o outputfilename
    outputfilename = None
    #-F check_for_fragments
    check_for_fragments = False
    #-I bonds_to_inactivate
    bonds_to_inactivate = ""
    #-Z inactivate_all_torsions
    inactivate_all_torsions = False
    #-g attach_nonbonded_fragments
    attach_nonbonded_fragments = False
    #-m mode
    mode = 'automatic'
    #-d dictionary
    dict = None


    if not  ligand_filename:
        print 'prepare_ligand4: ligand filename must be specified.'
        usage()
        sys.exit()

    mols = Read(ligand_filename)
    if verbose: print 'read ', ligand_filename
    mol = mols[0]
    if len(mols)>1:
        if verbose:
            print "more than one molecule in file"
        #use the one molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose:
                    print "mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms"
    coord_dict = {}
    for a in mol.allAtoms: coord_dict[a] = a.coords


    mol.buildBondsByDistance()
    if charges_to_add is not None:
        preserved = {}
        preserved_types = preserve_charge_types.split(',')
        for t in preserved_types:
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]



    if verbose:
        print "setting up LPO with mode=", mode,
        print "and outputfilename= ", outputfilename
        print "and check_for_fragments=", check_for_fragments
        print "and bonds_to_inactivate=", bonds_to_inactivate
    LPO = AD4LigandPreparation(mol, mode, repairs, charges_to_add,
                            cleanup, allowed_bonds, root,
                            outputfilename=outputfilename,
                            dict=dict, check_for_fragments=check_for_fragments,
                            bonds_to_inactivate=bonds_to_inactivate,
                            inactivate_all_torsions=inactivate_all_torsions,
                            attach_nonbonded_fragments=attach_nonbonded_fragments)
    #do something about atoms with too many bonds (?)
    #FIX THIS: could be peptide ligand (???)
    #          ??use isPeptide to decide chargeSet??
    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]
    #if verbose: print "returning ", mol.returnCode
    #bad_list = []
    #for a in mol.allAtoms:
    #    if a.coords!=coord_dict[a]: bad_list.append(a)
    #if len(bad_list):
    #    print len(bad_list), ' atom coordinates changed!'
    #    for a in bad_list:
    #        print a.name, ":", coord_dict[a], ' -> ', a.coords
    #else:
    #    if verbose: print "No change in atomic coordinates"
    #if mol.returnCode!=0:
    #    sys.stderr.write(mol.returnMsg+"\n")
    #sys.exit(mol.returnCode)








