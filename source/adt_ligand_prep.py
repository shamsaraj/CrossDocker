#!/usr/bin/env python
#
# Adapter around AutoDockTools' ligand preparation, for CrossDocker.
#
# Derived from AutoDockTools Utilities24/prepare_ligand4.py:
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_ligand4.py,v 1.5.4.1 2009/04/15 17:41:57 rhuey Exp $
# Copyright (c) Michel F. Sanner and TSRI. AutoDockTools is distributed under
# the MGLTools Software License Agreement -- see http://mgltools.scripps.edu
#
# Modifications by Jamal Shamsara: the command-line front end (option parsing
# and the usage text) has been removed, and the remaining call exposed as
# PL() with the options CrossDocker uses fixed as defaults.
#
# AutoDockTools is NOT bundled with CrossDocker. The MolKit and AutoDockTools
# imports below resolve against an MGLTools installation that you provide --
# see README.md.

from MolKit import Read

from AutoDockTools.MoleculePreparation import AD4LigandPreparation


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


    if not ligand_filename:
        raise ValueError('ligand filename must be specified')

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








