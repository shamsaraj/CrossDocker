#!/usr/bin/env python
#
# Adapter around AutoDockTools' receptor preparation, for CrossDocker.
#
# Derived from AutoDockTools Utilities24/prepare_receptor4.py:
# $Header: /opt/cvs/python/packages/share1.5/AutoDockTools/Utilities24/prepare_receptor4.py,v 1.11 2007/11/28 22:40:22 rhuey Exp $
# Copyright (c) Michel F. Sanner and TSRI. AutoDockTools is distributed under
# the MGLTools Software License Agreement -- see http://mgltools.scripps.edu
#
# Modifications by Jamal Shamsara: the command-line front end (option parsing
# and the usage text) has been removed, and the remaining call exposed as
# PR() with the options CrossDocker uses fixed as defaults.
#
# AutoDockTools is NOT bundled with CrossDocker. The MolKit and AutoDockTools
# imports below resolve against an MGLTools installation that you provide --
# see README.md.

from MolKit import Read
import MolKit.molecule
import MolKit.protein
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation


def PR (RF, OF):
    # initialize required parameters
    #-s: receptor
    receptor_filename =  RF

    # optional parameters
    verbose = None
    #-A: repairs to make: add bonds and/or hydrogens or checkhydrogens
    repairs = ''
    #-C default: add gasteiger charges 
    charges_to_add = 'gasteiger'
    #-p preserve charges on specific atom types
    preserve_charge_types=None
    #-U: cleanup by merging nphs_lps, nphs, lps, waters, nonstdres
    cleanup  = "nphs_lps_waters"####"nphs_lps_waters_nonstdres"
    #-o outputfilename
    outputfilename = OF
    #-m mode 
    mode = 'automatic'
    #-e delete every nonstd residue from each chain
    delete_single_nonstd_residues = None

    if not receptor_filename:
        raise ValueError('receptor filename must be specified')

    #what about nucleic acids???

    mols = Read(receptor_filename)
    if verbose: print 'read ', receptor_filename
    mol = mols[0]
    preserved = {}
    if charges_to_add is not None and preserve_charge_types is not None:
        preserved_types = preserve_charge_types.split(',') 
        if verbose: print "preserved_types=", preserved_types
        for t in preserved_types:
            if verbose: print 'preserving charges on type->', t
            if not len(t): continue
            ats = mol.allAtoms.get(lambda x: x.autodock_element==t)
            if verbose: print "preserving charges on ", ats.name
            for a in ats:
                if a.chargeSet is not None:
                    preserved[a] = [a.chargeSet, a.charge]

    if len(mols)>1:
        if verbose: print "more than one molecule in file"
        #use the molecule with the most atoms
        ctr = 1
        for m in mols[1:]:
            ctr += 1
            if len(m.allAtoms)>len(mol.allAtoms):
                mol = m
                if verbose: print "mol set to ", ctr, "th molecule with", len(mol.allAtoms), "atoms"
    mol.buildBondsByDistance()

    if verbose:
        print "setting up RPO with mode=", mode,
        print "and outputfilename= ", outputfilename
        print "charges_to_add=", charges_to_add
        print "delete_single_nonstd_residues=", delete_single_nonstd_residues

    RPO = AD4ReceptorPreparation(mol, mode, repairs, charges_to_add, 
                        cleanup, outputfilename=outputfilename,
                        preserved=preserved, 
                        delete_single_nonstd_residues=delete_single_nonstd_residues)    

    if charges_to_add is not None:
        #restore any previous charges
        for atom, chargeList in preserved.items():
            atom._charges[chargeList[0]] = chargeList[1]
            atom.chargeSet = chargeList[0]


# To execute this command type:
# prepare_receptor4.py -r pdb_file -o outputfilename -A checkhydrogens 

