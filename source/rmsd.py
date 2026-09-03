#!/usr/bin/env python
#
# Reduced from the "rmsd" project by Jimmy Charnley Kromann and Lars
# Bratholm: https://github.com/charnley/rmsd
#
# Modifications by Jamal Shamsara for CrossDocker: unused routines removed
# (Kabsch rotation, centroid fitting, coordinate output and the command-line
# interface), and rmsd2() added as the entry point CrossDocker calls.
#
# ---------------------------------------------------------------------------
# Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ---------------------------------------------------------------------------
#
# The full licence text is also in licenses/rmsd-BSD-2-Clause.txt.

import numpy as np
import re


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def get_coordinates(filename, ignore_hydrogens=False):
    """
    Get coordinates from filename.
    Get coordinates from a filename.xyz and return a vectorset with all the
    coordinates.
    This function has been written to parse XYZ files, but can easily be
    written to parse others.
    """
    f = open(filename, 'r')
    V = []
    atoms = []
    n_atoms = 0
    lines_read = 0

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(f.next())
    except ValueError:
        exit("Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    f.next()

    # Use the number of atoms to not read beyond the end of a file
    for line in f:

        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        numbers = re.findall(r'[-]?\d+\.\d*', line)
        numbers = [float(number) for number in numbers]

        # ignore hydrogens
        if ignore_hydrogens and atom.lower() == "h":
            continue

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) == 3:
            V.append(np.array(numbers))
            atoms.append(atom)
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

        lines_read += 1

    f.close()
    V = np.array(V)
    return atoms, V


def rmsd2(mol1, mol2):
    """Heavy-atom RMSD between two poses of the same molecule, in XYZ format.

    This is a direct, in-place RMSD: the coordinates are compared as they
    stand, with no superposition. That is what cross-docking needs -- the
    question is how far the docked pose sits from the crystallographic pose
    in the receptor's own frame, not how well the two shapes can be made to
    overlap.

    Both files must list their atoms in the same order, and the comparison
    is not symmetry-aware (equivalent atoms in, say, a phenyl ring are not
    matched up), which is the usual limitation of plain docking RMSD.
    """
    # Hydrogen positions are not determined by the docking, so compare heavy
    # atoms only.
    ignore_hydrogens = True

    atomsP, P = get_coordinates(mol1, ignore_hydrogens=ignore_hydrogens)
    atomsQ, Q = get_coordinates(mol2, ignore_hydrogens=ignore_hydrogens)

    return rmsd(P, Q)
