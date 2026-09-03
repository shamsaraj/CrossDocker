#!/usr/bin/env python
#
# center_of_mass.py -- mass-weighted centre of a PDB structure.
#
# Copyright (c) 2016-2026 Jamal Shamsara. MIT licensed (see LICENSE).
#
# Used by CrossDocker to place the centre of the Vina search box on the
# native ligand of each receptor.

# Standard atomic weights (IUPAC), restricted to the elements that occur in
# drug-like ligands and in the metal sites CrossDocker patches charges for.
# Add an entry here if a structure uses an element not listed.
ATOMIC_WEIGHTS = {
    "H": 1.008,     "B": 10.81,     "C": 12.011,    "N": 14.007,
    "O": 15.999,    "F": 18.998403, "NA": 22.989769, "MG": 24.305,
    "SI": 28.085,   "P": 30.973762, "S": 32.06,     "CL": 35.45,
    "K": 39.0983,   "CA": 40.078,   "MN": 54.938043, "FE": 55.845,
    "CU": 63.546,   "ZN": 65.38,    "SE": 78.971,   "BR": 79.904,
    "I": 126.90447,
}


def _element_of(line):
    """Return the element symbol for one PDB ATOM/HETATM line, uppercased.

    Prefers the element field (columns 77-78). Older files and some
    converters leave it blank, so fall back to the atom name (columns
    13-16), trying a two-letter symbol before a one-letter one so CL and
    BR are not read as carbon and boron.
    """
    element = line[76:78].strip().upper()
    if element in ATOMIC_WEIGHTS:
        return element

    name = line[12:16].strip().upper()
    alpha = "".join(ch for ch in name if ch.isalpha())
    if len(alpha) >= 2 and alpha[:2] in ATOMIC_WEIGHTS:
        return alpha[:2]
    if alpha[:1] in ATOMIC_WEIGHTS:
        return alpha[:1]

    raise ValueError(
        "Unknown element in PDB line (add it to ATOMIC_WEIGHTS):\n" + line.rstrip()
    )


def center_of_mass(pdbfile, include="ATOM,HETATM"):
    """Mass-weighted centre of the atoms in ``pdbfile``.

    ``include`` is a comma-separated list of PDB record names to read.
    Returns ``[x, y, z]`` in Angstroms, rounded to 3 decimals.
    """
    records = tuple(include.split(","))
    total_mass = 0.0
    weighted = [0.0, 0.0, 0.0]

    with open(pdbfile, "r") as pdb:
        for line in pdb:
            if not line.startswith(records):
                continue
            mass = ATOMIC_WEIGHTS[_element_of(line)]
            weighted[0] += mass * float(line[30:38])
            weighted[1] += mass * float(line[38:46])
            weighted[2] += mass * float(line[46:54])
            total_mass += mass

    if total_mass == 0.0:
        raise ValueError("No %s records found in %s" % (include, pdbfile))

    return [round(axis / total_mass, 3) for axis in weighted]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Mass-weighted centre of a PDB structure."
    )
    parser.add_argument("PDBfile")
    parser.add_argument(
        "-i", "--include", default="ATOM,HETATM", metavar="record-name",
        help='PDB records to include (default: "ATOM,HETATM")',
    )
    args = parser.parse_args()
    print(center_of_mass(args.PDBfile, include=args.include))
