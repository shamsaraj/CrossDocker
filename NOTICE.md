# Third-party components

CrossDocker's own code is MIT licensed (see [LICENSE](LICENSE)). It also uses,
and in some cases includes, work by others. Each is listed here with its
origin, its licence, and what was changed.

## Included in this repository

### `source/rmsd.py`
RMSD calculation, reduced from the **rmsd** project by Jimmy Charnley Kromann
and Lars Bratholm — <https://github.com/charnley/rmsd>.

- Licence: **BSD 2-Clause**, Copyright (c) 2013 Jimmy Charnley Kromann
  &lt;jimmy@charnley.dk&gt; &amp; Lars Bratholm. Full text:
  [licenses/rmsd-BSD-2-Clause.txt](licenses/rmsd-BSD-2-Clause.txt); the notice
  and disclaimer are also reproduced at the top of the file itself.
- Modifications: routines CrossDocker does not use were removed (Kabsch
  rotation, centroid fitting, coordinate output, and the command-line
  interface). `rmsd2()` was added as the entry point CrossDocker calls. The
  retained functions `rmsd()` and `get_coordinates()` are unmodified.

### `source/adt_ligand_prep.py`, `source/adt_receptor_prep.py`
Adapters around AutoDockTools' ligand and receptor preparation, derived from
`Utilities24/prepare_ligand4.py` and `Utilities24/prepare_receptor4.py`.

- Licence: Copyright (c) Michel F. Sanner and TSRI; AutoDockTools is
  distributed under the MGLTools Software License Agreement —
  <http://mgltools.scripps.edu>.
- Modifications: the command-line front end (option parsing and usage text)
  was removed and the remaining call exposed as `PL()` / `PR()` with the
  options CrossDocker uses fixed as defaults. The original CVS `$Header:`
  line is retained in each file as the provenance record.
- **AutoDockTools itself is not bundled.** These files import `MolKit` and
  `AutoDockTools` from an MGLTools installation that you provide — see the
  Prerequisites section of [README.md](README.md).

## Bundled in the v1.0 release archive

### `vina.exe`
AutoDock Vina, used as an external docking engine and bundled as a
convenience in the [v1.0 release](https://github.com/shamsaraj/CrossDocker/releases/tag/v1.0).
Distributed under its own separate licence — see `vina_license.rtf`.

## Note on the v1.0 release binary

`CrossDocker.exe` in the v1.0 release was built with py2exe in 2016 and embeds
a copy of the AutoDockTools modules it imported at build time. That binary
cannot practically be rebuilt on a current system, so it is kept as a legacy
artefact. The source in this repository reflects the current arrangement, in
which AutoDockTools is a user-supplied prerequisite rather than something
CrossDocker redistributes.

## Previously included

An earlier revision of this repository contained `source/center.py`, a
centre-of-mass helper taken from Sebastian Raschka's
[protein-science](https://github.com/rasbt/protein-science) repository, which
is GPL-3.0 licensed and therefore incompatible with CrossDocker's MIT licence.
It has been removed and replaced by `source/center_of_mass.py`, an independent
implementation. The two agree to the 3 decimal places CrossDocker uses.
