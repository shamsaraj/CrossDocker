# CrossDocker

Citation:

https://springerplus.springeropen.com/articles/10.1186/s40064-016-1972-4

Shamsara, J. (2016). "CrossDocker: a tool for performing cross-docking using Autodock Vina." SpringerPlus 5(1): 1.

```bibtex
@article{shamsara2016crossdocker,
  author  = {Shamsara, Jamal},
  title   = {{CrossDocker}: a tool for performing cross-docking using {Autodock Vina}},
  journal = {SpringerPlus},
  year    = {2016},
  volume  = {5},
  pages   = {1},
  doi     = {10.1186/s40064-016-1972-4}
}
```

Download:

The executable, bundled Vina binary, sample datasets, and standalone source
bundle are attached to the [v1.0 release](https://github.com/shamsaraj/CrossDocker/releases/tag/v1.0)
rather than committed in this repository (they're binaries/archives, so
Releases is a better fit than the git tree). The `source/` folder here has
the plain Python source for reference/review.

Installing:

The executable was built for windows. There is no need to installation; just extract contents of the CrossDocker.zip archieve in desired directory.

Prerequisites:

**Open Babel** 2.3 or higher is needed for successful execution of CrossDocker.
By default `babel.exe` is expected on the `PATH` environment variable; if it is
not, add it before running CrossDocker.

**MGLTools / AutoDockTools** must also be installed. CrossDocker does not bundle
it: `source/adt_ligand_prep.py` and `source/adt_receptor_prep.py` are thin
adapters that import `MolKit` and `AutoDockTools` from your own MGLTools
installation. Download it from <http://mgltools.scripps.edu> and make sure its
Python packages are importable by the interpreter running CrossDocker.

Usage:

Set the parameters defined in the config.txt file.
Run the CrossDocker.exe file from command line.
To test the CrossDocker the small_dataset can be used.

Tip: `set_vina_priority_high.bat` / `set_vina_priority_low.bat` set the
running vina.exe process to high or idle CPU priority (via `wmic ...
setpriority`), if you want docking to run faster or stay in the background
without slowing down other programs.

License:

CrossDocker's own code is MIT licensed (see LICENSE).

It also uses third-party work: an RMSD routine reduced from the BSD-2-Clause
[rmsd](https://github.com/charnley/rmsd) project, adapters derived from
AutoDockTools (Copyright Michel F. Sanner / TSRI, MGLTools Software License
Agreement), and a bundled copy of AutoDock Vina in the release archive under
its own separate licence (`vina_license.rtf`).

Each component, its licence and the modifications made are listed in
[NOTICE.md](NOTICE.md).

