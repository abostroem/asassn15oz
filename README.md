# Signatures of Circumstellar Interaction in the Type IIL Supernova ASASSN-15oz
This repository contains all of the automated scripts I used to produce the paper<sup>1</sup>: 
[Signatures of Circumstellar Interaction in the Type IIL Supernova ASASSN-15oz](https://ui.adsabs.harvard.edu//#abs/2019arXiv190109962A/abstract). 
All spectroscopic data of ASASSN-15oz is available on WISeREP (https://wiserep.weizmann.ac.il).
All photometric data of ASASSN-1oz is available in the SNDAVIS database (http://dark.physics.ucdavis.edu/sndavis/transient).
Some script utilize automated access to the SNDAVIS database that is private. If you are using
this part of the scripts, you will need to create alternate access to this data.

The code in this repository depends on the following third-party packages:
* extinction: https://extinction.readthedocs.io/en/latest/ (pip installable)
* utilities_az: https://github.com/abostroem/utilities (run setup.py)
     - side note: I turned this into a package late in the game, so if there is an import for spectroscopy,
     visualization, or supernova that isn't working, try importing it from utilities\_az
* line\_fitting: https://github.com/abostroem/line_fitting (not yet a package :(  )


## Approximate format of repository:
My workflow is generally to start in notebooks and then move to scripts when my code base
if stable. In general, the scripts supercede the analysis in the notebooks and the
scripts in paper/figures supercede everything else.
* notebooks: lab notebooks detailing all of the things I tried that didn't work and some of the things that did work
* code: scripts, usually python, usually analysis that was used in the paper
* paper: the tex documents, tables, and scripts to make the figures that finished in the paper
* data: **this is incomplete** I've only uploaded text files into the data directory. 
For fits files see two links above.

Feel free to contact me with questions:
abostroem@gmail.com

<sup>1</sup> Apologies in advance that this isn't as clean as I'd like it to be. Its my first paper and
attempt at reproducibility and I'm still learning ... a lot.