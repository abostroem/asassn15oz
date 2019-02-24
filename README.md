This repository contains all of the automated scripts I used to produce the paper: 
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
* line\_fitting: https://github.com/abostroem/line\_fitting (not yet a package :(  )

Feel free to contact me with questions:
abostroem@gmail.com