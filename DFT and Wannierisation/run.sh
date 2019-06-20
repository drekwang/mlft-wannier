#!/bin/bash
MATERIAL="NiO"
/home/dirac/tmcs/newc5616/qe_release_6.4/PW/src/pw.x < $MATERIAL.scf > scf.out &&
/home/dirac/tmcs/newc5616/qe_release_6.4/PW/src/pw.x < $MATERIAL.nscf > nscf.out &&
/home/dirac/tmcs/newc5616/wannier90-3.0.0/wannier90.x -pp $MATERIAL &&
/home/dirac/tmcs/newc5616/qe_release_6.4/PP/src/pw2wannier90.x < $MATERIAL.pw2wan > pw2wan.out &&
/home/dirac/tmcs/newc5616/wannier90-3.0.0/wannier90.x $MATERIAL

