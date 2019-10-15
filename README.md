# PSX            

Copyright 2019, University of Copenhagen

Martin Cramer Pedersen, mcpe@nbi.ku.dk

This file is part of PSX.

PSX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PSX is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PSX. If not, see <http://www.gnu.org/licenses/>.

If you use PSX in your work, please cite:

Pedersen, Wang, Tidemand, Martel, Lindorff-Larsen, & Arleth (2019)
Some journal XX(YY), ZZZ-WWW

## Table of Contents

 - Dependencies
 - Command line arguments
 - Example

## Dependencies

In order for the PSX to function correctly, the following must be available:
 - Python 2.7, http://www.python.org/download/ 
 - BioPython, https://biopython.org/
 - DSSP, https://swift.cmbi.umcn.nl/gv/dssp/
 - pdb2pqr, https://www.ics.uci.edu/~dock/pdb2pqr/userguide.html

These are readily installed from Ubuntu's repositories by running:

    sudo apt install python python-biopython dssp pdb2pqr

## Command line arguments

The full list of optional command line arguments is as follows:

| Argument | Description                                                                                        | Default |
| -------- | -------------------------------------------------------------------------------------------------- | ------- |
| -bc      | Coefficient for degree of burial in protection factor calculation, beta_c                          | 0.35    |
| -bh      | Coefficient for hydrogen bonds in protection factor calculation, beta_h                            | 2.0     |
| -dc      | CCut-off for degree of burial calculation in angstrom                                              | 6.5     |
| -dp      | Path to DSSP executable                                                                            | dssp    |
| -eh      | Cut-off for hydrogen bond energy in DSSP algorithm in kcal/mol                                     | -0.5    |
| -pc      | Probability cut-off above which backbone hydrogens will be exchanged                               | 0.5     |
| -pf      | If set to 1, protection factor calculation will be skipped and occupancies will be used as log(Pf) | 0       |
| -ph      | pH used to estimate exchange rates as measured by ph-meter                                         | 7.0     |
| -pp      | Path to pdb2pqr executable                                                                         | pdb2pqr |
| -pw      | If set to 1, the software will output warnings from the import of the given atomic coordinates     | 0       |
| -rd      | Reference data used in estimation of intrinsic exchange rates; options are poly and oligo          | poly    |
| -te      | Temperature used to estimate exchange rates in degrees Kelvin                                      | 293.15  |
| -ti      | Time elapsed since protein was exposed to the solvent in seconds                                   | 0.0     |

## Example

PSX is run via the command line. In order to run the software on the 2lyz structure of lyzosome under the default conditions, one simply needs to run e.g. (assuming the file is present in the working directory):

    python PSX.py 2lyz.pdb

in a suitable terminal. Should one want to input e.g. a pH of 8.0 for the calculation, the corresponding command would be:

    python PSX.py 2lyz.pdb -ph 8.0

The outputs of the software is described in the paper referenced above.
