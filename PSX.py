#!/usr/bin/env python

###########
# Imports #
###########
import sys
import argparse
import os
import subprocess

# Check for and import BioPython modules
try:
	from Bio.PDB.PDBParser import PDBParser
	from Bio.PDB import PDBIO

except:
	sys.stderr.write("Unable to import modules from BioPython. The module is necessary for the program to run. Visit https://biopython.org/ for details and installation instructions. Exiting.\n")
	sys.exit()

# Local imports
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "Code"))

from ProtectionFactor import *
from Deuteration import *
from CheckInput import *

########
# Main #
########
if __name__ == "__main__":
	sys.stdout.write("\n")
	sys.stdout.write("\n")

	# Print splash
	Description = "***************************************************************************\n" + \
                  "* PSX                                                                     *\n" + \
                  "*                                                                         *\n" + \
                  "* Software for estimating protium/deuterium exchange in proteins. Run the *\n" + \
                  "* software with the argument -h to see a list of options and arguments.   *\n" + \
                  "*                                                                         *\n" + \
                  "* If you use PSX in your research, please reference:                      *\n" + \
                  "* Pedersen, Wang, Tidemand, Martel, Lindorff-Larsen, & Arleth (2019)      *\n" + \
<<<<<<< HEAD
                  "* J. Appl. Crystallogr. 52, 1427-1436                                     *\n" + \
=======
                  "* J. Appl. Crystallgr. 52, 1427–1436                                      *\n" + \
>>>>>>> 9010dd5ac07c450f90f0c85dd21ce73903456ccb
                  "***************************************************************************\n"

	sys.stdout.write(Description)

	sys.stdout.write("\n")
	PrintBreak()

	# Assign CLI arguments
	Parser = argparse.ArgumentParser()

	# Path to the provided pdb-file
	Parser.add_argument("p", metavar = "Path to .pdb-file", help = "Path to atomic coordinates in .pdb-format")

	# Parameters used for calculating protection factors
	Parser.add_argument("-dc", metavar = "", help = "Cut-off for degree of burial calculation in angstrom, default is 6.5"                                        , default = "6.5"    , type = float)
	Parser.add_argument("-eh", metavar = "", help = "Cut-off for hydrogen bond energy in DSSP algorithm in kcal/mol, default is -0.5"                             , default = "-0.5"   , type = float)
	Parser.add_argument("-bc", metavar = "", help = "Coefficient for degree of burial in protection factor calculation, beta_c, default is 0.35"                  , default = "0.35"   , type = float)
	Parser.add_argument("-bh", metavar = "", help = "Coefficient for hydrogen bonds in protection factor calculation, beta_h, default is 2.0"                     , default = "2.0"    , type = float)

	# Parameters used for estimating intrinsic exchange rates
	Parser.add_argument("-te", metavar = "", help = "Temperature used to estimate exchange rates in degrees Kelvin, default is 293.15"                            , default = "293.15" , type = float)
	Parser.add_argument("-ph", metavar = "", help = "pH used to estimate exchange rates as measured by ph-meter, default is 7.0"                                  , default = "7.0"    , type = float)
	Parser.add_argument("-rd", metavar = "", help = "Reference data used in estimation of intrinsic exchange rates; options are poly and oligo, default is poly"  , default = "poly"   , type = CheckRDFormat)

	# Parameters used for H/D exchange
	Parser.add_argument("-pc", metavar = "", help = "Probability cut-off above which backbone hydrogens will be exchanged, default is 0.5"                        , default = "0.5"    , type = float)
	Parser.add_argument("-ti", metavar = "", help = "Time elapsed since protein was exposed to the solvent in seconds, default is 0.0"                            , default = "0.0"    , type = float)

	# Parameter used to locate DSSP and pdb2pqr
	Parser.add_argument("-dp", metavar = "", help = "Path to DSSP executable, default is dssp"                                                                    , default = "dssp"   , type = str)
	Parser.add_argument("-pp", metavar = "", help = "Path to pdb2pqr executable, default is pdb2pqr"                                                              , default = "pdb2pqr", type = str)

	# Parameter for outputting warnings in the import of the pdb
	Parser.add_argument("-pw", metavar = "", help = "If set to 1, the software will output warnings from the import of the given atomic coordinates, default is 0", default = "0"      , type = int)

	# Extract parameters
	Arguments = Parser.parse_args()

	PathToPDBFile     = Arguments.p
	DistanceCutoff    = Arguments.dc
	EnergyCutoff      = Arguments.eh
	Betac             = Arguments.bc
	Betah             = Arguments.bh
	Temperature       = Arguments.te
	pH                = Arguments.ph
	ReferenceData     = Arguments.rd
	DeuterationCutoff = Arguments.pc
	TimeElapsed       = Arguments.ti
	PathToDSSP        = Arguments.dp
	PathTopdb2pqr     = Arguments.pp
	OutputWarnings    = not bool(Arguments.pw)

	# Check for DSSP availability
	if not subprocess.call(PathToDSSP + " --version", shell = True, stdout = open(os.devnull, 'wb')):
		sys.stdout.write("DSSP located. \n")

	else:
		sys.stderr.write("Unable to locate and run DSSP; the provided path is {}. The software is necessary for the program to run. Visit https://swift.cmbi.umcn.nl/gv/dssp/ for details and installation instructions. Exiting.\n".format(PathToDSSP))
		sys.exit()

	# Check for pdb2pqr availability
	if not subprocess.call(PathTopdb2pqr + " --version", shell = True, stdout = open(os.devnull, 'wb')):
		sys.stdout.write("pdb2pqr located. \n")

	else:
		sys.stderr.write("Unable to locate and run pdb2pqr; the provided path is {}. The software is necessary for the program to run. Visit https://www.ics.uci.edu/~dock/pdb2pqr/userguide.html for details and installation instructions. Exiting.\n".format(PathTopdb2pqr))
		sys.exit()

	PrintBreak()

	# Extract filename
	PathToOutputFiles = os.path.splitext(PathToPDBFile)[0]

	# Open and parse the protein file (permissive read-in)
	InputParser      = PDBParser(PERMISSIVE = True, QUIET = OutputWarnings)
	ProteinStructure = InputParser.get_structure(PathToPDBFile, PathToPDBFile)

	# Check input
	CheckInput(ProteinStructure)

	# Run protection factor calculation
	IntrinsicExhangeRatesAndProtectionFactors(ProteinStructure, DistanceCutoff, Temperature, pH, ReferenceData, EnergyCutoff, PathToDSSP, Betac, Betah)

	# Output intermediate structure
	OutputParser = PDBIO()
	OutputParser.set_structure(ProteinStructure)
	OutputParser.save(PathToOutputFiles + ".prf.pdb")
	sys.stdout.write("Exported protection factors and intrinsic exchange rates to:\n{}.prf.pdb.\n".format(PathToOutputFiles))
	PrintBreak()

	# Average protection factors and assign the to each model
	AverageProtectionFactors(ProteinStructure)

	# Exchange
	ExchangeModels(ProteinStructure, PathToOutputFiles, PathTopdb2pqr, DeuterationCutoff, TimeElapsed, pH)

	# Output deuterated structure
	OutputParser = PDBIO()
	OutputParser.set_structure(ProteinStructure)
	OutputParser.save(PathToOutputFiles + ".deu.pdb")
	sys.stdout.write("Exported deuterated structure and extrinsic exchange rates to:\n{}.deu.pdb.\n".format(PathToOutputFiles))
	PrintBreak()

	# Clear
	sys.stdout.write("\n")
	sys.stdout.write("\n")
