###########
# Imports #
###########
import sys
import math

from CheckInput import *
from RemoveAndAddHydrogens import *

#######################
# Exchange dictionary #
#######################
def GenerateExchangeDictionary():
	ExchangeDictionary = {}

	# Backbone hydrogens - this should not be a tuple
	ExchangeDictionary["Backbone"] = "H"

	# N-terminal hydrogens
	ExchangeDictionary["NTerminal"] = ("H", "H2", "H3")

	# N-H in ASN, GLN, ARG, LYS, HIS, TRP
	# ARG
	ExchangeDictionary["ARG"] = ("HE", "HH11", "HH12", "HH21", "HH22")

	# ASN
	ExchangeDictionary["ASN"] = ("HD21", "HD22")

	# GLN
	ExchangeDictionary["GLN"] = ("HE21", "HE22")

	# HIS
	ExchangeDictionary["HIS"] = ("HD1", "HD2", "HE1", "HE2")

	# LYS
	ExchangeDictionary["LYS"] = ("HZ1", "HZ2", "HZ3")

	# TRP - only two of these should be occupied
	ExchangeDictionary["TRP"] = ("HE1", "HE2", "HE3")

	# O-H in SER, THR, TYR
	# SER
	ExchangeDictionary["SER"] = ("HG", )

	# THR
	ExchangeDictionary["THR"] = ("HG1", )

	# TYR
	ExchangeDictionary["TYR"] = ("HH", )

	# S-H in CYS
	# CYS
	ExchangeDictionary["CYS"] = ("HG", )

	return ExchangeDictionary

############
# Exchange #
############
def ExchangeSite(Residue, ID, Element = "D"):
	# Attempt to exchange the given hydrogen
	try:
		Residue[ID].element = Element
		return 1

	except:
		return 0

def DoExchange(ProteinModel, DeuterationCutoff, TimeElapsed):
	ExchangedBackboneHydrogens    = 0
	ExchangedNonBackboneHydrogens = 0
	ExchangedNTerminalHydrogens   = 0

	# Set up exchange dictionary
	ExchangeDictionary = GenerateExchangeDictionary()

	# Loop over all residues in the structure
	for Chain in ProteinModel:

		# Exchange N-terminal
		for ExchangeableSite in ExchangeDictionary["NTerminal"]:
			ExchangedNTerminalHydrogens += ExchangeSite(list(Chain.get_residues())[0], ExchangeableSite)
		
		# Exchange the rest of the chain
		for Residue in Chain:
			AminoAcidCode = Residue.get_resname()

			# Exchange rapidly-exchangeable hydrogens
			if AminoAcidCode in ExchangeDictionary.keys():

				for ExchangeableSite in ExchangeDictionary[AminoAcidCode]:
					ExchangedNonBackboneHydrogens += ExchangeSite(Residue, ExchangeableSite)

			# Deuterate backbone if it is above cutoff
			if ExchangeDictionary["Backbone"] in [Atom.name for Atom in list(Residue.get_atoms())]:
				LogProtectionFactor   = Residue["CA"].get_occupancy()
				IntrisicExchangeRate  = Residue["CA"].get_bfactor()
				ExtrinsicExchangeRate = IntrisicExchangeRate / math.exp(LogProtectionFactor)
				ProbabilityOfExchange = 1.0 - math.exp(-ExtrinsicExchangeRate * TimeElapsed)

				if ProbabilityOfExchange > DeuterationCutoff:
					ExchangedBackboneHydrogens += ExchangeSite(Residue, ExchangeDictionary["Backbone"])

				for Atom in Residue:

					if ExtrinsicExchangeRate != 0.0:
						Atom.set_occupancy(math.log(ExtrinsicExchangeRate))
					else:
						Atom.set_occupancy(0.0)

					Atom.set_bfactor(ProbabilityOfExchange)

			else:

				for Atom in Residue:
					Atom.set_occupancy(0.0)
					Atom.set_bfactor(0.0)

	# Return results
	return ExchangedBackboneHydrogens, ExchangedNonBackboneHydrogens, ExchangedNTerminalHydrogens

#################################################################
# Function used to run exchange algorithm on each protein model #
#################################################################
def ExchangeModels(ProteinStructure, PathToOutputFiles, PathTopdb2pqr, DeuterationCutoff, TimeElapsed, pH):
	InfoPrinted = False

	# Exchange each model
	for ProteinModel in ProteinStructure:

		# Add hydrogens
		NumberOfHydrogens = AddHydrogens(ProteinModel, PathToOutputFiles, PathTopdb2pqr, pH)

		if not InfoPrinted:
			sys.stdout.write("pdb2pqr estimates a total of " + str(NumberOfHydrogens) + " hydrogens in model " + str(ProteinModel.get_id()) + ".\n")
			PrintBreak()

		# Exchange the relevant hydrogens
		ExchangedBackboneHydrogens, ExchangedNonBackboneHydrogens, ExchangedNTerminalHydrogens = DoExchange(ProteinModel, DeuterationCutoff, TimeElapsed)

		if not InfoPrinted:
			sys.stdout.write("Exchanged " + str(ExchangedBackboneHydrogens) + " backbone hydrogens.\n")
			sys.stdout.write("Exchanged " + str(ExchangedNonBackboneHydrogens) + " non-backbone hydrogens.\n")
			sys.stdout.write("Exchanged " + str(ExchangedNTerminalHydrogens) + " N-terminal hydrogens.\n")
			sys.stdout.write("Exchanged " + str(ExchangedBackboneHydrogens + ExchangedNonBackboneHydrogens + ExchangedNTerminalHydrogens) + " hydrogens in total.\n")

			InfoPrinted = True

			PrintBreak()

	return
