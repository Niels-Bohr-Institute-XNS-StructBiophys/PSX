###########
# Imports #
###########
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB import DSSP
from Bio.PDB import Selection
from Bio.PDB import NeighborSearch
from Bio.PDB import PDBIO

import sys
import math

from RemoveAndAddHydrogens import *
from IntrinsicExchange import *
from CheckInput import *

####################
# Degree of burial #
####################
def CalculateDegreesOfBurial(ProteinModel, DistanceCutoff):
	# Set up atom search
	ListOfAtoms = Selection.unfold_entities(ProteinModel, "A")
	NeighborSearchOutput = NeighborSearch(ListOfAtoms)

	# Loop over all residues
	DegreesOfBurial = {}

	for Chain in ProteinModel:

		for Residue in Chain:

			if is_aa(Residue.get_resname(), standard = True):

				# Subtract one as the algorithm counts the atom itself
				BackboneNitrogen = Residue["N"]
				DegreesOfBurial[Chain, Residue] = len(NeighborSearchOutput.search(BackboneNitrogen.coord, DistanceCutoff)) - 1

	return DegreesOfBurial

##################
# Hydrogen bonds #
##################
def CalculateHydrogenBonds(ProteinModel, Filename, EnergyCutoff, PathToDSSP):
	# Run DSSP algorithm
	DSSPOutput = DSSP(ProteinModel, Filename, dssp = PathToDSSP)

	# Assign structure
	HydrogenBonds = {}
	TotalNumberOfHydrogenBonds = 0

	for Chain in ProteinModel:
		ChainID = Chain.get_id()
		
		for Residue in Chain:
			ResidueID = Residue.get_id()

			if is_aa(Residue.get_resname(), standard = True):
				HydrogenBonds[Chain, Residue] = 0

				try:
					DSSPEntry = DSSPOutput[(ChainID, ResidueID)]

					if float(DSSPEntry[7]) < EnergyCutoff:
						HydrogenBonds[Chain, Residue] += 1
						TotalNumberOfHydrogenBonds += 1

					if float(DSSPEntry[11]) < EnergyCutoff:
						HydrogenBonds[Chain, Residue] += 1
						TotalNumberOfHydrogenBonds += 1

				except:
					sys.stderr.write("No DSSP entry generated for amino acid residue " + str(ResidueID[1]) + ". Ignoring the residue. \n")

	sys.stdout.write(str(TotalNumberOfHydrogenBonds) + " backbone N-O hydrogen bonds.\n")

	return HydrogenBonds

######################
# Protection factors #
######################
def CalculateProtectionsFactors(DegreesOfBurial, HydrogenBonds, Betac, Betah):
	ProtectionFactors = {}

	for Entry in DegreesOfBurial.keys():
		ProtectionFactors[Entry] = math.exp(Betac * DegreesOfBurial[Entry] + Betah * HydrogenBonds[Entry])

	return ProtectionFactors

###################################################################
# Assign the logarithm of the protection factors to the occupancy #
###################################################################
def AssignProtectionFactors(ProteinModel, ProtectionsFactors):
	# Assign values to the relevant amino acid residues
	for Chain in ProteinModel:
		
		for Residue in Chain:

			if is_aa(Residue.get_resname(), standard = True):
				LogProtectionFactor = math.log(ProtectionsFactors[Chain, Residue])
			
				for Atom in Residue:
					Atom.set_occupancy(LogProtectionFactor)

	# Print range of protection factors
	LogProtectionFactors = [math.log(Value) for Value in ProtectionsFactors.values()]
	LowestLogProtectionFactor  = min(LogProtectionFactors)
	HighestLogProtectionFactor = max(LogProtectionFactors)

	sys.stdout.write("log(Protection factor) range from " + str(LowestLogProtectionFactor) + " to " + str(HighestLogProtectionFactor) + ".\n")

	return

#####################################################################################################
# Function for unwrapping a protein structure into models and running the algorithm on each of them #
#####################################################################################################
def IntrinsicExhangeRatesAndProtectionFactors(ProteinStructure, DistanceCutoff, Temperature, pH, ReferenceData, EnergyCutoff, PathToDSSP, Betac, Betah):
	PrintBreak()

	# Calculate protection factor
	for ProteinModel in ProteinStructure:
		sys.stdout.write("Model " + str(ProteinModel.get_id()) + " contains ")

		# Strip hydrogens from structure if they are there
		RemoveHydrogens(ProteinModel)

		# Output temporary PDB file for DSSP
		TemporaryProteinStructure = Bio.PDB.Structure.Structure("Temporary")
		TemporaryProteinStructure.add(ProteinModel)

		OutputParser = PDBIO()
		OutputParser.set_structure(TemporaryProteinStructure)
		OutputParser.save("Temp.pdb")

		# Get intrinsic exchange rates and assign them to the b-factors
		CalculateExchangeRates(ProteinModel, Temperature, pH, ReferenceData)

		# Calculate protection factors
		DegreesOfBurial    = CalculateDegreesOfBurial(ProteinModel, DistanceCutoff)
		HydrogenBonds      = CalculateHydrogenBonds(ProteinModel, "Temp.pdb", EnergyCutoff, PathToDSSP)
		ProtectionsFactors = CalculateProtectionsFactors(DegreesOfBurial, HydrogenBonds, Betac, Betah)

		# Assign the logarithm of the protection factors to the occupancy
		AssignProtectionFactors(ProteinModel, ProtectionsFactors)

		# Remove temporary PDB file
		os.remove("Temp.pdb")

	PrintBreak()

	return

###########################################################################
# Function for calculating the average protection factor for each residue #
###########################################################################
def AverageProtectionFactors(ProteinStructure):
	# Calculate the average protection factor
	AverageLogProtectionFactors = {}
	
	for ProteinModel in ProteinStructure:

		for Chain in ProteinModel:
			ChainID = Chain.get_id()

			for Residue in Chain:

				if is_aa(Residue.get_resname(), standard = True):
					LogProtectionFactors = [SecondProteinModel[ChainID][Residue.get_id()]["CA"].get_occupancy() for SecondProteinModel in ProteinStructure]

					AverageLogProtectionFactors[(Chain, Residue)] = sum(LogProtectionFactors) / float(len(LogProtectionFactors))

	# And assign them to each residue
	for ProteinModel in ProteinStructure:

		for Chain in ProteinModel:

			for Residue in Chain:

				if is_aa(Residue.get_resname(), standard = True):

					for Atom in Residue:
						Atom.set_occupancy(AverageLogProtectionFactors[(Chain, Residue)])

	return
