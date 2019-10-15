###########
# Imports #
###########
from Bio.PDB.Polypeptide import three_to_one
from Bio.PDB.Polypeptide import is_aa

import sys

###############################
# Function for checking input #
###############################
def CheckRDFormat(ReferenceData):

	if ReferenceData in ("poly", "oligo"):
		return ReferenceData

	else:
	    raise argparse.ArgumentTypeError("{} is not a valid argument for the reference data.\n".format(ReferenceData))

###############################################################################
# Function for printing breaks between the different segments of the software #
###############################################################################
def PrintBreak():
	sys.stdout.write("********************************************************************************\n")

	return

#################################
# Function for input validation #
#################################
def CheckInput(ProteinStructure):
	Sequences = {}
	SequencePrinted = False

	# Loop over all models in the structure
	for ProteinModel in ProteinStructure:
		NumberOfUnknownResidues = 0
		ModelID = ProteinModel.get_id()
		Sequences[ModelID] = {}

		if not SequencePrinted:
			sys.stdout.write("Model " + str(ModelID) + " contains the following chains:\n")

		# And loop over all chains in the given model
		for Chain in ProteinModel:
			ChainID = Chain.get_id()
			Sequences[ModelID][ChainID] = []

			# As well as all residues in the chain
			for Residue in Chain:
				ResidueName = Residue.get_resname()

				if is_aa(ResidueName, standard = True):
					Sequences[ModelID][ChainID].append(three_to_one(ResidueName))

				else:
					Sequences[ModelID][ChainID].append("?")
					NumberOfUnknownResidues += 1

			# And print the chain
			if not SequencePrinted:
				sys.stdout.write("Chain " + Chain.get_id() + ":\n" + str("".join(Sequences[ModelID][ChainID])) + "\n")

		if not SequencePrinted:
			SequencePrinted = True

		sys.stdout.write("Model " + str(ModelID) + " contains " + str(NumberOfUnknownResidues) + " unknown residues.\n")

	# Compare sequences to each other
	for Model1 in Sequences.keys():

		for Model2 in Sequences.keys():

			if Sequences[Model1] != Sequences[Model2]:
				sys.stderr.write("The sequences in models " +  str(Model1) + " and " + str(Model2) + " are different. Unfortunately, this prevents the software from continuing.\n")
				sys.exit()
