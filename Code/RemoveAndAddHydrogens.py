import os
import subprocess
import sys
import Bio.PDB.Atom

########################
# Removal of hydrogens #
########################
def RemoveHydrogens(ProteinModel):
	# Loop over entire structure
	for Chain in list(ProteinModel.get_chains()):

		for Residue in list(Chain.get_residues()):

			if Residue.get_resname() not in ("HOH", "WAT"):
				# And get rid of hydrogens
				for Atom in list(Residue.get_atoms()):

					if Atom.element == "H":
						Residue.detach_child(Atom.get_id())

			else:
				Chain.detach_child(Residue.get_id())

		if len(list(Chain.get_residues())) == 0:
			ProteinModel.detach_child(Chain.get_id())

	return

########################
# Importing .pqr-files #
########################
def ImportPQR(Filename):
	Hydrogens = []

	# Read file
	File = open(Filename)
	Data = File.readlines()
	File.close()

	# Search out all hydrogens
	for Line in Data:

		if Line[0:4] == "ATOM":
			AtomName = Line[12:16].lstrip().rstrip()

			if AtomName[0] == "H":
				ChainID   = Line[21]
				ResidueID = int(Line[22:26])

				x = float(Line[30:38])
				y = float(Line[38:46])
				z = float(Line[46:54])

				Hydrogens.append((AtomName, ResidueID, ChainID, (x, y, z)))

	return Hydrogens

####################
# Adding hydrogens #
####################
def AddHydrogensToModel(ProteinModel, Hydrogens):
	# Loop over hydrogens and add them appropriately
	for Entry in Hydrogens:
		AtomName  = Entry[0]
		ResidueID = Entry[1]
		ChainID   = Entry[2]
		x, y, z   = Entry[3]
		AtomID    = len(list(ProteinModel[ChainID][ResidueID].get_atoms()))

		ProteinModel[ChainID][ResidueID].add(Bio.PDB.Atom.Atom(AtomName, (x, y, z), 1.0, 1.0, " ", AtomName, AtomID, "H"))

	return

#########################
# Addition of hydrogens #
#########################
def AddHydrogens(ProteinModel, PathToOutputFiles, PathTopdb2pqr, pH):
	pD = pH + 0.4

	# Run pdb2pqr on the structure
	subprocess.call("{} {}.prf.pdb --ff=amber --chain --nodebump --ph-calc-method=propka --with-ph={} Temp.pqr".format(PathTopdb2pqr, PathToOutputFiles, pD), shell = True, stdout = open(os.devnull, 'wb'))

	# Get coordinates, names and residues of hydrogens from generated file
	Hydrogens = ImportPQR("Temp.pqr")

	# Add hydrogens to the model
	AddHydrogensToModel(ProteinModel, Hydrogens)

	# And remove the temporary file
	os.remove("Temp.pqr")
	os.remove("Temp.propka")

	return len(Hydrogens)
