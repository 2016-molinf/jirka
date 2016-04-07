# TODO: from SDF import also name of the molecule

from django.core.exceptions import ObjectDoesNotExist
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from django_rdkit import models

class Molecule(models.Model):
    """
    Represents one molecule.
    """

    # fields, which can be calculated on save():
    rdmol = models.MolField()
    internal_id = models.CharField(max_length=32, db_index=True)
    image = models.TextField()
    mw = models.FloatField(db_index=True)
    sum_formula = models.CharField(max_length=32, db_index=True)
    fingerprint = models.CharField(max_length=1024, db_index=True)
    inchi = models.TextField(db_index=True)
    inchi_key = models.CharField(max_length=27, db_index=True)
    name = models.TextField(db_index=True, null=True)
    smiles = models.TextField(db_index=True)
    created = models.DateTimeField(auto_now_add=True)

    # excluded molecules SMILES (they cause rdKit stuck)
    EXCLUDED_MOLECULES = ["C", "CH3", "CH4", "[CH3]", "[C]", "[CH4]"]

    def __str__(self):
        return "Molecule ({id}): '{name}', formula: '{formula}'".format(id=self.internal_id, name=self.name, formula=self.sum_formula)

    def save(self, smiles=None, molfile=None, rdmol=None, inchi=None, update=False, *args, **kwargs):
        if not update:
            if molfile:
                mol = AllChem.MolFromMolBlock(molfile)
            elif smiles:
                mol = AllChem.MolFromSmiles(smiles)
            elif rdmol:
                mol = rdmol
            elif inchi:
                mol = AllChem.MolFromInchi(inchi)

            if mol:
                inchi = AllChem.MolToInchi(mol)
                smiles = AllChem.MolToSmiles(mol)

                if inchi and Molecule.objects.filter(inchi=inchi).count() == 0 and len(inchi) > 1:
                    self.inchi = inchi

                    self.mw = float("{0:.2f}".format(AllChem.CalcExactMolWt(mol)))
                    self.sum_formula = AllChem.CalcMolFormula(mol)
                    self.fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=1024).ToBitString()
                    self.inchi_key = AllChem.InchiToInchiKey(self.inchi)
                    self.molfile = AllChem.MolToMolBlock(mol)
                    self.smiles = smiles
                    self.rdmol = mol

                    # generating SVG image
                    if self.smiles not in self.EXCLUDED_MOLECULES:
                        binMol = AllChem.Mol(self.rdmol.ToBinary())

                        if not binMol.GetNumConformers():
                            rdDepictor.Compute2DCoords(self.rdmol)

                        drawer = rdMolDraw2D.MolDraw2DSVG(100, 100)
                        drawer.DrawMolecule(self.rdmol)
                        drawer.FinishDrawing()
                        svg = drawer.GetDrawingText().replace('svg:', '')

                        # remove first line containg XML meta information
                        self.image_svg = "\n".join(svg.split("\n")[1:]).strip()
                    else:
                        self.image_svg = None

                    super(Molecule, self).save(*args, **kwargs)
                else:
                    raise self.MoleculeExistsInDatabase(smiles)
            else:
                raise self.MoleculeCreationError
        else:
            super(Molecule, self).save(*args, **kwargs)

    class MoleculeExistsInDatabase(Exception):
        def __init__(self, smiles):
            super(Exception, self).__init__(smiles)
            self.smiles = smiles
            self.message = "Cannot add the molecule: it already exists in database."

    class MoleculeCreationError(Exception):
        def __init__(self):
            super(Exception, self).__init__()
            self.message = "Cannot add the molecule: check your structure (valence etc.)."