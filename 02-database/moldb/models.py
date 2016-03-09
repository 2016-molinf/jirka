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

    def __str__(self):
        return "Molecule ({id}): '{name}', formula: '{formula}'".format(id=self.internal_id, name=self.name, formula=self.sum_formula)

    def save(self, smiles=None, molfile=None, rdmol=None, *args, **kwargs):
        if molfile:
            mol = Chem.MolFromMolBlock(molfile)
        elif smiles:
            mol = Chem.MolFromSmiles(smiles)
        elif rdmol:
            mol = rdmol

        if mol:
            smiles = Chem.MolToSmiles(mol)

            if smiles and Molecule.objects.filter(smiles=smiles).count() == 0 and len(smiles) > 1:
                self.smiles = smiles

                try:
                    last_id = Molecule.objects.latest('id').id + 1
                except ObjectDoesNotExist:
                    last_id = 1

                self.internal_id = "MI-J-{}".format(last_id)

                # generating SVG image
                binMol = Chem.Mol(mol.ToBinary())

                if not binMol.GetNumConformers():
                    rdDepictor.Compute2DCoords(mol)

                drawer = rdMolDraw2D.MolDraw2DSVG(100, 100)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText().replace('svg:', '')
                # remove first line containg XML meta information
                self.image = "\n".join(svg.split("\n")[1:]).strip()

                self.mw = float("{0:.2f}".format(AllChem.CalcExactMolWt(mol)))
                self.sum_formula = AllChem.CalcMolFormula(mol)
                self.fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 4, nBits=1024).ToBitString()
                self.inchi = Chem.MolToInchi(mol)
                self.inchi_key = AllChem.InchiToInchiKey(self.inchi)

                self.rdmol = mol

                super(Molecule, self).save(*args, **kwargs)
            else:
                raise self.MoleculeExistsInDatabase(smiles)
        else:
            raise self.MoleculeCreationError

    class MoleculeExistsInDatabase(Exception):
        def __init__(self, smiles):
            super(Exception, self).__init__(smiles)
            self.smiles = smiles
            self.message = "Cannot add the molecule: it already exists in database."

    class MoleculeCreationError(Exception):
        def __init__(self):
            super(Exception, self).__init__()
            self.message = "Cannot add the molecule: check your structure (valence etc.)."