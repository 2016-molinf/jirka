from django.db import models
from django.core.exceptions import ObjectDoesNotExist
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

class Molecule(models.Model):
    """
    Represents one molecule.
    """

    # fields, which can be calculated on save():
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

    def save(self, molfile, *args, **kwargs):
        mol = Chem.MolFromMolBlock(molfile)

        if mol:
            smiles = Chem.MolToSmiles(mol)
            if smiles and Molecule.objects.filter(smiles=smiles).count() == 0 and len(smiles) > 1:
                self.smiles = smiles

                try:
                    last_id = Molecule.objects.latest('id').id
                except ObjectDoesNotExist:
                    last_id = 1

                self.internal_id = "MI-J-{}".format(last_id)

                # generating SVG image
                binMol = Chem.Mol(mol.ToBinary())

                if not binMol.GetNumConformers():
                    rdDepictor.Compute2DCoords(mol)

                drawer = rdMolDraw2D.MolDraw2DSVG(250, 250)
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

                super(Molecule, self).save(*args, **kwargs)
            else:
                raise self.MoleculeExistsInDatabase
        else:
            raise self.MoleculeCreationError

    class MoleculeExistsInDatabase(Exception):
        pass

    class MoleculeCreationError(Exception):
        pass