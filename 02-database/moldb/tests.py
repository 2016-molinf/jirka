from django.test import TestCase
from .models import Molecule
import json

class AddMoleculeTestCase(TestCase):
    def setUp(self):
        self.propane = """Molecule from ChemDoodle Web Components

http://www.ichemlabs.com
  3  2  0  0  0  0            999 V2000
   -0.9330   -0.2500    0.0000 C   0  0  0  0  0  0
    0.0670   -0.2500    0.0000 C   0  0  0  0  0  0
    0.9330    0.2500    0.0000 C   0  0  0  0  0  0
  1  21  0     0  0
  2  31  0     0  0
M  END"""

    def test_add(self):
        mol = Molecule()
        mol.save(self.propane)

        self.assertEqual(mol.internal_id, "MI-J-1")
        self.assertEqual(mol.smiles, "CCC")
        self.assertEqual(mol.sum_formula, "C3H8")
        self.assertEqual(mol.mw, 44.06)
        self.assertEqual(mol.inchi, "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3")
        self.assertEqual(mol.inchi_key, "ATUOYWHBWRKTHZ-UHFFFAOYSA-N")

    def test_add_api(self):
        response = self.client.post(path="/api/addMolecule", data={"molfile": self.propane})
        self.assertEqual(response.status_code, 200)

    def test_converter_api(self):
        response = self.client.post(path="/api/molConverter", data={"data": self.propane,
                                                                "format_from": "molfile",
                                                               "format_to": "smiles"})
        self.assertEqual(response.status_code, 200)
        data = json.loads(str(response.content, encoding="utf-8"))
        self.assertEqual(data["data"], "CCC")

        response = self.client.post(path="/api/molConverter", data={"data": "NotValidMolfile",
                                                                "format_from": "molfile",
                                                                "format_to": "smiles"})
        self.assertEqual(response.status_code, 200)
        data = json.loads(str(response.content, encoding="utf-8"))
        self.assertTrue(data["error"])

