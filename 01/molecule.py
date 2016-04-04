import re

ATOM_WEIGHTS_FILE = "atom_weights.txt"

def _initStdMolWeights():
    stdMolWeights = {}

    with open(ATOM_WEIGHTS_FILE, mode="r", encoding="utf-8") as f:
        for line in f:
            cols = line.split()
            if cols[3] != "-":
                stdMolWeights[cols[1]] = float(re.sub(r"\(\d+\)", "",
                                                        cols[3].replace("[", "")
                                                        .replace(",", "")))
    return stdMolWeights

class Molecule(object):
    __pattern_sum_formula = re.compile(r"([A-Z][a-z]*)(\d*)")
    __stdMolWeights = _initStdMolWeights()

    def __init__(self, sum_formula):
        self.__init_atoms(sum_formula)

    def __init_atoms(self, sum_formula):
        self.__atoms = {}
        self.__atoms_count = 0

        sum_formula_control = ""

        # iterate over atom groups ("atom", "count")
        for group in re.findall(self.__pattern_sum_formula, sum_formula):
            if group[0] not in self.__atoms.keys():
                self.__atoms[group[0]] = {"count": 0}

            if group[1] == "":
                self.__atoms[group[0]]["count"] = 1
            else:
                self.__atoms[group[0]]["count"] += int(group[1])

            sum_formula_control += group[0] + group[1]

        # check summary formula
        if sum_formula == sum_formula_control:
            self.__sum_formula = sum_formula
            self.__atoms_count = sum([atom["count"] for atom in self.__atoms.values()])

            # calculate molecular weight
            self.__mol_weight = 0.0
            for atom in self.__atoms.keys():
                self.__mol_weight += self.__stdMolWeights[atom] * self.__atoms[atom]["count"]
        else:
            raise self.SummaryFormulaException("You have error in you summary formula. Check it! Allowed format is for example 'HClNa2F3'.")

    class SummaryFormulaException(Exception):
        pass

    def getAtomCount(self, atom_symbol=None):
        if atom_symbol is not None:
            if atom_symbol in self.__atoms.keys():
                return self.__atoms[atom_symbol]["count"]
            else:
                return 0
        else:
            return self.__atoms_count

    def getMolWeight(self):
        return self.__mol_weight

    def getSumFormula(self):
        return self.__sum_formula

    def getStdMolWeights(self):
        return self.__stdMolWeights

    def getAtoms(self):
        return self.__atoms

mol = Molecule("H2O")

assert mol.getAtomCount() == 3
assert mol.getAtomCount("H") == 2
assert mol.getAtomCount("O") == 1
assert -1e-1 < mol.getMolWeight() - 18 < 1e-1
assert mol.getSumFormula() == "H2O"