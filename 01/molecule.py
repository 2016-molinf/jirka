import re

class Molecule(object):
    __pattern_formula = re.compile(r"([A-Z][a-z]|[A-Z])(\d+)?")

    def __init__(self, sum_formula):
        if re.match(self.__pattern_formula, sum_formula):
            self.__sum_formula = sum_formula
            self.__initStdMolWeights()
        else:
            raise Exception("Incorrect summary formula format!")

    def __initStdMolWeights(self):
        self.__stdMolWeights = {}

        with open("atom_weights.txt", mode="r", encoding="utf-8") as f:
            for line in f:
                cols = line.split()
                if cols[3] != "-":
                    self.__stdMolWeights[cols[1]] = float(re.sub(r"\(\d+\)", "",
                                                            cols[3].replace("[", "").replace(",", "")))

    def atom_count(self, atom_symbol=None):
        if atom_symbol is not None:
            pattern = re.compile(r"({})(\d+)?".format(atom_symbol))
            try:
                count = re.findall(pattern, self.__sum_formula)[0][1]
                if count == "":
                    return 1
                else:
                    return int(re.findall(pattern, self.__sum_formula)[0][1])
            except:
                return 0

        count = 0
        for s in re.findall(self.__pattern_formula, self.__sum_formula):
            if s[1] == "":
                count += 1
            else:
                count += int(s[1])

        return count

    def mol_weight(self):
        weight = 0
        for s in re.findall(self.__pattern_formula, self.__sum_formula):
            if s[1] == "":
                weight += self.__stdMolWeights[s[0]]
            else:
                weight += self.__stdMolWeights[s[0]] * int(s[1])

        return weight

    def sum_formula(self):
        return self.__sum_formula

    def getStdMolWeights(self):
        return self.__stdMolWeights

mol = Molecule("H2O")

assert mol.atom_count() == 3
assert mol.atom_count("H") == 2
assert mol.atom_count("O") == 1
assert -1e-1 < mol.mol_weight() - 18 < 1e-1
assert mol.sum_formula() == "H2O"