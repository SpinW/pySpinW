from CifFile import ReadCif, CifFile

def load_mcif(filename: str):
    data: CifFile = ReadCif(filename)
    return data

if __name__ == "__main__":
    # data = load_mcif("../../example_structures/1.2_CuSe2O5.mcif")
    data = load_mcif("../../example_structures/1.2_CuSe2O5-extended.mcif")
    # data = load_mcif("../../example_structures/1100231.cif")

