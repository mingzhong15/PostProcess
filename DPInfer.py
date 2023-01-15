from deepmd import DeepPotential
from deepmd.utils.data import DeepmdData

dp = DeepPotential(DIR+'cp.pb')
tmap = dp.get_type_map()

data = DeepmdData(TEST_DIR, type_map=tmap)

numb_test = 100

def _prepare_data(data, numb_test):

    data.add("energy", 1, atomic=False, must=False, high_prec=True)
    data.add("force", 3, atomic=True, must=False, high_prec=False)
    data.add("virial", 9, atomic=False, must=False, high_prec=False)

    if dp.get_dim_fparam() > 0:
        data.add(
            "fparam", dp.get_dim_fparam(), atomic=False, must=True, high_prec=False
        )

    test_data = data.get_test()
    natoms = len(test_data["type"][0])
    nframes = test_data["box"].shape[0]
    numb_test = min(nframes, numb_test)

    coord = test_data["coord"][:numb_test].reshape([numb_test, -1])
    box = test_data["box"][:numb_test]

    if dp.has_efield:
        efield = test_data["efield"][:numb_test].reshape([numb_test, -1])
    else:
        efield = None
    if not data.pbc:
        box = None

    atype = test_data["type"][0]

    if dp.get_dim_fparam() > 0:
        fparam = test_data["fparam"][:numb_test]
    else:
        fparam = None

    return coord, box, atype


coords = coord
cells = box
atom_types = atype
