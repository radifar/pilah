from rdkit import Chem
import pytest

from pilah.protonate import (
    process_ligand,
    get_atoms_from_pattern,
    process_protein,
    ionizable_AA_Smarts,
    AA_modifier
)

from data import added_records, ionization_data, pkai_data
from helpers import assert_protonated_mol

config = dict(
        input = "tests/data/6hsh.pdb",
        ligand_id = "GOK",
        include_metal = "yes",
        protein_out = "protein_6hsh.pdb",
        ligand_out = "GOK.pdb",
        ligand_image = "GOK.png",
        image_size = "large",
        ligand_smiles = "Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO",
        ph = 7.4
    )

def test_get_atoms_from_pattern(ionizable_peptide):
    # only test ionizable moiety pattern
    # only ionizable amino acids
    
    carboxylate_indices = get_atoms_from_pattern(
        ionizable_peptide,
        ionizable_AA_Smarts["carboxylate"]
    )
    imidazole_indices = get_atoms_from_pattern(
        ionizable_peptide,
        ionizable_AA_Smarts["imidazole"]
    )
    thiol_indices = get_atoms_from_pattern(
        ionizable_peptide,
        ionizable_AA_Smarts["thiol"]
    )
    phenol_indices = get_atoms_from_pattern(
        ionizable_peptide,
        ionizable_AA_Smarts["phenol"]
    )
    guanidinium_indices = get_atoms_from_pattern(
        ionizable_peptide,
        ionizable_AA_Smarts["guanidinium"]
    )
    amonium_indices = get_atoms_from_pattern(
        ionizable_peptide,
        ionizable_AA_Smarts["amonium"]
    )

    assert len(carboxylate_indices) == 3
    assert len(imidazole_indices) == 1
    assert len(thiol_indices) == 1
    assert len(phenol_indices) == 1
    assert len(guanidinium_indices) == 2  # the pattern match the guanidinium twice
    assert len(amonium_indices) == 2

@pytest.mark.parametrize(
        "rdkit_mol, non_organic_num",
        [
            ("cr_apo_mbs", 1),    # 1 CR
            ("zb_ADAR1", 3),      # 2 CD 1 NI
            ("methanobactin", 6), # 6 CU
            ("calmodulin", 5)     # 1 NA 1 MG 2 ZN 1 AS
        ]
)
def test_get_atoms_from_pattern_non_organic(rdkit_mol, non_organic_num, request):
    rdkit_mol = request.getfixturevalue(rdkit_mol)
    non_organic_indices = get_atoms_from_pattern(
        rdkit_mol,
        ionizable_AA_Smarts["non-organic"]
    )

    assert len(non_organic_indices) == non_organic_num

# Use fixture for AA_modifier object to keep the object
# pristine for each test
@pytest.fixture
def ionizable_AA_1j3f(cr_apo_mbs):
    return AA_modifier(cr_apo_mbs, pkai_data.result_1j3f)

@pytest.fixture
def ionizable_AA_1xmk(zb_ADAR1):
    return AA_modifier(zb_ADAR1, pkai_data.result_1xmk)

@pytest.fixture
def ionizable_AA_2xji(methanobactin):
    return AA_modifier(methanobactin, pkai_data.result_2xji)

@pytest.fixture
def ionizable_AA_3ucy(calmodulin):
    return AA_modifier(calmodulin, pkai_data.result_3ucy)

@pytest.mark.parametrize(
        "ionizable_AA, pdb_id, i_records_length",
        [
            ("ionizable_AA_1j3f", "1j3f", 58),
            ("ionizable_AA_1xmk", "1xmk", 31),
            ("ionizable_AA_2xji", "2xji", 24),
            ("ionizable_AA_3ucy", "3ucy", 26)
        ]
)
def test_AA_modifier_init(ionizable_AA, pdb_id, i_records_length, request):
    ionizable_AA = request.getfixturevalue(ionizable_AA)
    i_records = ionizable_AA.ionization_records

    added_records_pdb = added_records.added_records[pdb_id]

    assert added_records_pdb <= i_records.keys()
    assert len(i_records) == i_records_length

    for key in added_records_pdb:
        record = i_records[key]
        assert record[0] == 14.0
        assert record[2] == "Positive"

@pytest.mark.parametrize(
        "pH, pT, pHpT",
        [
            (7.4, 1.0, "pH7.4"),
            (3.0, 1.0, "pH3"),
            (11.0, 1.0, "pH11"),
            (3.0, 2.0, "pH3pT2"),
            (11.0, 2.0, "pH11pT2")
        ]
)
@pytest.mark.parametrize(
        "ionizable_AA, pdb_id",
        [
            ("ionizable_AA_1j3f", "1j3f"),
            ("ionizable_AA_1xmk", "1xmk"),
            ("ionizable_AA_2xji", "2xji"),
            ("ionizable_AA_3ucy", "3ucy")
        ]
)
def test_AA_modifier_ionize(pH, pT, pHpT, ionizable_AA, pdb_id, request):
    ionizable_AA = request.getfixturevalue(ionizable_AA)
    i_records_charge = ionization_data.i_records_charge[pdb_id]

    ionizable_AA.ionize_aa(pH, pT)

    records = ionizable_AA.ionization_records
    for key, value in i_records_charge.items():
        assert value[pHpT] == records[key][2]

@pytest.mark.parametrize(
        "pH, pT",
        [
            (7.4, 1.0),
            (3.0, 1.0),
            (11.0, 1.0),
            (3.0, 2.0),
            (11.0, 2.0,)
        ]
)
@pytest.mark.parametrize(
        "ionizable_AA, pdb_id",
        [
            ("ionizable_AA_1j3f", "1j3f"),
            ("ionizable_AA_1xmk", "1xmk"),
            ("ionizable_AA_2xji", "2xji"),
            ("ionizable_AA_3ucy", "3ucy")
        ]
)
def test_AA_modifier_get_protonated_mol(pH, pT, ionizable_AA, pdb_id, request):
    ionizable_AA = request.getfixturevalue(ionizable_AA)
    i_records_ids = ionization_data.i_records_charge[pdb_id].keys()

    ionizable_AA.ionize_aa(pH, pT)

    whole_i_records = ionizable_AA.ionization_records
    protonated_mol = ionizable_AA.get_protonated_mol()
    
    assert_protonated_mol(protonated_mol, whole_i_records, i_records_ids)

def test_AA_modifier_get_protonated_mol_pdbqt(ionizable_AA_2xji):
    i_records_ids = ionization_data.i_records_charge["2xji_pdbqt"].keys()
    ionizable_AA_2xji.no_tyr_ionization = True
    ionizable_AA_2xji.ionize_aa(11.0, 1.0)

    whole_i_records = ionizable_AA_2xji.ionization_records
    protonated_mol = ionizable_AA_2xji.get_protonated_mol()
    
    assert_protonated_mol(protonated_mol, whole_i_records, i_records_ids)

@pytest.mark.parametrize(
        "ligand_pdb_block, num_atoms",
        [
            ("pdb_block_ligand_HUX", 21),
            ("pdb_block_ligand_G39", 20),
            ("pdb_block_ligand_GOK", 29),
        ]
)
def test_process_ligand_no_smiles(ligand_pdb_block, num_atoms, request):
    data = {}
    
    ligand_pdb_block = request.getfixturevalue(ligand_pdb_block)
    ligand, ligand_Hs = process_ligand(data, ligand_pdb_block)

    assert ligand.GetNumAtoms() == num_atoms
    assert ligand.GetNumAtoms() == ligand.GetNumHeavyAtoms()
    assert ligand_Hs.GetNumAtoms() > ligand_Hs.GetNumHeavyAtoms()

@pytest.mark.parametrize(
        "ligand_pdb_block, smiles, dimorphite_smiles",
        # SMILES generated using dimorphite-dl in pilah/src/pilah/dimorphite_dl directory
        # using the following command:
        # python dimorphite_dl.py --smiles "CCC1=C[C@@H]2Cc3c(c(c4ccc(cc4n3)Cl)N)[C@@H](C2)C1" --min_ph 7.4 --max_ph 7.4 --pka_precision 0
        # python dimorphite_dl.py --smiles "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O" --min_ph 7.4 --max_ph 7.4 --pka_precision 0
        # python dimorphite_dl.py --smiles "Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO" --min_ph 7.4 --max_ph 7.4 --pka_precision 0
        [
            ("pdb_block_ligand_HUX",
             "CCC1=C[C@@H]2Cc3c(c(c4ccc(cc4n3)Cl)N)[C@@H](C2)C1",
             "CCC1=C[C@@H]2Cc3nc4cc(Cl)ccc4c(N)c3[C@H](C1)C2"),
            ("pdb_block_ligand_G39",
             "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)O",
             "CCC(CC)O[C@@H]1C=C(C(=O)[O-])C[C@H]([NH3+])[C@H]1NC(C)=O"),
            ("pdb_block_ligand_GOK",
             "Cn1cc(c2c1cccc2)CNCC3CCN(CC3)c4ncc(cn4)C(=O)NO",
             "Cn1cc(C[NH2+]CC2CCN(c3ncc(C(=O)NO)cn3)CC2)c2ccccc21"),
        ]
)
def test_process_ligand_with_smiles(ligand_pdb_block, smiles, dimorphite_smiles, request):
    data = {"ligand_smiles": smiles}
    ligand_pdb_block = request.getfixturevalue(ligand_pdb_block)

    ligand, ligand_Hs = process_ligand(data, ligand_pdb_block)

    assert ligand_Hs.GetNumAtoms() > ligand_Hs.GetNumHeavyAtoms()
    assert Chem.MolToSmiles(ligand) == dimorphite_smiles

@pytest.mark.parametrize(
        "protein_pdb_block",
        [
            ("pdb_block_protein_1e66"),
            ("pdb_block_protein_5nzn"),
            ("pdb_block_protein_6hsh"),
        ]
)
def test_process_protein(protein_pdb_block, request):
    data = {"protein_out": "protein.mol2"}
    protein_pdb_block = request.getfixturevalue(protein_pdb_block)

    ionized_mol, whole_i_records = process_protein(data, protein_pdb_block)
    ids = whole_i_records.keys()

    assert_protonated_mol(ionized_mol, whole_i_records, ids)
