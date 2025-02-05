pyinstaller src/pilah/pilah.py \
--name "pilah" \
--onefile \
--add-data="${CONDA_PREFIX}/share/openbabel/*/*:Openbabel/data/" \
--add-binary="${CONDA_PREFIX}/lib/openbabel/*/mol2format.so:Openbabel/" \
--add-binary="${CONDA_PREFIX}/lib/openbabel/*/pdbformat.so:Openbabel/" \
--runtime-hook=pyi_rth_obdata.py \
--add-data "src/pilah/dimorphite_dl/site_substructures.smarts:pilah/dimorphite_dl" \
--add-data "src/pilah/pKAI/models/*.pt:pilah/pKAI/models" \
--add-data "src/pilah/meeko/data/*.*:pilah/meeko/data" \
--add-data "src/pilah/meeko/data/*.*:meeko/data" \
--clean

rm -rf dist/examples/
cp -r examples/ dist/