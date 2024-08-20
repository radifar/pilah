pyinstaller src/pilah/pilah.py \
--name "pilah" \
--onefile \
--add-data "src/pilah/dimorphite_dl/site_substructures.smarts:pilah/dimorphite_dl" \
--add-data "src/pilah/pKAI/models/*.pt:pilah/pKAI/models" \
--add-data "src/pilah/meeko/data/*.*:pilah/meeko/data" \
--add-data "src/pilah/meeko/data/*.*:meeko/data" \
--clean