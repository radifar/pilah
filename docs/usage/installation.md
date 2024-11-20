# Installation

Download the binary release from the [release page](https://github.com/radifar/pilah/releases) then extract it.
In the extracted folder you will find an executable file named with `pilah` and
an example folder which contain configuration and input files for the following
usage guidelines.

You can execute PiLAH directly from the extraction folder as such:

`./pilah --help`

It is highly recommended to copy PiLAH to one of your executable path.
You can check your available executable path by running the command `echo $PATH`.
Some path that you can use are `/home/$USER/.local/bin` or `usr/local/bin`.

After copying PiLAH to your executable path, go to PiLAH's example directory and run
PiLAH again:

`pilah --help`

If it runs properly, the next step is to run PiLAH using one of the provided example.
Still in the example directory run PiLAH with the following command:

`pilah run config_pdb_hux.txt`

Normally it will take between 5-15 seconds for PiLAH to complete the protein ligand extraction and protonation.
After PiLAH finished the processing, PiLAH will generate four files.
One for the isolated protein molecule, `protein_1e66.pdb`.
Second one for the extracted ligand molecule, `HUX.pdb`.
Third one for the ligand image in 2D representation, `HUX.png`.
The last one for the log file which is named with the date and time when PiLAH was executed, such as `log_20240815_172113.txt`

Now PiLAH is ready to be used on your system.