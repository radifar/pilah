import os
import sys


pyinstaller_temp_dir = sys._MEIPASS
ob_temp_dir = os.path.join(pyinstaller_temp_dir, "Openbabel")

os.environ["BABEL_DATADIR"] = os.path.join(ob_temp_dir, 'data')
os.environ["BABEL_LIBDIR"] = ob_temp_dir
os.environ["PATH"] = ob_temp_dir + ";" + os.environ["PATH"]