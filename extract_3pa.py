#!/usr/bin/python3

import os
import sys
import glob
import math
import csv
import re
from scipy.constants import fine_structure, c, physical_constants


assert len(sys.argv) == 3, "ext_3pa_to_table.py OUTPUT.csv WORKDIR_PATH"
OUTFILE = sys.argv[1]
WORKDIR_PATH = sys.argv[2]
FILES_PATH = WORKDIR_PATH + '*.out'
OUTFILE_PATH = WORKDIR_PATH + OUTFILE


# Math constants
PI = math.pi
HART = 0.0367493 # (1 eV ~ 0.036.. Hartree)
GAMMA = float(0.00367493) # in Hartrees
LIGHTVEL = c*100 # in cm^-1
BOHRAD = physical_constants['Bohr radius'][0]*100 # in cm


def get_outfiles(DIR):
	"""Gets the file paths from which the 3PA info will be extracted.

	:param DIR: The directory of the output files.
	:return: ``file_paths``, a list of file paths.
	"""

	file_paths = glob.glob(DIR, recursive=True)
	return file_paths


def get_3pa(outfiles):
	"""Gets the 3PA data for each excited state for each file.

	Note: just collects linear polarization.
	:param outfiles: A list of file paths.
	:return: ``dict_3pa``, a dictionary of 3pa features.
	"""

	sigma_a = (4*(PI**3) * (fine_structure*BOHRAD**8)) / (3*(LIGHTVEL**2))

	rex_bas = 'Basis set:'
	rex_3pa = '\d+\s+Linear\s+\d+'
	dict_3pa = []
	for file in outfiles:
		file_data = []
		bas = None
		with open(f'{file}', 'r') as file_out:
			for line in file_out:
				if re.search(rex_bas, line):
					bas = line.split()
					if len(bas) > 2:
						bas = bas[2]
					else:
						bas = None
				if re.search(rex_3pa, line):
					file_data.append(line.split())


		for state in file_data:
			# Excitation energy in Hartree
			omega_h = '{:.3f}'.format(float(state[2])*HART)
			# Excitation energy in eV
			omega_ev = '{:.2f}'.format(float(state[2]))
			sigma_b = (((float(omega_h)/3)**3)/GAMMA)
			deltaf = '{:.0f}'.format(float(state[4]))
			deltag = '{:.0f}'.format(float(state[5]))
			delta3p = float(state[6])
			sigma3pa = (float(sigma_a)*sigma_b*float(delta3p))
			dict_3pa.append(
				{
					'File': file,
					'Basis Set': bas,
					'State (S$\\_{n}$)': int(state[1]),
					'$\\omega$ (Hartree)': float(omega_h),
					'$\\omega$ (eV)': float(omega_ev),
					'$\\delta_{f}$ (au)': float(deltaf),
					'$\\delta_{g}$ (au)': float(deltag),
					'$\\delta_{3PA}$ (au)': float(delta3p),
					'$\\sigma^{3PA}$ (NA)': float(sigma3pa),
				}
			)

	return dict_3pa


def get_cvsfile(OUTFILE_NAMEPATH, dict_3pa):
	"""Writes a file listing the 3PA features found in each file.

	:param OUTFILE_NAMEPATH: ``WORKDIR_PATH`` + the name of the csv file to be returned.
	:type OUTFILE_NAME: str
	:param dict_3pa: 3PA features.
	:type dict_3pa: list of dictionaries.
	:return: a .csv file (``OUTFILE_NAMEPATH``) listing 3PA features, including exc. energies and 3PA cross-sections.
	"""

	def head():
		head_print = "\nThe 3PA equation employed here:\n"
		head_print += "\n       (4 * pi**3 * alpha * a_0**8 * (omega/3)**3 * delta)\n"
		head_print += "sigma =  --------------------------------------------------- \n"
		head_print += "                      (3 * c_0**2 * gamma) \n"
		head_print += "\nin units of cm^6 * s^2 * photon^-2 \n"
		head_print += "\n"
		head_print += "The value of the constants used in this script: \n"
		head_print += f"Light speed (c_0): {LIGHTVEL} \n"
		head_print += f"Bohr radius (a_0): {BOHRAD} \n"
		head_print += f"Fine structure constant (alpha): {fine_structure} \n"
		print(head_print)

	try:
		with open(OUTFILE_NAMEPATH, mode='w') as outfile:
			writer = csv.DictWriter(outfile, dict_3pa[0].keys(), lineterminator='\n',)
			writer.writeheader()
			writer.writerows(dict_3pa)

		head()
		print('\nFinished extracting 3PA data.\n')

	except IndexError:
		print(f'\nNo 3PA info to extract from {WORKDIR_PATH}\n')



# Get files list.
outfiles = get_outfiles(FILES_PATH)

# Get 3PA features.
data_3pa = get_3pa(outfiles)

# Write csv file.
get_cvsfile(OUTFILE_PATH, data_3pa)

