import os
from xml.etree import ElementTree as ET
import json

script_dir = os.path.dirname(os.path.abspath(__file__))

Hartree2eV = 27.2114

# FROM ANALYZING THE OPTIMIZED GROUND STATE IN PART 2
fname = os.path.join(script_dir, 'relax/pwscf.save/data-file-schema.xml')

xmlData = ET.parse(fname)
root = xmlData.getroot()
gs_etot = float(root.find('output').find('total_energy').find('etot').text) * Hartree2eV

print(f'Ground state total energy: {gs_etot:.6f} eV')


# FROM ANALYZING THE OPTIMIZED GROUND STATE IN PART 2
fname = os.path.join(script_dir, 'dscf/pwscf.save/data-file-schema.xml')
xmlData = ET.parse(fname)
root = xmlData.getroot()
dscf_etot = float(root.find('output').find('total_energy').find('etot').text) * Hartree2eV

print(f'ΔSCF state total energy: {dscf_etot:.6f} eV')


# NOTE:
# - total energy of the optimized excited state can be decomposed into two components:
#   1.) ground state total energy at this optimized excited state geometry
#   2.) excitation energy at the optimized excited state geometry

# -- LET'S FIND (1.) FIRST
fname = os.path.join(script_dir, 'tddft/final_geo/pwscf.save/data-file-schema.xml')

xmlData = ET.parse(fname)
root = xmlData.getroot()
# THIS IS (1.) THE GROUND STATE TOTAL ENERGY AT OPTIMIZED EXCITED STATE GEOMETRY
es_etot1 = float(root.find('output').find('total_energy').find('etot').text) * Hartree2eV


# -- NOW LET'S FIND (2.)
fname = os.path.join(script_dir, 'tddft/final_geo/west.wbse.save/wbse.json')

with open(fname, 'r') as f:
    j = json.load(f)

which_es = int(j['input']['wbse_control']['forces_state'])
# THIS IS (2.) THE EXCITATION ENERGY AT THE OPTIMIZED GEOMETRY
es_etot2 = float(j['exec']['davitr'][-1]['ev'][which_es - 1]) * Hartree2eV / 2.0

print(f"GS total energy @ optimized ES geo: {es_etot1:.6f} eV")
print(f"EXCITATION ENERGY from GS to ES @ optimized ES geo: {es_etot2:.6f} eV")
es_etot = es_etot1 + es_etot2
print(f'Excited state total energy: {es_etot:.6f} eV')

# NOTE: The adiabatic excitation energy is the difference between the total energies
# of the excited and ground states at their respective optimized geometries:
aee = es_etot - gs_etot
dscf_dif = dscf_etot - gs_etot
print(f'Adiabatic excitation energy (TDDFT predicted ZPL): {aee:.6f} eV')
print(f'(ΔSCF energy) - (GS total energy): {dscf_dif:.6f} eV')