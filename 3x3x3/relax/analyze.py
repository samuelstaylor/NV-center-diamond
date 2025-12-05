from xml.etree import ElementTree as ET

Hartree2eV = 27.2114

fname = 'pwscf.save/data-file-schema.xml'

xmlData = ET.parse(fname)
root = xmlData.getroot()
gs_etot = float(root.find('output').find('total_energy').find('etot').text) * Hartree2eV

print(f'Ground state total energy: {gs_etot:.6f} eV')