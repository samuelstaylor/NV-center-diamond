#!/usr/bin/env python3
"""
merge_qe_header.py

Merge QE header (template) into phonopy-generated supercell-*.in files.

- Expects pw-gs-dft.in (template) in the current directory.
- Expects supercell-*.in produced by phonopy in the current directory.
- Produces run-supercell-XXX.in files (ready to feed to pw.x).
"""

from pathlib import Path
import re

# filenames (change if you used a different template name)
TEMPLATE = Path("es-dft-pw.in") # NOTE: CHANGED THIS TO BE FILE NAME
SUPERCELL_GLOB = "displacements/supercell-*.in"
OUT_PREFIX = ""

if not TEMPLATE.exists():
    raise SystemExit(f"Template {TEMPLATE} not found. Place your pw-gs-dft.in in this folder.")

# read template
tmpl_text = TEMPLATE.read_text().splitlines()

# find start of ATOMIC_SPECIES in template (case-insensitive)
idx_as = None
for i, ln in enumerate(tmpl_text):
    if ln.strip().upper().startswith("ATOMIC_SPECIES"):
        idx_as = i
        break

if idx_as is None:
    raise SystemExit("Could not find 'ATOMIC_SPECIES' in template. Ensure pw-gs-dft.in is well-formed.")

# header is everything before ATOMIC_SPECIES (we will let phonopy-written file supply ATOMIC_SPECIES/CELL/ATOMIC_POSITIONS)
header_lines = tmpl_text[:idx_as]

# modify header: ensure calculation = 'scf' and tprnfor = .TRUE. present
header_text = "\n".join(header_lines)

# Replace calculation namelist setting if present in &CONTROL
# common forms: calculation = 'relax' or calculation='relax'
header_text = re.sub(r"(calculation\s*=\s*)'[^']+'", r"\1'scf'", header_text, flags=re.IGNORECASE)

# Ensure tprnfor is present in &CONTROL: if not, add after &CONTROL line block
if "tprnfor" not in header_text.lower():
    # insert tprnfor = .TRUE. just after &CONTROL line if &CONTROL exists
    m = re.search(r"(&CONTROL\b.*?/)", header_text, flags=re.IGNORECASE | re.DOTALL)
    if m:
        # insert tprnfor inside the control block before the closing /
        block = m.group(1)
        block_new = re.sub(r"(/\s*$)", "   tprnfor = .TRUE.\n\\1", block, flags=re.MULTILINE)
        header_text = header_text.replace(block, block_new)
    else:
        # fallback: append a minimal &CONTROL block
        header_text = "&CONTROL\n   calculation = 'scf'\n   tprnfor = .TRUE.\n/\n" + header_text

# Ensure we have &ELECTRONS and &SYSTEM blocks present is assumed; otherwise user should fix template.

# find all supercell files
sc_files = sorted(Path(".").glob(SUPERCELL_GLOB))
if not sc_files:
    raise SystemExit("No files matching supercell-*.in found in current directory.")

# default K_POINTS block to append if supercell file lacks K_POINTS
kpoints_block = "K_POINTS gamma"

created = []
for sc in sc_files:
    body = sc.read_text()  # this contains ATOMIC_SPECIES, CELL_PARAMETERS, ATOMIC_POSITIONS, maybe K_POINTS
    # Check if `K_POINTS` already present in body
    if re.search(r"^\s*K_POINTS\b", body, flags=re.IGNORECASE | re.MULTILINE):
        merged = header_text + "\n" + body
    else:
        merged = header_text + "\n" + body + "\n" + kpoints_block

    out_name = f"{OUT_PREFIX}{sc.name}"
    Path(out_name).write_text(merged)
    created.append(out_name)

print(f"Created {len(created)} merged input files (prefix='{OUT_PREFIX}'):")
for p in created[:3]:
    print(" ", p)
if len(created) > 10:
    print("  ...")
    print(" ", created[-1])
print("Done. Review one file (e.g. supercell-001.in) before launching jobs.")
