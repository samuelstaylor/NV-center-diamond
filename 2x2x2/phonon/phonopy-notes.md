NOTE: gs-dft-pw.in uses the generated supercell geometry. other input params referenced from pypl tutorial:

generates: 

supercell-001.in … supercell-078.in

Each file is:

→ one supercell where one atom has been displaced by +δ or −δ along x, y, or z.

Specifically:

Phonopy picks a symmetry-inequivalent atom.

Applies a displacement like:

+0.01 Å in x

−0.01 Å in x

+0.01 Å in y

−0.01 Å in y

+0.01 Å in z

−0.01 Å in z

But because symmetry reduces the required displacements, one atom may require only ±x, or only +x/+y, etc.

Why 78 files?

Phonopy reduced the number of displacements from the full 63 atoms × 3 directions × 2 signs = 378 possibilities.

Your structure has a LOT of symmetry (almost cubic, periodic C–N network), so symmetry reduction drastically cuts down the number of unique displacements.

Thus: 78 supercells must be computed.

Each one of these should be run with pw.x, then the forces extracted for phonopy.


Here’s *exactly* what that command does and what **each flag** means:

---

# ✅ Full Meaning of the Command

```bash
phonopy --qe -d --dim="1 1 1" -c gs-dft-pw.in -v
```

This command tells **phonopy** to:

**Generate supercells with finite displacements for Quantum ESPRESSO**,
using the QE input file `gs-dft-pw.in`,
with a supercell size of 1×1×1,
and print verbose information.

---

# 📌 Breaking down each flag

### **1. `--qe`**

Tells phonopy that your input file (`-c gs-dft-pw.in`) is a **Quantum ESPRESSO PWscf input**.

Phonopy will parse:

* CELL_PARAMETERS
* ATOMIC_POSITIONS
* ATOMIC_SPECIES
* ibrav
* nat / ntyp
* etc.

And it will create QE-format displaced supercell files.

Without this flag, phonopy assumes VASP by default.

---

### **2. `-d`**

Short for **`--displacement`**.

This tells phonopy to:

👉 generate finite-displacement supercells
👉 AND create `phonopy_disp.yaml` describing the displacement pattern

Without `-d`, phonopy would not produce the `supercell-*` files.

---

### **3. `--dim="1 1 1"`**

This sets the **supercell multiplication matrix** (the supercell size):

```
a' = 1 × a
b' = 1 × b
c' = 1 × c
```

Meaning:

* no enlarged supercell
* you are doing phonons in the *primitive* cell
* phonopy will generate displacements **only inside the primitive cell**

If you used `--dim="2 2 2"`, you'd generate an 8× larger supercell for more accurate phonons.

---

### **4. `-c gs-dft-pw.in`**

This tells phonopy which file contains the **primitive cell structure**.

Equivalent to:

```
--cell gs-dft-pw.in
```

Phonopy reads:

* structure
* atomic species
* positions
* (optionally) magnetization info

It **does not** read &CONTROL, &ELECTRONS, etc for physics — only the structure matters.

---

### **5. `-v`**

Verbose mode.

Outputs more information to the terminal, such as:

* symmetry operations found
* displacement patterns
* atomic mapping between symmetry groups
* generated filenames

---

# 🔍 Summary Table

| Flag              | Meaning                                              |
| ----------------- | ---------------------------------------------------- |
| `--qe`            | Read & write QE PWscf format                         |
| `-d`              | Generate displaced supercells                        |
| `--dim="1 1 1"`   | Set supercell replication (here primitive cell only) |
| `-c gs-dft-pw.in` | Input QE structure file                              |
| `-v`              | Verbose output                                       |

---

# 🧠 What the command *produces*

You should get:

* `phonopy_disp.yaml` → displacement map
* `supercell.in` → undisplaced primitive supercell (1×1×1)
* `supercell-001.in`, `supercell-002.in`, … → displaced structures

You already saw this.

---

If you'd like, I can also tell you **how many displacements you *should* get theoretically** for a 63-atom primitive cell with C₃ᵥ-like symmetry.
