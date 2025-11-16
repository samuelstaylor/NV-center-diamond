import math

# generate_occupations.py
def write_occupations(
    nbnd=152,
    nup_occ=128,  # number of occupied bands for spin-up
    ndown_occ=126,  # number of occupied bands for spin-down
    excitations_down=None,  # list of (remove, add) tuples for ΔSCF excitations
    filename="occupations.out",
):
    """
    Generate Quantum ESPRESSO OCCUPATIONS card for spin-polarized calculation.

    Parameters
    ----------
    nbnd : int
        Total number of bands.
    nup_occ : int
        Number of occupied bands in the spin-up channel.
    ndown_occ : int
        Number of occupied bands in the spin-down channel.
    excitations_down : list of tuples
        Each tuple defines an excitation as (empty_band, filled_band)
        Example: [(1022, 1023)] will empty band 1022 and fill band 1023 for spin down.
    filename : str
        Output filename.
    """

    print(f"Generating: '{filename}' with the following parameters:")
    print(f"- nbnd={nbnd}\n- nup_occ={nup_occ}\n- ndown_occ={ndown_occ}\n- excitations_down={excitations_down}")
    
    # initialize occupations
    occ_up = [1.0 if i < nup_occ else 0.0 for i in range(nbnd)]
    occ_down = [1.0 if i < ndown_occ else 0.0 for i in range(nbnd)]

    # apply ΔSCF excitations to spin-down
    if excitations_down:
        for empty_band, filled_band in excitations_down:
            # convert from 1-indexed to 0-indexed
            empty_band_idx = empty_band - 1
            filled_band_idx = filled_band - 1
            if empty_band_idx < nbnd and filled_band_idx < nbnd:
                occ_down[empty_band_idx] = 0.0
                occ_down[filled_band_idx] = 1.0
            else:
                print(f"Warning: excitation indices ({empty_band}, {filled_band}) exceed nbnd={nbnd}")

    # write to file
    with open(filename, "w") as f:
        f.write("OCCUPATIONS\n")

        # --- Write spin-up occupations ---
        for i in range(0, len(occ_up), 10):
            line = " ".join(str(int(x)) for x in occ_up[i : i + 10])
            f.write(line + "\n")

        # --- Add a blank line between spin channels ---
        f.write("\n")

        # --- Write spin-down occupations ---
        for i in range(0, len(occ_down), 10):
            line = " ".join(str(int(x)) for x in occ_down[i : i + 10])
            f.write(line + "\n")

    print(f"- Occupations file written to {filename}")
    print(f"   Spin-up: {sum(occ_up)} electrons")
    print(f"   Spin-down: {sum(occ_down)} electrons")
    print(f"   Total: {sum(occ_up) + sum(occ_down)} electrons\n")

    # --- Print a concise summary of occupied/unoccupied bands ---
    def summarize_bands(occ_list, label):
        occupied = [i + 1 for i, occ in enumerate(occ_list) if occ == 1.0]
        unoccupied = [i + 1 for i, occ in enumerate(occ_list) if occ == 0.0]

        # Helper to compactly represent ranges (e.g., [1–128,130])
        def compact_ranges(indices):
            if not indices:
                return "none"
            ranges = []
            start = prev = indices[0]
            for i in indices[1:]:
                if i == prev + 1:
                    prev = i
                else:
                    ranges.append(f"{start}–{prev}" if start != prev else f"{start}")
                    start = prev = i
            ranges.append(f"{start}–{prev}" if start != prev else f"{start}")
            return ",".join(ranges)

        occ_str = compact_ranges(occupied)
        unocc_str = compact_ranges(unoccupied)
        print(f"{label} occupied bands [{occ_str}] | unoccupied bands [{unocc_str}]")

    summarize_bands(occ_up, "Spin-up")
    summarize_bands(occ_down, "Spin-down")


# Example usage:
# 1040 total bands, 1024 up, 1022 down
# Excitation: spin-down electron from band 1022 → 1023
if __name__ == "__main__":
    
    nat = 63
    # for NV^- center in diamond
    #            n_carbon*4 + valence electrons in N-
    n_electrons = (nat - 1) * 4 + 6
    homo = math.ceil(n_electrons / 2)
    n_excited_states = 25
    nbnd = homo + n_excited_states
    nup_occ = homo + 1
    ndown_occ = homo - 1

    write_occupations(
        nbnd=nbnd,
        nup_occ=nup_occ,
        ndown_occ=ndown_occ,
        # [(i,j)] = excitation from the i-th band to the j-th band.
        excitations_down=[(ndown_occ, ndown_occ + 1)],
    )
