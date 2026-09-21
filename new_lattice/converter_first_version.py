import math
from pathlib import Path
import csv

# ============================================================
# SETTINGS
# ============================================================
INPUT = "lattice.csv"
OUTPUT_PY = "run_muons.py"
BEAM_INPUT = "BeamInput.txt"

def get_peak_field(cur, r_out, r_in, L):
    """Return the on-axis peak field B0 [T].

    Parameters
    ----------
    cur : float
        Current density [A/m^2].
    r_out, r_in : float
        Outer/inner solenoid radii [m].
    L : float
        Solenoid length [m]. Kept in the signature to match the
        original helper, although it is not used by this expression.
    """
    mu0 = 4 * math.pi * 1e-7
    return mu0 * cur * (r_out - r_in)


def fmt(value):
    """Compact, reproducible float formatting for generated Python."""
    return f"{float(value):.15g}"


def solenoid_to_rftrack(name, center, length, rmin, rout, current_density, nsheets):
    """Generate one RF-Track Solenoid using a converter-calculated B0."""
    center = float(center)
    length = float(length)
    rmin = float(rmin)
    rout = float(rout)
    current_density = float(current_density)
    nsheets = int(float(nsheets))

    # Excel: A/mm^2 -> SI: A/m^2
    current_density_am2 = current_density * 1e6
    B0 = get_peak_field(current_density_am2, rout, rmin, length)

    return f"""# ---- rft.Solenoid: {name} ----
{name} = rft.Solenoid({fmt(length)}, {fmt(B0)}, {fmt(rmin)}, {fmt(rout)}, {nsheets})
{name}.set_aperture({fmt(rmin)})
V1.add({name}, 0.0, 0.0, {fmt(center)}, 'center')

"""


def absorber_to_rftrack(name, center, length, density):
    """Generate one variable-density hydrogen absorber.

    RF-Track constructor used here (as agreed for this lattice):
        rft.Absorber(length, 890.5, 1, 1.007, density)
    """
    center = float(center)
    length = float(length)
    density = float(density)

    return f"""# ---- rft.Absorber: {name} ----
{name} = rft.Absorber({fmt(length)}, 890.5, 1, 1.007, {fmt(density)})
V1.add({name}, 0.0, 0.0, {fmt(center)}, 'center')

"""


def cavity_to_rftrack(name, center, num_cav, length_cav, length_drift,
                       gradient, frequency, phase):
    """Generate a group of individual RF-Track Pillbox_Cavity objects.

    Assumptions for this first converter version:
      * center is the centre of the FIRST cavity.
      * length_drift is the empty edge-to-edge drift between cavities.
        Therefore centre-to-centre spacing = length_cav + length_drift.
      * RF-Track phase = input lattice phase - 90 degrees.
    """
    center = float(center)
    num_cav = int(float(num_cav))
    length_cav = float(length_cav)
    length_drift = float(length_drift)
    gradient_vm = float(gradient) * 1e6  # MV/m -> V/m
    frequency_hz = float(frequency) * 1e6  # MHz -> Hz
    phase_rftrack = float(phase) - 90.0

    blocks = []
    for i in range(num_cav):
        cavity_name = f"{name}_{i + 1}"
        z_pos = center + i * (length_cav + length_drift)

        blocks.append(f"""# ---- rft.Pillbox_Cavity: {cavity_name} ----
{cavity_name} = rft.Pillbox_Cavity({fmt(gradient_vm)}, {fmt(frequency_hz)}, {fmt(length_cav)}, 1)
{cavity_name}.set_phid({fmt(phase_rftrack)})
V1.add({cavity_name}, 0.0, 0.0, {fmt(z_pos)}, 'center')

""")

    return "".join(blocks)

# ============================================================
# GENERATED run_muons.py: HEADER
# ============================================================
def make_header():
    return f'''# Auto-generated from {INPUT} by {Path(__file__).name}
# First version: cavity gap length is interpreted as an edge-to-edge drift.

import os
import numpy as np
import RF_Track as rft

# ==========================
# Part 1: Bunch creation
# ==========================
mass = rft.muonmass
Q = 1
charge = 5.3e12

M = np.loadtxt({BEAM_INPUT!r})

x = M[:, 0]
y = M[:, 1]
t = M[:, 6] * rft.ns
Px = M[:, 3]
Py = M[:, 4]
Pz = M[:, 5]

xp = Px / Pz * 1e3
yp = Py / Pz * 1e3
P = np.sqrt(Px**2 + Py**2 + Pz**2)

B0 = rft.Bunch6d(mass, charge, Q, np.column_stack([x, xp, y, yp, t, P]))
B0.set_lifetime(rft.muonlifetime)

L = rft.Lattice()
V1 = rft.Volume()

# ==========================
# Part 2: Lattice elements
# ==========================
'''


def make_final_part(s1):
    return f'''# ==========================
# Final setup
# ==========================
V1.verbosity = 2
V1.dt_mm = 1
V1.cfx_dt_mm = 10
V1.odeint_algorithm = 'rk2'
V1.odeint_epsabs = 1e-6

C = rft.Drift(0.1)
C.set_aperture(0.0)
V1.add(C, 0, 0, -0.001, 'exit')

V1.set_s0(0.0)
V1.set_s1({fmt(s1)})

# Reference particle / RF autophase
M0 = B0.get_phase_space()
B_refp = rft.Bunch6d(mass, charge, Q, np.mean(M0, axis=0))
B1_refp = V1.autophase(B_refp)

rf = V1.get_rf_elements()
for i, cav in enumerate(rf):
    print(f"cavity {{i + 1}} t0 = {{cav.get_t0()}}")

V1.set_tt_nsteps(250)
L.append(V1)

# ==========================
# Tracking
# ==========================
print("Tracking...")

B1 = L.track(B0)
M_lost = L.get_lost_particles()

next_bunch = B1.get_phase_space()
S = L[0].get_bunch_at_tt_screens()

os.makedirs("Phase_space_Bunch6d", exist_ok=True)

for i, B_i in enumerate(S):
    M_i = B_i.get_phase_space()
    info = B_i.get_info()

    other_info = [
        i + 1,
        info.S,
        info.transmission,
        info.emitt_4d,
        info.emitt_z,
        info.mean_K,
        info.sigma_t,
        info.emitt_6d,
    ]

    np.savetxt(
        f"Phase_space_Bunch6d/other_info_{{i + 1:03d}}.txt",
        [other_info],
        delimiter="\\t",
    )

    np.savetxt(
        f"Phase_space_Bunch6d/M_{{i + 1:03d}}.txt",
        M_i,
        delimiter="\\t",
    )

# ==========================
# Transport table
# ==========================
T = L.get_transport_table(
    "%S %mean_K %sigma_x %sigma_y %sigma_t %emitt_x %emitt_y "
    "%emitt_4d %emitt_z %N %beta_x %beta_y %alpha_x %alpha_y"
)

np.savetxt("transport_table.txt", T, fmt="%s")
'''


# ============================================================
# MAIN CONVERSION
# ============================================================
def main():
    rows = []
    with open(INPUT, newline="", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if not row.get("Element") or not row.get("Name"):
                continue
            rows.append(row)

    elements_output = []

    for row in rows:
        element = str(row["Element"]).strip().upper()
        name = str(row["Name"]).strip()

        if element == "SOLENOID":
            elements_output.append(
                solenoid_to_rftrack(
                    name=name,
                    center=row["Centre (m)"],
                    length=row["Length (m)"],
                    rmin=row["Inner radius (m)"],
                    rout=row["Outer radius (m)"],
                    current_density=row["Current density (A/mm^2)"],
                    nsheets=row["nSheets"],
                )
            )

        elif element == "ABSORBER":
            elements_output.append(
                absorber_to_rftrack(
                    name=name,
                    center=row["Centre (m)"],
                    length=row["Length (m)"],
                    density=row["Material density"],
                )
            )

        elif element == "CAVITY":
            elements_output.append(
                cavity_to_rftrack(
                    name=name,
                    center=row["Centre (m)"],
                    num_cav=row["Num cav"],
                    length_cav=row["cavlength"],
                    length_drift=row["gap length"],
                    gradient=row["Gradient (Mv/m)"],
                    frequency=row["Frequency (MHz)"],
                    phase=row["Phase (deg)"],
                )
            )

        else:
            raise ValueError(f"Unknown element type: {element!r} for {name!r}")

    # Small margin after the last physical element.
    s1 = 14

    with open(OUTPUT_PY, "w", encoding="utf-8") as f:
        f.write(make_header())
        for block in elements_output:
            f.write(block)
        f.write(make_final_part(s1))

    print(f"Conversion complete: {INPUT} -> {OUTPUT_PY}")
    print(f"Generated {len(rows)} lattice rows; V1.set_s1({s1:.6g})")


if __name__ == "__main__":
    main()

