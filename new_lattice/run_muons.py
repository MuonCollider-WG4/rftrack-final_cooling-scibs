# Auto-generated from lattice.csv by converter_first_version.py
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

M = np.loadtxt('BeamInput.txt')

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
# ---- rft.Solenoid: L1 ----
L1 = rft.Solenoid(0.5, 5.02654824574367, 0.1, 0.3, 10)
L1.set_aperture(0.1)
V1.add(L1, 0.0, 0.0, 0, 'center')

# ---- rft.Solenoid: MU1 ----
MU1 = rft.Solenoid(0.1, 5.48362033528598, 0.1, 0.3, 10)
MU1.set_aperture(0.1)
V1.add(MU1, 0.0, 0.0, 0.4, 'center')

# ---- rft.Solenoid: MU2 ----
MU2 = rft.Solenoid(0.1, 10.0872186362005, 0.1, 0.3, 10)
MU2.set_aperture(0.1)
V1.add(MU2, 0.0, 0.0, 0.6, 'center')

# ---- rft.Solenoid: MU3 ----
MU3 = rft.Solenoid(0.1, 1.78551790198511, 0.1, 0.3, 10)
MU3.set_aperture(0.1)
V1.add(MU3, 0.0, 0.0, 0.8, 'center')

# ---- rft.Solenoid: MU4 ----
MU4 = rft.Solenoid(0.1, 21.6929409223974, 0.1, 0.3, 10)
MU4.set_aperture(0.1)
V1.add(MU4, 0.0, 0.0, 1, 'center')

# ---- rft.Solenoid: H1corrU1 ----
H1corrU1 = rft.Solenoid(0.012, 70.6311709935979, 0.03, 0.135, 10)
H1corrU1.set_aperture(0.03)
V1.add(H1corrU1, 0.0, 0.0, 1.106, 'center')

# ---- rft.Solenoid: H1corrU2 ----
H1corrU2 = rft.Solenoid(0.012, 63.9043928037314, 0.03, 0.125, 10)
H1corrU2.set_aperture(0.03)
V1.add(H1corrU2, 0.0, 0.0, 1.118, 'center')

# ---- rft.Solenoid: H1corrU3 ----
H1corrU3 = rft.Solenoid(0.012, 57.177614613865, 0.03, 0.115, 10)
H1corrU3.set_aperture(0.03)
V1.add(H1corrU3, 0.0, 0.0, 1.13, 'center')

# ---- rft.Solenoid: H1 ----
H1 = rft.Solenoid(1.228, 40.3606691391988, 0.03, 0.09, 10)
H1.set_aperture(0.03)
V1.add(H1, 0.0, 0.0, 1.75, 'center')

# ---- rft.Absorber: ABS1 ----
ABS1 = rft.Absorber(1.1914, 890.5, 1, 1.007, 0.048492925)
V1.add(ABS1, 0.0, 0.0, 1.75, 'center')

# ---- rft.Solenoid: H1corrD3 ----
H1corrD3 = rft.Solenoid(0.012, 57.177614613865, 0.03, 0.115, 10)
H1corrD3.set_aperture(0.03)
V1.add(H1corrD3, 0.0, 0.0, 2.37, 'center')

# ---- rft.Solenoid: H1corrD2 ----
H1corrD2 = rft.Solenoid(0.012, 63.9043928037314, 0.03, 0.125, 10)
H1corrD2.set_aperture(0.03)
V1.add(H1corrD2, 0.0, 0.0, 2.382, 'center')

# ---- rft.Solenoid: H1corrD1 ----
H1corrD1 = rft.Solenoid(0.012, 70.6311709935979, 0.03, 0.135, 10)
H1corrD1.set_aperture(0.03)
V1.add(H1corrD1, 0.0, 0.0, 2.394, 'center')

# ---- rft.Solenoid: MD1 ----
MD1 = rft.Solenoid(0.1, 9.11242374145182, 0.1, 0.3, 10)
MD1.set_aperture(0.1)
V1.add(MD1, 0.0, 0.0, 2.6, 'center')

# ---- rft.Solenoid: MD2 ----
MD2 = rft.Solenoid(0.1, 0.405175263722511, 0.1, 0.3, 10)
MD2.set_aperture(0.1)
V1.add(MD2, 0.0, 0.0, 2.8, 'center')

# ---- rft.Solenoid: MD3 ----
MD3 = rft.Solenoid(0.1, 7.62926908339193, 0.1, 0.3, 10)
MD3.set_aperture(0.1)
V1.add(MD3, 0.0, 0.0, 3, 'center')

# ---- rft.Solenoid: MD4 ----
MD4 = rft.Solenoid(0.1, -2.17870733753954, 0.1, 0.3, 10)
MD4.set_aperture(0.1)
V1.add(MD4, 0.0, 0.0, 3.2, 'center')

# ---- rft.Solenoid: MD5 ----
MD5 = rft.Solenoid(0.1, -7.43052212784349, 0.1, 0.3, 10)
MD5.set_aperture(0.1)
V1.add(MD5, 0.0, 0.0, 3.4, 'center')

# ---- rft.Solenoid: MD6 ----
MD6 = rft.Solenoid(0.1, -3.15980258627361, 0.1, 0.3, 10)
MD6.set_aperture(0.1)
V1.add(MD6, 0.0, 0.0, 3.6, 'center')

# ---- rft.Solenoid: L2 ----
L2 = rft.Solenoid(0.5, -5.02654824574367, 0.1, 0.3, 10)
L2.set_aperture(0.1)
V1.add(L2, 0.0, 0.0, 4.05, 'center')

# ---- rft.Solenoid: MU5 ----
MU5 = rft.Solenoid(0.1, -6.25120176522286, 0.1, 0.3, 10)
MU5.set_aperture(0.1)
V1.add(MU5, 0.0, 0.0, 4.45, 'center')

# ---- rft.Solenoid: MU6 ----
MU6 = rft.Solenoid(0.1, -4.3045446499543, 0.1, 0.3, 10)
MU6.set_aperture(0.1)
V1.add(MU6, 0.0, 0.0, 4.65, 'center')

# ---- rft.Solenoid: MU7 ----
MU7 = rft.Solenoid(0.1, -4.59186590251255, 0.1, 0.3, 10)
MU7.set_aperture(0.1)
V1.add(MU7, 0.0, 0.0, 4.85, 'center')

# ---- rft.Solenoid: MU8 ----
MU8 = rft.Solenoid(0.1, -10.1163472219041, 0.1, 0.3, 10)
MU8.set_aperture(0.1)
V1.add(MU8, 0.0, 0.0, 5.05, 'center')

# ---- rft.Solenoid: H2corrU1 ----
H2corrU1 = rft.Solenoid(0.012, -70.6311709935979, 0.03, 0.135, 10)
H2corrU1.set_aperture(0.03)
V1.add(H2corrU1, 0.0, 0.0, 5.156, 'center')

# ---- rft.Solenoid: H2corrU2 ----
H2corrU2 = rft.Solenoid(0.012, -63.9043928037314, 0.03, 0.125, 10)
H2corrU2.set_aperture(0.03)
V1.add(H2corrU2, 0.0, 0.0, 5.168, 'center')

# ---- rft.Solenoid: H2corrU3 ----
H2corrU3 = rft.Solenoid(0.012, -57.177614613865, 0.03, 0.115, 10)
H2corrU3.set_aperture(0.03)
V1.add(H2corrU3, 0.0, 0.0, 5.18, 'center')

# ---- rft.Solenoid: H2 ----
H2 = rft.Solenoid(0.428, -40.3606691391988, 0.03, 0.09, 10)
H2.set_aperture(0.03)
V1.add(H2, 0.0, 0.0, 5.4, 'center')

# ---- rft.Absorber: ABS2 ----
ABS2 = rft.Absorber(0.355140399, 890.5, 1, 1.007, 0.039871555)
V1.add(ABS2, 0.0, 0.0, 5.4, 'center')

# ---- rft.Solenoid: H2corrD3 ----
H2corrD3 = rft.Solenoid(0.012, -57.177614613865, 0.03, 0.115, 10)
H2corrD3.set_aperture(0.03)
V1.add(H2corrD3, 0.0, 0.0, 5.62, 'center')

# ---- rft.Solenoid: H2corrD2 ----
H2corrD2 = rft.Solenoid(0.012, -63.9043928037314, 0.03, 0.125, 10)
H2corrD2.set_aperture(0.03)
V1.add(H2corrD2, 0.0, 0.0, 5.632, 'center')

# ---- rft.Solenoid: H2corrD1 ----
H2corrD1 = rft.Solenoid(0.012, -70.6311709935979, 0.03, 0.135, 10)
H2corrD1.set_aperture(0.03)
V1.add(H2corrD1, 0.0, 0.0, 5.644, 'center')

# ---- rft.Solenoid: MD7 ----
MD7 = rft.Solenoid(0.1, -10.0146464083098, 0.2, 0.4, 10)
MD7.set_aperture(0.2)
V1.add(MD7, 0.0, 0.0, 5.89, 'center')

# ---- rft.Solenoid: MD8 ----
MD8 = rft.Solenoid(0.1, -4.10626566353427, 0.2, 0.4, 10)
MD8.set_aperture(0.2)
V1.add(MD8, 0.0, 0.0, 6.42, 'center')

# ---- rft.Solenoid: MD9 ----
MD9 = rft.Solenoid(0.1, -4.78548709831925, 0.2, 0.4, 10)
MD9.set_aperture(0.2)
V1.add(MD9, 0.0, 0.0, 6.95, 'center')

# ---- rft.Pillbox_Cavity: Acc_1_1 ----
Acc_1_1 = rft.Pillbox_Cavity(7000000, 70000000, 0.25, 1)
Acc_1_1.set_phid(-52)
V1.add(Acc_1_1, 0.0, 0.0, 6.95, 'center')

# ---- rft.Pillbox_Cavity: Acc_1_2 ----
Acc_1_2 = rft.Pillbox_Cavity(7000000, 70000000, 0.25, 1)
Acc_1_2.set_phid(-52)
V1.add(Acc_1_2, 0.0, 0.0, 7.48, 'center')

# ---- rft.Pillbox_Cavity: Acc_1_3 ----
Acc_1_3 = rft.Pillbox_Cavity(7000000, 70000000, 0.25, 1)
Acc_1_3.set_phid(-52)
V1.add(Acc_1_3, 0.0, 0.0, 8.01, 'center')

# ---- rft.Pillbox_Cavity: Acc_1_4 ----
Acc_1_4 = rft.Pillbox_Cavity(7000000, 70000000, 0.25, 1)
Acc_1_4.set_phid(-52)
V1.add(Acc_1_4, 0.0, 0.0, 8.54, 'center')

# ---- rft.Solenoid: MD10 ----
MD10 = rft.Solenoid(0.1, -3.06832543251915, 0.2, 0.4, 10)
MD10.set_aperture(0.2)
V1.add(MD10, 0.0, 0.0, 7.48, 'center')

# ---- rft.Solenoid: MD11 ----
MD11 = rft.Solenoid(0.1, -5.17226960098055, 0.2, 0.4, 10)
MD11.set_aperture(0.2)
V1.add(MD11, 0.0, 0.0, 8.01, 'center')

# ---- rft.Solenoid: MD12 ----
MD12 = rft.Solenoid(0.1, 2.55345377846888, 0.2, 0.4, 10)
MD12.set_aperture(0.2)
V1.add(MD12, 0.0, 0.0, 8.54, 'center')

# ---- rft.Pillbox_Cavity: Acc_1_1_1 ----
Acc_1_1_1 = rft.Pillbox_Cavity(2500000, 70000000, 0.25, 1)
Acc_1_1_1.set_phid(-90)
V1.add(Acc_1_1_1, 0.0, 0.0, 8.805, 'center')

# ---- rft.Pillbox_Cavity: Acc_1_1_2 ----
Acc_1_1_2 = rft.Pillbox_Cavity(2500000, 70000000, 0.25, 1)
Acc_1_1_2.set_phid(-90)
V1.add(Acc_1_1_2, 0.0, 0.0, 9.335, 'center')

# ---- rft.Pillbox_Cavity: Acc_1_1_3 ----
Acc_1_1_3 = rft.Pillbox_Cavity(2500000, 70000000, 0.25, 1)
Acc_1_1_3.set_phid(-90)
V1.add(Acc_1_1_3, 0.0, 0.0, 9.865, 'center')

# ---- rft.Solenoid: MD13 ----
MD13 = rft.Solenoid(0.1, 7.99846058482129, 0.2, 0.4, 10)
MD13.set_aperture(0.2)
V1.add(MD13, 0.0, 0.0, 9.07, 'center')

# ---- rft.Solenoid: L3 ----
L3 = rft.Solenoid(0.5, 5.02654824574367, 0.1, 0.3, 10)
L3.set_aperture(0.1)
V1.add(L3, 0.0, 0.0, 9.85, 'center')

# ---- rft.Solenoid: MU9 ----
MU9 = rft.Solenoid(0.1, 6.3339509186211, 0.1, 0.3, 10)
MU9.set_aperture(0.1)
V1.add(MU9, 0.0, 0.0, 10.25, 'center')

# ---- rft.Solenoid: MU10 ----
MU10 = rft.Solenoid(0.1, 3.16892574620689, 0.1, 0.3, 10)
MU10.set_aperture(0.1)
V1.add(MU10, 0.0, 0.0, 10.45, 'center')

# ---- rft.Solenoid: MU11 ----
MU11 = rft.Solenoid(0.1, 5.12119961444984, 0.1, 0.3, 10)
MU11.set_aperture(0.1)
V1.add(MU11, 0.0, 0.0, 10.65, 'center')

# ---- rft.Solenoid: MU12 ----
MU12 = rft.Solenoid(0.1, 7.84396745804606, 0.1, 0.3, 10)
MU12.set_aperture(0.1)
V1.add(MU12, 0.0, 0.0, 10.85, 'center')

# ---- rft.Solenoid: H3corrU1 ----
H3corrU1 = rft.Solenoid(0.012, 70.6311709935979, 0.03, 0.135, 10)
H3corrU1.set_aperture(0.03)
V1.add(H3corrU1, 0.0, 0.0, 10.956, 'center')

# ---- rft.Solenoid: H3corrU2 ----
H3corrU2 = rft.Solenoid(0.012, 63.9043928037314, 0.03, 0.125, 10)
H3corrU2.set_aperture(0.03)
V1.add(H3corrU2, 0.0, 0.0, 10.968, 'center')

# ---- rft.Solenoid: H3corrU3 ----
H3corrU3 = rft.Solenoid(0.012, 57.177614613865, 0.03, 0.115, 10)
H3corrU3.set_aperture(0.03)
V1.add(H3corrU3, 0.0, 0.0, 10.98, 'center')

# ---- rft.Solenoid: H3 ----
H3 = rft.Solenoid(0.328, 40.3606691391988, 0.03, 0.09, 10)
H3.set_aperture(0.03)
V1.add(H3, 0.0, 0.0, 11.15, 'center')

# ---- rft.Absorber: ABS3 ----
ABS3 = rft.Absorber(0.3, 890.5, 1, 1.007, 0.032892196)
V1.add(ABS3, 0.0, 0.0, 11.15, 'center')

# ---- rft.Solenoid: H3corrD3 ----
H3corrD3 = rft.Solenoid(0.012, 57.177614613865, 0.03, 0.115, 10)
H3corrD3.set_aperture(0.03)
V1.add(H3corrD3, 0.0, 0.0, 11.32, 'center')

# ---- rft.Solenoid: H3corrD2 ----
H3corrD2 = rft.Solenoid(0.012, 63.9043928037314, 0.03, 0.125, 10)
H3corrD2.set_aperture(0.03)
V1.add(H3corrD2, 0.0, 0.0, 11.332, 'center')

# ---- rft.Solenoid: H3corrD1 ----
H3corrD1 = rft.Solenoid(0.012, 70.6311709935979, 0.03, 0.135, 10)
H3corrD1.set_aperture(0.03)
V1.add(H3corrD1, 0.0, 0.0, 11.344, 'center')

# ---- rft.Solenoid: MD14 ----
MD14 = rft.Solenoid(0.1, 6.463538418697, 0.1, 0.3, 10)
MD14.set_aperture(0.1)
V1.add(MD14, 0.0, 0.0, 11.55, 'center')

# ---- rft.Solenoid: MD15 ----
MD15 = rft.Solenoid(0.1, -0.0934921149851344, 0.1, 0.3, 10)
MD15.set_aperture(0.1)
V1.add(MD15, 0.0, 0.0, 11.75, 'center')

# ---- rft.Solenoid: MD16 ----
MD16 = rft.Solenoid(0.1, 11.2617059902758, 0.1, 0.3, 10)
MD16.set_aperture(0.1)
V1.add(MD16, 0.0, 0.0, 11.95, 'center')

# ---- rft.Solenoid: MD17 ----
MD17 = rft.Solenoid(0.1, -10.4326038245717, 0.1, 0.3, 10)
MD17.set_aperture(0.1)
V1.add(MD17, 0.0, 0.0, 12.15, 'center')

# ==========================
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
V1.set_s1(13.0)

# Reference particle / RF autophase
M0 = B0.get_phase_space()
B_refp = rft.Bunch6d(mass, charge, Q, np.mean(M0, axis=0))
B1_refp = V1.autophase(B_refp)

rf = V1.get_rf_elements()
for i, cav in enumerate(rf):
    print(f"cavity {i + 1} t0 = {cav.get_t0()}")

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
        f"Phase_space_Bunch6d/other_info_{i + 1:03d}.txt",
        [other_info],
        delimiter="\t",
    )

    np.savetxt(
        f"Phase_space_Bunch6d/M_{i + 1:03d}.txt",
        M_i,
        delimiter="\t",
    )

# ==========================
# Transport table
# ==========================
T = L.get_transport_table(
    "%S %mean_K %sigma_x %sigma_y %sigma_t %emitt_x %emitt_y "
    "%emitt_4d %emitt_z %N %beta_x %beta_y %alpha_x %alpha_y"
)

np.savetxt("transport_table.txt", T, fmt="%s")
