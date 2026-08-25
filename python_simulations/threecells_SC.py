## PYTHON INTERFACE
import numpy as np
import RF_Track as rft
import os

# ==========================
# Part 1: Bunch creation
# ==========================
mass = rft.muonmass
Q = 1
charge = 5.3e12

Nx = 32 # Horizontal mesh points, for SC and IBS.
Ny = 32 # Vertical mesh points, for SC and IBS.
Nz = 32 # Longitudinal mesh points, for SC and IBS.
# Ncol = 200 # IBS parameter to limit the maximum number of collisiosn per kick. If there were more, an analytical formula computes the further ones.
SC = rft.SpaceCharge_PIC_FreeSpace(Nx, Ny, Nz) # SC
# IBS = rft.IntraBeamScattering(Nx, Ny, Nz, Ncol) # IBS
rft.SC_engine = SC # Adding Space Charge


M = np.loadtxt('../050625BeamInput_10000_4.04124T.txt')

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
# Part 2: Function
# ==========================
def get_peak_field(cur, r_out, r_in, L):
    mu0 = 4 * np.pi * 1e-7
    return mu0 * cur * (r_out - r_in)

# ==========================
# Example elements
# ==========================
# ---- rft.Solenoid: L_0 ----
L_0 = rft.Solenoid(1, get_peak_field(17000000.0, 0.3, 0.1, 1), 0.1, 0.3, 5)
L_0.set_aperture(0.1)
V1.add(L_0, 0.0, 0.0, 0, 'center')
# ---- rft.Solenoid: MU_0 ----
MU_0 = rft.Solenoid(0.2, get_peak_field(34079073.34, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MU_0.set_aperture(0.1)
V1.add(MU_0, 0.0, 0.0, 1, 'center')
# ---- rft.Solenoid: MU_1 ----
MU_1 = rft.Solenoid(0.2, get_peak_field(3321975.541, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MU_1.set_aperture(0.1)
V1.add(MU_1, 0.0, 0.0, 1.6, 'center')
# ---- rft.Solenoid: MU_2 ----
MU_2 = rft.Solenoid(0.2, get_peak_field(51652979.870000005, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MU_2.set_aperture(0.1)
V1.add(MU_2, 0.0, 0.0, 1.9, 'center')
# ---- rft.Solenoid: HU_0 ----
HU_0 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5)
HU_0.set_aperture(0.03)
V1.add(HU_0, 0.0, 0.0, 2.106, 'center')
# ---- rft.Solenoid: HU_1 ----
HU_1 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5)
HU_1.set_aperture(0.03)
V1.add(HU_1, 0.0, 0.0, 2.118, 'center')
# ---- rft.Solenoid: HU_2 ----
HU_2 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5)
HU_2.set_aperture(0.03)
V1.add(HU_2, 0.0, 0.0, 2.13, 'center')
# ---- rft.Solenoid: H_0 ----
H_0 = rft.Solenoid(1.51736269, get_peak_field(535299999.99999994, 0.09, 0.03, 1.51736269), 0.03, 0.09, 5)
H_0.set_aperture(0.03)
V1.add(H_0, 0.0, 0.0, 2.894681345, 'center')
# ---- rft.Absorber: abs_0 ----
abs_0 = rft.Absorber(1.27736269, 'liquid_hydrogen')
abs_0.set_aperture(0.025)
V1.add(abs_0, 0.0, 0.0, 2.894681345, 'center')

HD_2 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5)
HD_2.set_aperture(0.03)
V1.add(HD_2, 0.0, 0.0, 3.65936269, 'center')

# ---- rft.Solenoid: HD_1 ----
HD_1 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5)
HD_1.set_aperture(0.03)
V1.add(HD_1, 0.0, 0.0, 3.67136269, 'center')
# ---- rft.Solenoid: HD_0 ----
HD_0 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5)
HD_0.set_aperture(0.03)
V1.add(HD_0, 0.0, 0.0, 3.68336269, 'center')
# ---- rft.Solenoid: MD_0 ----
MD_0 = rft.Solenoid(0.2, get_peak_field(29111657.75, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_0.set_aperture(0.1)
V1.add(MD_0, 0.0, 0.0, 3.88936269, 'center')

rf1_1 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
rf1_1.set_phid(-90.0)
rf1_1.set_aperture(0.16)
rf1_1.set_t0(1789.179061)
V1.add(rf1_1, 0.0, 0.0, 4.16436269, 'center')

# ---- rft.Solenoid: MD_1 ----
MD_1 = rft.Solenoid(0.2, get_peak_field(7770668.698, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_1.set_aperture(0.1)
V1.add(MD_1, 0.0, 0.0, 4.43936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_2 ----
rf1_2 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
rf1_2.set_phid(-90.0)
rf1_2.set_aperture(0.16)
rf1_2.set_t0(2454.045840)
V1.add(rf1_2, 0.0, 0.0, 4.71436269, 'center')

# ---- rft.Solenoid: MD_2 ----
MD_2 = rft.Solenoid(0.2, get_peak_field(2078073.5140000002, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_2.set_aperture(0.1)
V1.add(MD_2, 0.0, 0.0, 4.98936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_3 ----
rf1_3 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
rf1_3.set_phid(-90.0)
rf1_3.set_aperture(0.16)
rf1_3.set_t0(120.955517)
V1.add(rf1_3, 0.0, 0.0, 5.26436269, 'center')

# ---- rft.Solenoid: MD_3 ----
MD_3 = rft.Solenoid(0.2, get_peak_field(10499009.54, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_3.set_aperture(0.1)
V1.add(MD_3, 0.0, 0.0, 5.53936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_4 ----
rf1_4 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
rf1_4.set_phid(-90.0)
rf1_4.set_aperture(0.16)
rf1_4.set_t0(785.765173)
V1.add(rf1_4, 0.0, 0.0, 5.81436269, 'center')

# ---- rft.Solenoid: MD_4 ----
MD_4 = rft.Solenoid(0.2, get_peak_field(-7495575.152, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_4.set_aperture(0.1)
V1.add(MD_4, 0.0, 0.0, 6.08936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_5 ----
rf1_5 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
rf1_5.set_phid(-10.0)
rf1_5.set_aperture(0.16)
rf1_5.set_t0(7445.973727)
V1.add(rf1_5, 0.0, 0.0, 6.36436269, 'center')

MD_5 = rft.Solenoid(0.2, get_peak_field(7953739.713, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_5.set_aperture(0.1)
V1.add(MD_5, 0.0, 0.0, 6.63936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_6 ----
rf1_6 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
rf1_6.set_phid(-10.0)
rf1_6.set_aperture(0.16)
rf1_6.set_t0(8098.525725)
V1.add(rf1_6, 0.0, 0.0, 6.91436269, 'center')

# ---- rft.Solenoid: MD_6 ----
MD_6 = rft.Solenoid(0.2, get_peak_field(-11816292.66, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_6.set_aperture(0.1)
V1.add(MD_6, 0.0, 0.0, 7.18936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_7 ----
rf1_7 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
rf1_7.set_phid(-10.0)
rf1_7.set_aperture(0.16)
rf1_7.set_t0(8742.313534)
V1.add(rf1_7, 0.0, 0.0, 7.46436269, 'center')

# ---- rft.Solenoid: MD_7 ----
MD_7 = rft.Solenoid(0.2, get_peak_field(-2474372.499, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_7.set_aperture(0.1)
V1.add(MD_7, 0.0, 0.0, 7.73936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_8 ----
rf1_8 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
rf1_8.set_phid(-10.0)
rf1_8.set_aperture(0.16)
rf1_8.set_t0(9295.232106)
V1.add(rf1_8, 0.0, 0.0, 8.01436269, 'center')

MD_8 = rft.Solenoid(0.2, get_peak_field(-9754001.591, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_8.set_aperture(0.1)
V1.add(MD_8, 0.0, 0.0, 8.28936269, 'center')
# ---- rft.Pillbox_Cavity: rf1_9 ----
rf1_9 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
rf1_9.set_phid(-10.0)
rf1_9.set_aperture(0.16)
V1.add(rf1_9, 0.0, 0.0, 8.56436269, 'center')

# ---- rft.Solenoid: MD_9 ----
MD_9 = rft.Solenoid(0.2, get_peak_field(-38405798.160000004, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_9.set_aperture(0.1)
V1.add(MD_9, 0.0, 0.0, 8.83936269, 'center')
# ---- rft.Solenoid: HU_3 ----
HU_3 = rft.Solenoid(0.012, get_peak_field(-535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5)
HU_3.set_aperture(0.03)
V1.add(HU_3, 0.0, 0.0, 9.04536269, 'center')
# ---- rft.Solenoid: HU_4 ----
HU_4 = rft.Solenoid(0.012, get_peak_field(-535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5)
HU_4.set_aperture(0.03)
V1.add(HU_4, 0.0, 0.0, 9.05736269, 'center')
# ---- rft.Solenoid: HU_5 ----
HU_5 = rft.Solenoid(0.012, get_peak_field(-535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5)
HU_5.set_aperture(0.03)
V1.add(HU_5, 0.0, 0.0, 9.06936269, 'center')
# ---- rft.Solenoid: H_1 ----
H_1 = rft.Solenoid(1.56518469, get_peak_field(-535299999.99999994, 0.09, 0.03, 1.56518469), 0.03, 0.09, 5)
H_1.set_aperture(0.03)
V1.add(H_1, 0.0, 0.0, 9.857955035, 'center')

# ---- rft.Absorber: rft.Abs_1 ----
abs_1 = rft.Absorber(1.32518469, 'liquid_hydrogen')
abs_1.set_aperture(0.025)
V1.add(abs_1, 0.0, 0.0, 9.857955035, 'center')

# ---- rft.Solenoid: HD_5 ----
HD_5 = rft.Solenoid(0.012, get_peak_field(-535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5)
HD_5.set_aperture(0.03)
V1.add(HD_5, 0.0, 0.0, 10.64654738, 'center')

# ---- rft.Solenoid: HD_4 ----
HD_4 = rft.Solenoid(0.012, get_peak_field(-535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5)
HD_4.set_aperture(0.03)
V1.add(HD_4, 0.0, 0.0, 10.65854738, 'center')

# ---- rft.Solenoid: HD_3 ----
HD_3 = rft.Solenoid(0.012, get_peak_field(-535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5)
HD_3.set_aperture(0.03)
V1.add(HD_3, 0.0, 0.0, 10.67054738, 'center')

# ---- rft.Solenoid: MD_10 ----
MD_10 = rft.Solenoid(0.2, get_peak_field(-28752919.330000002, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_10.set_aperture(0.1)
V1.add(MD_10, 0.0, 0.0, 10.87654738, 'center')

# ---- rft.Pillbox_Cavity: rf2_1 ----
rf2_1 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
rf2_1.set_phid(-90.0)
rf2_1.set_aperture(0.16)
V1.add(rf2_1, 0.0, 0.0, 11.15154738, 'center')

# ---- rft.Solenoid: MD_11 ----
MD_11 = rft.Solenoid(0.2, get_peak_field(-5101607.766, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_11.set_aperture(0.1)
V1.add(MD_11, 0.0, 0.0, 11.42654738, 'center')

# ---- rft.Pillbox_Cavity: rf2_2 ----
rf2_2 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
rf2_2.set_phid(-90.0)
rf2_2.set_aperture(0.16)
V1.add(rf2_2, 0.0, 0.0, 11.70154738, 'center')

# ---- rft.Solenoid: MD_12 ----
MD_12 = rft.Solenoid(0.2, get_peak_field(-10181668.370000001, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_12.set_aperture(0.1)
V1.add(MD_12, 0.0, 0.0, 11.97654738, 'center')

# ---- rft.Pillbox_Cavity: rf2_3 ----
rf2_3 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
rf2_3.set_phid(-90.0)
rf2_3.set_aperture(0.16)
V1.add(rf2_3, 0.0, 0.0, 12.25154738, 'center')

# ---- rft.Solenoid: MD_13 ----
MD_13 = rft.Solenoid(0.2, get_peak_field(-5363069.454, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_13.set_aperture(0.1)
V1.add(MD_13, 0.0, 0.0, 12.52654738, 'center')
# ---- rft.Pillbox_Cavity: rf2_4 ----
rf2_4 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
rf2_4.set_phid(-90.0)
rf2_4.set_aperture(0.16)
V1.add(rf2_4, 0.0, 0.0, 12.80154738, 'center')

# ---- rft.Solenoid: MD_14 ----
MD_14 = rft.Solenoid(0.2, get_peak_field(7639264.017, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_14.set_aperture(0.1)
V1.add(MD_14, 0.0, 0.0, 13.07654738, 'center')
# ---- rft.Pillbox_Cavity: rf2_5 ----
rf2_5 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
rf2_5.set_phid(-10.0)
rf2_5.set_aperture(0.16)
V1.add(rf2_5, 0.0, 0.0, 13.35154738, 'center')

# ---- rft.Solenoid: MD_15 ----
MD_15 = rft.Solenoid(0.2, get_peak_field(-7899052.061, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_15.set_aperture(0.1)
V1.add(MD_15, 0.0, 0.0, 13.62654738, 'center')

# ---- rft.Pillbox_Cavity: rf2_6 ----
rf2_6 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
rf2_6.set_phid(-10.0)
rf2_6.set_aperture(0.16)
V1.add(rf2_6, 0.0, 0.0, 13.90154738, 'center')

# ---- rft.Solenoid: MD_16 ----
MD_16 = rft.Solenoid(0.2, get_peak_field(5739264.072, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_16.set_aperture(0.1)
V1.add(MD_16, 0.0, 0.0, 14.17654738, 'center')
# ---- rft.Pillbox_Cavity: rf2_7 ----
rf2_7 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
rf2_7.set_phid(-10.0)
rf2_7.set_aperture(0.16)
V1.add(rf2_7, 0.0, 0.0, 14.45154738, 'center')

# ---- rft.Solenoid: MD_17 ----
MD_17 = rft.Solenoid(0.2, get_peak_field(11250204.07, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_17.set_aperture(0.1)
V1.add(MD_17, 0.0, 0.0, 14.72654738, 'center')
# ---- rft.Pillbox_Cavity: rf2_8 ----
rf2_8 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
rf2_8.set_phid(-10.0)
rf2_8.set_aperture(0.16)
V1.add(rf2_8, 0.0, 0.0, 15.00154738, 'center')

# ---- rft.Solenoid: MD_18 ----
MD_18 = rft.Solenoid(0.2, get_peak_field(5818381.93, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_18.set_aperture(0.1)
V1.add(MD_18, 0.0, 0.0, 15.27654738, 'center')

# ---- rft.Pillbox_Cavity: rf2_9 ----
rf2_9 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
rf2_9.set_phid(-10.0)
rf2_9.set_aperture(0.16)
V1.add(rf2_9, 0.0, 0.0, 15.55154738, 'center')

# ---- rft.Solenoid: MD_19 ----
MD_19 = rft.Solenoid(0.2, get_peak_field(33814699.760000005, 0.3, 0.1, 0.2), 0.1, 0.3, 5)
MD_19.set_aperture(0.1)
V1.add(MD_19, 0.0, 0.0, 15.82654738, 'center')

# ---- rft.Solenoid: HU_6 ----
HU_6 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5)
HU_6.set_aperture(0.03)
V1.add(HU_6, 0.0, 0.0, 16.03254738, 'center')

# ---- rft.Solenoid: HU_7 ----
HU_7 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5)
HU_7.set_aperture(0.03)
V1.add(HU_7, 0.0, 0.0, 16.04454738, 'center')

# ---- rft.Solenoid: HU_8 ----
HU_8 = rft.Solenoid(0.012, get_peak_field(535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5)
HU_8.set_aperture(0.03)
V1.add(HU_8, 0.0, 0.0, 16.05654738, 'center')

# ---- rft.rft.Solenoid: H_2 ----
H_2 = rft.Solenoid(1.34521782, get_peak_field(535299999.99999994, 0.09, 0.03, 1.34521782), 0.03, 0.09, 5)
H_2.set_aperture(0.03)
V1.add(H_2, 0.0, 0.0, 16.73515629, 'center')
# ---- rft.Absorber: abs_2 ----
abs_2 = rft.Absorber(1.10521782, 'liquid_hydrogen')
abs_2.set_aperture(0.025)
V1.add(abs_2, 0.0, 0.0, 16.73515629, 'center')


# ==========================
# Final setup
# ==========================
V1.verbosity = 2
V1.dt_mm = 1
V1.cfx_dt_mm = 10
V1.sc_dt_mm = 10
V1.odeint_algorithm = 'rk2'
V1.odeint_epsabs = 1e-6

C = rft.Drift(0.1)
C.set_aperture(0.0)
V1.add(C, 0, 0, -0.001, 'exit')

V1.set_s0(0.0)
V1.set_s1(17.4)

# Reference particle
M0 = B0.get_phase_space()
B_refp = rft.Bunch6d(mass, charge, Q, np.mean(M0, axis=0))

B1_refp = V1.autophase(B_refp)

rf = V1.get_rf_elements()
for i, cav in enumerate(rf):
    print(f'cavity {i+1} t0 = {cav.get_t0()}')

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

os.makedirs("Phase_space_Bunch6d_SC", exist_ok=True)

for i, B_i in enumerate(S):
    M_i = B_i.get_phase_space()
    info = B_i.get_info()

    other_info = [
        i+1,
        info.S,
        info.transmission,
        info.emitt_4d,
        info.emitt_z,
        info.mean_K,
        info.sigma_t,
        info.emitt_6d
    ]

    np.savetxt(f"Phase_space_Bunch6d_SC/other_info_{i+1:03d}.txt",
               [other_info], delimiter="\t")

    np.savetxt(f"Phase_space_Bunch6d_SC/M_{i+1:03d}.txt",
               M_i, delimiter="\t")

# ==========================
# Transport table
# ==========================
T = L.get_transport_table(
    "%S %mean_K %sigma_x %sigma_y %sigma_t %emitt_x %emitt_y %emitt_4d %emitt_z %N %beta_x %beta_y %alpha_x %alpha_y"
)

np.savetxt("transport_table_3cells_SC.txt", T, fmt="%s")