import numpy as np
import matplotlib.pyplot as plt
import os

# ==========================
# Paths
# ==========================
folder_0 = "Phase_space_Bunch6d"
folder_SC = "Phase_space_Bunch6d_SC"
folder_SC_IBS = "Phase_space_Bunch6d_SC_IBS"

cut = 2.5
n_steps = 249

# ==========================
# Preallocate
# ==========================
mean_t0 = np.full(n_steps, np.nan)
mean_t_SC = np.full(n_steps, np.nan)
mean_t_SC_IBS = np.full(n_steps, np.nan)

emitt4d = np.full(n_steps, np.nan)
emitt4d_SC = np.full(n_steps, np.nan)
emitt4d_SC_IBS = np.full(n_steps, np.nan)

emittz = np.full(n_steps, np.nan)
emittz_SC = np.full(n_steps, np.nan)
emittz_SC_IBS = np.full(n_steps, np.nan)

ES_arr = np.full(n_steps, np.nan)
ES_SC_arr = np.full(n_steps, np.nan)
ES_SC_IBS_arr = np.full(n_steps, np.nan)

length_arr = np.full(n_steps, np.nan)
length_SC_arr = np.full(n_steps, np.nan)
length_SC_IBS_arr = np.full(n_steps, np.nan)

mean_K_arr = np.full(n_steps, np.nan)
mean_K_SC_arr = np.full(n_steps, np.nan)
mean_K_SC_IBS_arr = np.full(n_steps, np.nan)

sigma_t_arr = np.full(n_steps, np.nan)
sigma_t_SC_arr = np.full(n_steps, np.nan)
sigma_t_SC_IBS_arr = np.full(n_steps, np.nan)

N_arr = np.full(n_steps, np.nan)
N_SC_arr = np.full(n_steps, np.nan)
N_SC_IBS_arr = np.full(n_steps, np.nan)

struct_M = [None]*n_steps
struct_M_SC = [None]*n_steps
struct_M_SC_IBS = [None]*n_steps

allX, allY = [], []

mass = 105.658  # MeV/c^2 if muon (adjust if RF_Track gives different)

# ==========================
# Main loop
# ==========================
for i in range(n_steps):

    f0 = f"{folder_0}/M_{i+1:03d}.txt"
    f1 = f"{folder_SC}/M_{i+1:03d}.txt"
    f2 = f"{folder_SC_IBS}/M_{i+1:03d}.txt"

    o0 = f"{folder_0}/other_info_{i+1:03d}.txt"
    o1 = f"{folder_SC}/other_info_{i+1:03d}.txt"
    o2 = f"{folder_SC_IBS}/other_info_{i+1:03d}.txt"

    if not all(os.path.exists(f) for f in [f0, f1, f2, o0, o1, o2]):
        print(f"Skipping i={i} (missing file)")
        continue

    info = np.loadtxt(o0)
    info_SC = np.loadtxt(o1)
    info_IBS = np.loadtxt(o2)

    length_arr[i] = info[1]
    N_arr[i] = info[2]
    mean_K_arr[i] = info[5]
    sigma_t_arr[i] = info[6]

    length_SC_arr[i] = info_SC[1]
    N_SC_arr[i] = info_SC[2]
    mean_K_SC_arr[i] = info_SC[5]
    sigma_t_SC_arr[i] = info_SC[6]

    length_SC_IBS_arr[i] = info_IBS[1]
    N_SC_IBS_arr[i] = info_IBS[2]
    mean_K_SC_IBS_arr[i] = info_IBS[5]
    sigma_t_SC_IBS_arr[i] = info_IBS[6]

    N0 = 5.3e12

    # ==========================
    # Load phase space
    # ==========================
    M = np.loadtxt(f0)
    M_SC = np.loadtxt(f1)
    M_IBS = np.loadtxt(f2)

    struct_M[i] = M
    struct_M_SC[i] = M_SC
    struct_M_SC_IBS[i] = M_IBS

    allX.extend([M[:,0], M_SC[:,0], M_IBS[:,0]])
    allY.extend([M[:,1], M_SC[:,1], M_IBS[:,1]])

    # ==========================
    # Physics variables
    # ==========================
    x, xp = M[:,0], M[:,1]
    y, yp = M[:,2], M[:,3]
    t = M[:,4] * 3.335641e-9
    P = M[:,5]

    x_SC, xp_SC = M_SC[:,0], M_SC[:,1]
    y_SC, yp_SC = M_SC[:,2], M_SC[:,3]
    t_SC = M_SC[:,4] * 3.335641e-9
    P_SC = M_SC[:,5]

    x_IBS, xp_IBS = M_IBS[:,0], M_IBS[:,1]
    y_IBS, yp_IBS = M_IBS[:,2], M_IBS[:,3]
    t_IBS = M_IBS[:,4] * 3.335641e-9 # ms (from mm/c to ms)
    P_IBS = M_IBS[:,5] # MeV

    E = np.sqrt(P**2 + mass**2)
    E_SC = np.sqrt(P_SC**2 + mass**2)
    E_IBS = np.sqrt(P_IBS**2 + mass**2)

    ES_arr[i] = np.std(E) # MeV
    ES_SC_arr[i] = np.std(E_SC)
    ES_SC_IBS_arr[i] = np.std(E_IBS)

    # ==========================
    # Cuts
    # ==========================
    def cut_mask(x, xp, y, yp, t, E):
        return (
            (np.abs(x - np.mean(x)) < cut*np.std(x)) &
            (np.abs(xp - np.mean(xp)) < cut*np.std(xp)) &
            (np.abs(y - np.mean(y)) < cut*np.std(y)) &
            (np.abs(yp - np.mean(yp)) < cut*np.std(yp)) &
            (np.abs(t - np.mean(t)) < cut*np.std(t)) &
            (np.abs(E - np.mean(E)) < cut*np.std(E))
        )

    mask = cut_mask(x,xp,y,yp,t,E)
    mask_SC = cut_mask(x_SC,xp_SC,y_SC,yp_SC,t_SC,E_SC)
    mask_IBS = cut_mask(x_IBS,xp_IBS,y_IBS,yp_IBS,t_IBS,E_IBS)

    X = np.column_stack([x[mask], xp[mask], y[mask], yp[mask]])
    X_SC = np.column_stack([x_SC[mask_SC], xp_SC[mask_SC], y_SC[mask_SC], yp_SC[mask_SC]])
    X_IBS = np.column_stack([x_IBS[mask_IBS], xp_IBS[mask_IBS], y_IBS[mask_IBS], yp_IBS[mask_IBS]])

    if len(X) < 10 or len(X_SC) < 10 or len(X_IBS) < 10:
        continue

    Z = np.column_stack([t[mask], E[mask]])
    Z_SC = np.column_stack([t_SC[mask_SC], E_SC[mask_SC]])
    Z_IBS = np.column_stack([t_IBS[mask_IBS], E_IBS[mask_IBS]])

    emittz[i] = np.sqrt(np.abs(np.linalg.det(np.cov(Z.T))))
    emittz_SC[i] = np.sqrt(np.abs(np.linalg.det(np.cov(Z_SC.T))))
    emittz_SC_IBS[i] = np.sqrt(np.abs(np.linalg.det(np.cov(Z_IBS.T))))

    emitt4d[i] = np.abs(np.linalg.det(np.cov(X.T)))**0.25 * (np.mean(P[mask]) / mass)
    emitt4d_SC[i] = np.abs(np.linalg.det(np.cov(X_SC.T)))**0.25 * (np.mean(P_SC[mask_SC]) / mass)
    emitt4d_SC_IBS[i] = np.abs(np.linalg.det(np.cov(X_IBS.T)))**0.25 * (np.mean(P_IBS[mask_IBS]) / mass)

# ==========================
# PLOTTING (no Octave needed)
# ==========================

plt.figure(figsize=(9,6))
plt.plot(length_arr, emitt4d, 'k')
plt.plot(length_SC_arr, emitt4d_SC, 'b')
plt.plot(length_SC_IBS_arr, emitt4d_SC_IBS, 'r')
plt.xlabel("z [m]", fontsize=18)
plt.ylabel("4D emittance [mm * mrad]", fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(["No SC", "SC", "SC + IBS"], fontsize=18)
plt.grid()
plt.savefig("emitt4d.png", dpi=200)
plt.savefig("emitt4d.svg")

plt.figure(figsize=(9,6))
plt.plot(length_arr, N_arr/N0*100, 'k')
plt.plot(length_SC_arr, N_SC_arr/N0*100, 'b')
plt.plot(length_SC_IBS_arr, N_SC_IBS_arr/N0*100, 'r')
plt.xlabel("z [m]", fontsize=18)
plt.ylabel("Transmission [%]", fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(["No SC", "SC", "SC + IBS"], fontsize=18)
plt.grid()
plt.savefig("transmission.png", dpi=200)
plt.savefig("transmission.svg")

plt.figure(figsize=(9,6))
plt.plot(length_arr, emittz * 1e6, 'k')
plt.plot(length_SC_arr, emittz_SC * 1e6, 'b')
plt.plot(length_SC_IBS_arr, emittz_SC_IBS * 1e6, 'r')
plt.xlabel("z [m]", fontsize=18)
plt.ylabel("Longitudinal emittance [eV * ms]", fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(["No SC", "SC", "SC + IBS"], fontsize=18)
plt.grid()
plt.savefig("emittz.png", dpi=200)
plt.savefig("emittz.svg")

plt.figure(figsize=(9,6))
plt.plot(length_arr, mean_K_arr, 'k')
plt.plot(length_SC_arr, mean_K_SC_arr, 'b')
plt.plot(length_SC_IBS_arr, mean_K_SC_IBS_arr, 'r')
plt.xlabel("z [m]", fontsize=18)
plt.ylabel("Kinetic energy [MeV]", fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(["No SC", "SC", "SC + IBS"], fontsize=18)
plt.grid()
plt.savefig("kinetic_energy.png", dpi=200)
plt.savefig("kinetic_energy.svg")


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

# -----------------------
# Top plot: bunch length
# -----------------------
ax1.plot(length_arr, sigma_t_arr, 'k')
ax1.plot(length_SC_arr, sigma_t_SC_arr, 'b')
ax1.plot(length_SC_IBS_arr, sigma_t_SC_IBS_arr, 'r')

ax1.set_ylabel("Bunch length [mm/c]", fontsize=18)
ax1.legend(["No SC", "SC", "SC + IBS"], fontsize=14)
ax1.grid()

ax1.tick_params(axis='both', labelsize=16)

# -----------------------
# Bottom plot: energy spread
# -----------------------
ax2.plot(length_arr, ES_arr, 'k')
ax2.plot(length_SC_arr, ES_SC_arr, 'b')
ax2.plot(length_SC_IBS_arr, ES_SC_IBS_arr, 'r')

ax2.set_xlabel("z [m]", fontsize=18)
ax2.set_ylabel("Energy spread [MeV]", fontsize=18)
ax2.legend(["No SC", "SC", "SC + IBS"], fontsize=14)
ax2.grid()

ax2.tick_params(axis='both', labelsize=16)

# -----------------------
# X-axis range (IMPORTANT)
# -----------------------
ax2.set_xlim(0, 17.5)

# Optional: make layout tight for paper
plt.tight_layout()

# Save high-quality outputs
plt.savefig("combined_plot.png", dpi=300, bbox_inches="tight")
plt.savefig("combined_plot.svg", bbox_inches="tight")

plt.show()
