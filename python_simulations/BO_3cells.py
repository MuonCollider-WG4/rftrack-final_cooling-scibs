"""Constrained multi-objective Bayesian optimization of three cooling cells.

RF-Track simulation has 24 optimized solenoid fields.
Two physics outputs are treated as equally important optimization objectives:

    1. minimize transverse 4D emittance,
    2. minimize longitudinal emittance.

Transmission - outcome constraint >= 0.986.

BoTorch acquisition --> maximize.
Two objectives: 1) objective_4d = -emitt4d / 164.4, 2) objective_z  = -emittz  / 1.2

Candidate selection uses constrained qLogNEHVI (log noisy expected hypervolume improvement), which learns and improves the feasible Pareto front.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
import torch
import RF_Track as rft

try:
    # Preferred numerically stable acquisition in current BoTorch.
    from botorch.acquisition.multi_objective.logei import (
        qLogNoisyExpectedHypervolumeImprovement as QNEHVI,
    )
    ACQUISITION_NAME = "qLogNEHVI"
except ImportError:  # compatibility with older BoTorch releases
    from botorch.acquisition.multi_objective.monte_carlo import (
        qNoisyExpectedHypervolumeImprovement as QNEHVI,
    )
    ACQUISITION_NAME = "qNEHVI"

from botorch.acquisition.multi_objective.objective import IdentityMCMultiOutputObjective
from botorch.fit import fit_gpytorch_mll
from botorch.models import ModelListGP, SingleTaskGP
from botorch.models.transforms import Normalize, Standardize
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.utils.multi_objective.pareto import is_non_dominated
from gpytorch.kernels import MaternKernel, ScaleKernel
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.mlls.sum_marginal_log_likelihood import SumMarginalLogLikelihood
from torch.quasirandom import SobolEngine

DTYPE = torch.double
DEVICE = torch.device("cpu")
torch.set_default_dtype(DTYPE)

# -----------------------------------------------------------------------------
# Problem definition
# -----------------------------------------------------------------------------
PARAMETER_NAMES = [
    "B_L_0", "B_MU_0", "B_MU_1", "B_MU_2",
    "B_MD_0", "B_MD_1", "B_MD_2", "B_MD_3", "B_MD_4", "B_MD_5",
    "B_MD_6", "B_MD_7", "B_MD_8", "B_MD_9", "B_MD_10", "B_MD_11",
    "B_MD_12", "B_MD_13", "B_MD_14", "B_MD_15", "B_MD_16", "B_MD_17",
    "B_MD_18", "B_MD_19",
]

START_RAW = np.array([
    4.272566009, 8.565005316, 0.8349035164, 12.98180977,
    7.31655761, 1.952982056, 0.5222768388, 2.638688899,
    -1.883843507, 1.99899282, -2.969758257, -0.6218776372,
    -2.451447979, -9.652429868, -7.226396811, -1.282173878,
    -2.558932364, -1.347886368, 1.919956457, -1.985248314,
    1.442434388, 2.827484677, 1.462318874, 8.498560988,
], dtype=float)

MAX_ABS_FIELD_T = 12.0
START = np.clip(START_RAW, -MAX_ABS_FIELD_T, MAX_ABS_FIELD_T)

# Fixed for this first optimization attempt
B_HUD_0 = 70.63117099
B_HUD_1 = 63.9043928
B_HUD_2 = 57.17761461
B_H = 40.36066914

# Targets stated in BO_3cells.py.
TARGET_EMITT4D = 164.4   # um
TARGET_EMITTZ = 1.2      # eV*ms
MIN_TRANSMISSION = 0.986

# RF-Track input.
DEFAULT_BEAM_FILE = Path(__file__).resolve().parent / "../050625BeamInput_10000_4.04124T.txt"

mass = rft.muonmass
Q = 1
charge = 5.3e12
SC = rft.SpaceCharge_PIC_FreeSpace(32, 32, 32)
rft.SC_engine = SC

@dataclass
class Metrics:
    emitt4d: float
    emittz: float
    transmission: float
    n_surviving_after_cut: int
    runtime_s: float
    failed: bool = False
    error: str = ""


def load_input_matrix(beam_file: Path) -> np.ndarray:
    beam_file = Path(beam_file).expanduser().resolve()
    if not beam_file.exists():
        raise FileNotFoundError(
            f"Beam input file not found: {beam_file}"
        )
    return np.loadtxt(beam_file)


def make_initial_bunch(input_matrix: np.ndarray):
    x = input_matrix[:, 0]
    y = input_matrix[:, 1]
    t = input_matrix[:, 6] * rft.ns
    Px = input_matrix[:, 3]
    Py = input_matrix[:, 4]
    Pz = input_matrix[:, 5]
    xp = Px / Pz * 1e3
    yp = Py / Pz * 1e3
    P = np.sqrt(Px**2 + Py**2 + Pz**2)
    B0 = rft.Bunch6d(mass, charge, Q, np.column_stack([x, xp, y, yp, t, P]))
    B0.set_lifetime(rft.muonlifetime)
    return B0


def cut_mask(x, xp, y, yp, t, E):
    return (
        (np.abs(x - np.mean(x)) < 3 * np.std(x))
        & (np.abs(xp - np.mean(xp)) < 3 * np.std(xp))
        & (np.abs(y - np.mean(y)) < 3 * np.std(y))
        & (np.abs(yp - np.mean(yp)) < 3 * np.std(yp))
        & (np.abs(t - np.mean(t)) < 3 * np.std(t))
        & (np.abs(E - np.mean(E)) < 3 * np.std(E))
    )


def evaluate_fields(fields_t: np.ndarray, input_matrix: np.ndarray) -> Metrics:
    """Run the expensive RF-Track simulation for one 24-D field vector."""
    t0 = time.perf_counter()
    fields_t = np.asarray(fields_t, dtype=float).reshape(-1)
    if fields_t.size != len(PARAMETER_NAMES):
        raise ValueError(f"Expected {len(PARAMETER_NAMES)} fields, got {fields_t.size}.")
    if np.any(np.abs(fields_t) > MAX_ABS_FIELD_T + 1e-12):
        raise ValueError("Candidate violates hard |B| <= 12 T constraint.")

    (
        B_L_0, B_MU_0, B_MU_1, B_MU_2,
        B_MD_0, B_MD_1, B_MD_2, B_MD_3, B_MD_4, B_MD_5,
        B_MD_6, B_MD_7, B_MD_8, B_MD_9, B_MD_10, B_MD_11,
        B_MD_12, B_MD_13, B_MD_14, B_MD_15, B_MD_16, B_MD_17,
        B_MD_18, B_MD_19,
    ) = fields_t.tolist()

    B0 = make_initial_bunch(input_matrix)
    L = rft.Lattice()
    V1 = rft.Volume()

    # ---- rft.Solenoid: L_0 ----
    L_0 = rft.Solenoid(1, B_L_0, 0.1, 0.3, 5)
    L_0.set_aperture(0.1)
    V1.add(L_0, 0.0, 0.0, 0, 'center')

    # ---- rft.Solenoid: MU_0 ----
    MU_0 = rft.Solenoid(0.2, B_MU_0, 0.1, 0.3, 5)
    MU_0.set_aperture(0.1)
    V1.add(MU_0, 0.0, 0.0, 1, 'center')

    # ---- rft.Solenoid: MU_1 ----
    MU_1 = rft.Solenoid(0.2, B_MU_1, 0.1, 0.3, 5)
    MU_1.set_aperture(0.1)
    V1.add(MU_1, 0.0, 0.0, 1.6, 'center')

    # ---- rft.Solenoid: MU_2 ----
    MU_2 = rft.Solenoid(0.2, B_MU_2, 0.1, 0.3, 5)
    MU_2.set_aperture(0.1)
    V1.add(MU_2, 0.0, 0.0, 1.9, 'center')

    # ---- rft.Solenoid: HU_0 ----
    HU_0 = rft.Solenoid(0.012, B_HUD_0, 0.03, 0.135, 5)
    HU_0.set_aperture(0.03)
    V1.add(HU_0, 0.0, 0.0, 2.106, 'center')

    # ---- rft.Solenoid: HU_1 ----
    HU_1 = rft.Solenoid(0.012, B_HUD_1, 0.03, 0.125, 5)
    HU_1.set_aperture(0.03)
    V1.add(HU_1, 0.0, 0.0, 2.118, 'center')

    # ---- rft.Solenoid: HU_2 ----
    HU_2 = rft.Solenoid(0.012, B_HUD_2, 0.03, 0.115, 5)
    HU_2.set_aperture(0.03)
    V1.add(HU_2, 0.0, 0.0, 2.13, 'center')

    # ---- rft.Solenoid: H_0 ----
    H_0 = rft.Solenoid(1.51736269, B_H, 0.03, 0.09, 5)
    H_0.set_aperture(0.03)
    V1.add(H_0, 0.0, 0.0, 2.894681345, 'center')

    # ---- rft.Absorber: abs_0 ----
    abs_0 = rft.Absorber(1.27736269, 'liquid_hydrogen')
    abs_0.set_aperture(0.025)
    V1.add(abs_0, 0.0, 0.0, 2.894681345, 'center')

    HD_2 = rft.Solenoid(0.012, B_HUD_2, 0.03, 0.115, 5)
    HD_2.set_aperture(0.03)
    V1.add(HD_2, 0.0, 0.0, 3.65936269, 'center')

    # ---- rft.Solenoid: HD_1 ----
    HD_1 = rft.Solenoid(0.012, B_HUD_1, 0.03, 0.125, 5)
    HD_1.set_aperture(0.03)
    V1.add(HD_1, 0.0, 0.0, 3.67136269, 'center')

    # ---- rft.Solenoid: HD_0 ----
    HD_0 = rft.Solenoid(0.012, B_HUD_0, 0.03, 0.135, 5)
    HD_0.set_aperture(0.03)
    V1.add(HD_0, 0.0, 0.0, 3.68336269, 'center')

    # ---- rft.Solenoid: MD_0 ----
    MD_0 = rft.Solenoid(0.2, B_MD_0, 0.1, 0.3, 5)
    MD_0.set_aperture(0.1)
    V1.add(MD_0, 0.0, 0.0, 3.88936269, 'center')

    rf1_1 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
    rf1_1.set_phid(-90.0)
    rf1_1.set_aperture(0.16)
    rf1_1.set_t0(1789.179061)
    V1.add(rf1_1, 0.0, 0.0, 4.16436269, 'center')

    # ---- rft.Solenoid: MD_1 ----
    MD_1 = rft.Solenoid(0.2, B_MD_1, 0.1, 0.3, 5)
    MD_1.set_aperture(0.1)
    V1.add(MD_1, 0.0, 0.0, 4.43936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_2 ----
    rf1_2 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
    rf1_2.set_phid(-90.0)
    rf1_2.set_aperture(0.16)
    rf1_2.set_t0(2454.045840)
    V1.add(rf1_2, 0.0, 0.0, 4.71436269, 'center')

    # ---- rft.Solenoid: MD_2 ----
    MD_2 = rft.Solenoid(0.2, B_MD_2, 0.1, 0.3, 5)
    MD_2.set_aperture(0.1)
    V1.add(MD_2, 0.0, 0.0, 4.98936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_3 ----
    rf1_3 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
    rf1_3.set_phid(-90.0)
    rf1_3.set_aperture(0.16)
    rf1_3.set_t0(120.955517)
    V1.add(rf1_3, 0.0, 0.0, 5.26436269, 'center')

    # ---- rft.Solenoid: MD_3 ----
    MD_3 = rft.Solenoid(0.2, B_MD_3, 0.1, 0.3, 5)
    MD_3.set_aperture(0.1)
    V1.add(MD_3, 0.0, 0.0, 5.53936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_4 ----
    rf1_4 = rft.Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1)
    rf1_4.set_phid(-90.0)
    rf1_4.set_aperture(0.16)
    rf1_4.set_t0(785.765173)
    V1.add(rf1_4, 0.0, 0.0, 5.81436269, 'center')

    # ---- rft.Solenoid: MD_4 ----
    MD_4 = rft.Solenoid(0.2, B_MD_4, 0.1, 0.3, 5)
    MD_4.set_aperture(0.1)
    V1.add(MD_4, 0.0, 0.0, 6.08936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_5 ----
    rf1_5 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
    rf1_5.set_phid(-10.0)
    rf1_5.set_aperture(0.16)
    rf1_5.set_t0(7445.973727)
    V1.add(rf1_5, 0.0, 0.0, 6.36436269, 'center')

    MD_5 = rft.Solenoid(0.2, B_MD_5, 0.1, 0.3, 5)
    MD_5.set_aperture(0.1)
    V1.add(MD_5, 0.0, 0.0, 6.63936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_6 ----
    rf1_6 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
    rf1_6.set_phid(-10.0)
    rf1_6.set_aperture(0.16)
    rf1_6.set_t0(8098.525725)
    V1.add(rf1_6, 0.0, 0.0, 6.91436269, 'center')

    # ---- rft.Solenoid: MD_6 ----
    MD_6 = rft.Solenoid(0.2, B_MD_6, 0.1, 0.3, 5)
    MD_6.set_aperture(0.1)
    V1.add(MD_6, 0.0, 0.0, 7.18936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_7 ----
    rf1_7 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
    rf1_7.set_phid(-10.0)
    rf1_7.set_aperture(0.16)
    rf1_7.set_t0(8742.313534)
    V1.add(rf1_7, 0.0, 0.0, 7.46436269, 'center')

    # ---- rft.Solenoid: MD_7 ----
    MD_7 = rft.Solenoid(0.2, B_MD_7, 0.1, 0.3, 5)
    MD_7.set_aperture(0.1)
    V1.add(MD_7, 0.0, 0.0, 7.73936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_8 ----
    rf1_8 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
    rf1_8.set_phid(-10.0)
    rf1_8.set_aperture(0.16)
    rf1_8.set_t0(9295.232106)
    V1.add(rf1_8, 0.0, 0.0, 8.01436269, 'center')

    MD_8 = rft.Solenoid(0.2, B_MD_8, 0.1, 0.3, 5)
    MD_8.set_aperture(0.1)
    V1.add(MD_8, 0.0, 0.0, 8.28936269, 'center')

    # ---- rft.Pillbox_Cavity: rf1_9 ----
    rf1_9 = rft.Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1)
    rf1_9.set_phid(-10.0)
    rf1_9.set_aperture(0.16)
    V1.add(rf1_9, 0.0, 0.0, 8.56436269, 'center')

    # ---- rft.Solenoid: MD_9 ----
    MD_9 = rft.Solenoid(0.2, B_MD_9, 0.1, 0.3, 5)
    MD_9.set_aperture(0.1)
    V1.add(MD_9, 0.0, 0.0, 8.83936269, 'center')

    # ---- rft.Solenoid: HU_3 ----
    HU_3 = rft.Solenoid(0.012, -B_HUD_0, 0.03, 0.135, 5)
    HU_3.set_aperture(0.03)
    V1.add(HU_3, 0.0, 0.0, 9.04536269, 'center')

    # ---- rft.Solenoid: HU_4 ----
    HU_4 = rft.Solenoid(0.012, -B_HUD_1, 0.03, 0.125, 5)
    HU_4.set_aperture(0.03)
    V1.add(HU_4, 0.0, 0.0, 9.05736269, 'center')

    # ---- rft.Solenoid: HU_5 ----
    HU_5 = rft.Solenoid(0.012, -B_HUD_2, 0.03, 0.115, 5)
    HU_5.set_aperture(0.03)
    V1.add(HU_5, 0.0, 0.0, 9.06936269, 'center')

    # ---- rft.Solenoid: H_1 ----
    H_1 = rft.Solenoid(1.56518469, -B_H, 0.03, 0.09, 5)
    H_1.set_aperture(0.03)
    V1.add(H_1, 0.0, 0.0, 9.857955035, 'center')

    # ---- rft.Absorber: rft.Abs_1 ----
    abs_1 = rft.Absorber(1.32518469, 'liquid_hydrogen')
    abs_1.set_aperture(0.025)
    V1.add(abs_1, 0.0, 0.0, 9.857955035, 'center')

    # ---- rft.Solenoid: HD_5 ----
    HD_5 = rft.Solenoid(0.012, -B_HUD_2, 0.03, 0.115, 5)
    HD_5.set_aperture(0.03)
    V1.add(HD_5, 0.0, 0.0, 10.64654738, 'center')

    # ---- rft.Solenoid: HD_4 ----
    HD_4 = rft.Solenoid(0.012, -B_HUD_1, 0.03, 0.125, 5)
    HD_4.set_aperture(0.03)
    V1.add(HD_4, 0.0, 0.0, 10.65854738, 'center')

    # ---- rft.Solenoid: HD_3 ----
    HD_3 = rft.Solenoid(0.012, -B_HUD_0, 0.03, 0.135, 5)
    HD_3.set_aperture(0.03)
    V1.add(HD_3, 0.0, 0.0, 10.67054738, 'center')

    # ---- rft.Solenoid: MD_10 ----
    MD_10 = rft.Solenoid(0.2, B_MD_10, 0.1, 0.3, 5)
    MD_10.set_aperture(0.1)
    V1.add(MD_10, 0.0, 0.0, 10.87654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_1 ----
    rf2_1 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
    rf2_1.set_phid(-90.0)
    rf2_1.set_aperture(0.16)
    V1.add(rf2_1, 0.0, 0.0, 11.15154738, 'center')

    # ---- rft.Solenoid: MD_11 ----
    MD_11 = rft.Solenoid(0.2, B_MD_11, 0.1, 0.3, 5)
    MD_11.set_aperture(0.1)
    V1.add(MD_11, 0.0, 0.0, 11.42654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_2 ----
    rf2_2 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
    rf2_2.set_phid(-90.0)
    rf2_2.set_aperture(0.16)
    V1.add(rf2_2, 0.0, 0.0, 11.70154738, 'center')

    # ---- rft.Solenoid: MD_12 ----
    MD_12 = rft.Solenoid(0.2, B_MD_12, 0.1, 0.3, 5)
    MD_12.set_aperture(0.1)
    V1.add(MD_12, 0.0, 0.0, 11.97654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_3 ----
    rf2_3 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
    rf2_3.set_phid(-90.0)
    rf2_3.set_aperture(0.16)
    V1.add(rf2_3, 0.0, 0.0, 12.25154738, 'center')

    # ---- rft.Solenoid: MD_13 ----
    MD_13 = rft.Solenoid(0.2, B_MD_13, 0.1, 0.3, 5)
    MD_13.set_aperture(0.1)
    V1.add(MD_13, 0.0, 0.0, 12.52654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_4 ----
    rf2_4 = rft.Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1)
    rf2_4.set_phid(-90.0)
    rf2_4.set_aperture(0.16)
    V1.add(rf2_4, 0.0, 0.0, 12.80154738, 'center')

    # ---- rft.Solenoid: MD_14 ----
    MD_14 = rft.Solenoid(0.2, B_MD_14, 0.1, 0.3, 5)
    MD_14.set_aperture(0.1)
    V1.add(MD_14, 0.0, 0.0, 13.07654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_5 ----
    rf2_5 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
    rf2_5.set_phid(-10.0)
    rf2_5.set_aperture(0.16)
    V1.add(rf2_5, 0.0, 0.0, 13.35154738, 'center')

    # ---- rft.Solenoid: MD_15 ----
    MD_15 = rft.Solenoid(0.2, B_MD_15, 0.1, 0.3, 5)
    MD_15.set_aperture(0.1)
    V1.add(MD_15, 0.0, 0.0, 13.62654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_6 ----
    rf2_6 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
    rf2_6.set_phid(-10.0)
    rf2_6.set_aperture(0.16)
    V1.add(rf2_6, 0.0, 0.0, 13.90154738, 'center')

    # ---- rft.Solenoid: MD_16 ----
    MD_16 = rft.Solenoid(0.2, B_MD_16, 0.1, 0.3, 5)
    MD_16.set_aperture(0.1)
    V1.add(MD_16, 0.0, 0.0, 14.17654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_7 ----
    rf2_7 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
    rf2_7.set_phid(-10.0)
    rf2_7.set_aperture(0.16)
    V1.add(rf2_7, 0.0, 0.0, 14.45154738, 'center')

    # ---- rft.Solenoid: MD_17 ----
    MD_17 = rft.Solenoid(0.2, B_MD_17, 0.1, 0.3, 5)
    MD_17.set_aperture(0.1)
    V1.add(MD_17, 0.0, 0.0, 14.72654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_8 ----
    rf2_8 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
    rf2_8.set_phid(-10.0)
    rf2_8.set_aperture(0.16)
    V1.add(rf2_8, 0.0, 0.0, 15.00154738, 'center')

    # ---- rft.Solenoid: MD_18 ----
    MD_18 = rft.Solenoid(0.2, B_MD_18, 0.1, 0.3, 5)
    MD_18.set_aperture(0.1)
    V1.add(MD_18, 0.0, 0.0, 15.27654738, 'center')

    # ---- rft.Pillbox_Cavity: rf2_9 ----
    rf2_9 = rft.Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1)
    rf2_9.set_phid(-10.0)
    rf2_9.set_aperture(0.16)
    V1.add(rf2_9, 0.0, 0.0, 15.55154738, 'center')

    # ---- rft.Solenoid: MD_19 ----
    MD_19 = rft.Solenoid(0.2, B_MD_19, 0.1, 0.3, 5)
    MD_19.set_aperture(0.1)
    V1.add(MD_19, 0.0, 0.0, 15.82654738, 'center')

    # ---- rft.Solenoid: HU_6 ----
    HU_6 = rft.Solenoid(0.012, B_HUD_0, 0.03, 0.135, 5)
    HU_6.set_aperture(0.03)
    V1.add(HU_6, 0.0, 0.0, 16.03254738, 'center')

    # ---- rft.Solenoid: HU_7 ----
    HU_7 = rft.Solenoid(0.012, B_HUD_1, 0.03, 0.125, 5)
    HU_7.set_aperture(0.03)
    V1.add(HU_7, 0.0, 0.0, 16.04454738, 'center')

    # ---- rft.Solenoid: HU_8 ----
    HU_8 = rft.Solenoid(0.012, B_HUD_2, 0.03, 0.115, 5)
    HU_8.set_aperture(0.03)
    V1.add(HU_8, 0.0, 0.0, 16.05654738, 'center')

    # ---- rft.rft.Solenoid: H_2 ----
    H_2 = rft.Solenoid(1.34521782, B_H, 0.03, 0.09, 5)
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

    M0 = B0.get_phase_space()
    B_refp = rft.Bunch6d(mass, charge, Q, np.mean(M0, axis=0))
    B1_refp = V1.autophase(B_refp)
    V1.set_tt_nsteps(250)
    L.append(V1)
    B1 = L.track(B0)

    M = B1.get_phase_space()
    x, xp = M[:, 0], M[:, 1]
    y, yp = M[:, 2], M[:, 3]
    t = M[:, 4] * 3.335641e-9
    P = M[:, 5]
    E = np.sqrt(P**2 + mass**2)
    mask = cut_mask(x, xp, y, yp, t, E)

    X_phase = np.column_stack([x[mask], xp[mask], y[mask], yp[mask]])
    Z_phase = np.column_stack([t[mask], E[mask]])

    transmission = float(B1.get_info().transmission / charge)
    emitt4d = float(
        np.abs(np.linalg.det(np.cov(X_phase.T))) ** 0.25
        * (np.mean(P[mask]) / mass)
    )
    emittz = float(np.sqrt(np.abs(np.linalg.det(np.cov(Z_phase.T)))) * 1e6)

    if not np.all(np.isfinite([emitt4d, emittz, transmission])):
        raise RuntimeError("Non-finite objective metric returned by RF-Track.")

    return Metrics(
        emitt4d=emitt4d,
        emittz=emittz,
        transmission=transmission,
        runtime_s=time.perf_counter() - t0,
    )


def safe_evaluate(fields_t: np.ndarray, input_matrix: np.ndarray) -> Metrics:
    """Convert simulation failures into a finite bad BO observation."""
    try:
        return evaluate_fields(fields_t, input_matrix)
    except Exception as exc:
        # Keep failures so they can safely enter GP training.
        return Metrics(
            emitt4d=10.0 * TARGET_EMITT4D,
            emittz=10.0 * TARGET_EMITTZ,
            transmission=0.0,
            runtime_s=0.0,
            failed=True,
            error=f"{type(exc).__name__}: {exc}",
        )


# -----------------------------------------------------------------------------
# Bayesian optimization helpers
# -----------------------------------------------------------------------------
def make_search_bounds(local_half_width_t: float) -> torch.Tensor:
    """Hard [-12, 12] T bounds intersected with a local box around START."""
    lower = np.maximum(-MAX_ABS_FIELD_T, START - local_half_width_t)
    upper = np.minimum(MAX_ABS_FIELD_T, START + local_half_width_t)
    return torch.tensor(np.vstack([lower, upper]), dtype=DTYPE, device=DEVICE)


def initial_design(
    bounds: torch.Tensor,
    n_initial: int,
    sobol_radius_t: float,
    seed: int,
) -> torch.Tensor:
    """Start exactly at supplied settings, then add local Sobol perturbations."""
    if n_initial < 1:
        raise ValueError("n_initial must be >= 1")

    center = torch.tensor(START, dtype=DTYPE, device=DEVICE)
    if n_initial == 1:
        return center.unsqueeze(0)

    global_lo, global_hi = bounds
    lo = torch.maximum(global_lo, center - sobol_radius_t)
    hi = torch.minimum(global_hi, center + sobol_radius_t)

    sobol = SobolEngine(dimension=len(PARAMETER_NAMES), scramble=True, seed=seed)
    u = sobol.draw(n_initial - 1).to(dtype=DTYPE, device=DEVICE)
    perturb = lo + (hi - lo) * u
    return torch.cat([center.unsqueeze(0), perturb], dim=0)


def metrics_to_model_outputs(m: Metrics) -> torch.Tensor:
    """Convert physical metrics to the three GP outputs.

    Outputs 0 and 1 are maximized by BoTorch. Negating the normalized
    emittances therefore corresponds to minimizing the physical emittances.

    Output 2 is a constraint value. BoTorch treats constraint <= 0 as
    feasible, so MIN_TRANSMISSION - transmission <= 0 means feasible.
    """
    return torch.tensor(
        [
            -m.emitt4d / TARGET_EMITT4D,
            -m.emittz / TARGET_EMITTZ,
            MIN_TRANSMISSION - m.transmission,
        ],
        dtype=DTYPE,
        device=DEVICE,
    )


def build_model(
    train_X: torch.Tensor,
    train_Y: torch.Tensor,
    bounds: torch.Tensor,
) -> ModelListGP:
    """Fit three independent ARD Matern-5/2 GPs.

    GP 0: normalized negative transverse emittance
    GP 1: normalized negative longitudinal emittance
    GP 2: transmission constraint value
    """
    d = train_X.shape[-1]
    models = []
    for output_idx in range(3):
        covar_module = ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=d))
        gp = SingleTaskGP(
            train_X=train_X,
            train_Y=train_Y[:, output_idx : output_idx + 1],
            input_transform=Normalize(d=d, bounds=bounds),
            outcome_transform=Standardize(m=1),
            covar_module=covar_module,
        )
        models.append(gp)

    model = ModelListGP(*models)
    mll = SumMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    return model


def feasible_mask(Y: torch.Tensor) -> torch.Tensor:
    """Observed feasibility: constraint output <= 0."""
    return Y[:, 2] <= 0.0


def pareto_indices(Y: torch.Tensor) -> torch.Tensor:
    """Indices of the observed feasible non-dominated solutions."""
    feasible = feasible_mask(Y)
    feasible_idx = torch.where(feasible)[0]
    if feasible_idx.numel() == 0:
        return torch.empty(0, dtype=torch.long, device=Y.device)
    nd_local = is_non_dominated(Y[feasible, :2])
    return feasible_idx[nd_local]


def choose_trust_center(X: torch.Tensor, Y: torch.Tensor) -> torch.Tensor:
    """Choose a smooth-search center without scalarizing the BO objectives.

    If feasible Pareto points exist, choose the Pareto point closest to the
    ideal objective corner (0, 0). Since objectives are target-normalized,
    Euclidean distance gives the two emittances equal scale for this *search
    center only*. The acquisition itself remains genuinely multi-objective.

    If no feasible point exists yet, center on the point with the smallest
    transmission-constraint violation.
    """
    pidx = pareto_indices(Y)
    if pidx.numel() > 0:
        pareto_obj = Y[pidx, :2]
        distance_to_ideal = torch.linalg.vector_norm(pareto_obj, dim=1)
        return X[pidx[torch.argmin(distance_to_ideal)]]

    least_infeasible_idx = torch.argmin(Y[:, 2])
    return X[least_infeasible_idx]


def trust_region_bounds(
    global_bounds: torch.Tensor,
    center_x: torch.Tensor,
    radius_t: float,
) -> torch.Tensor:
    """Restrict each next proposal to a smooth change around a chosen center."""
    lo = torch.maximum(global_bounds[0], center_x - radius_t)
    hi = torch.minimum(global_bounds[1], center_x + radius_t)
    tiny = torch.tensor(1e-9, dtype=DTYPE, device=DEVICE)
    hi = torch.maximum(hi, lo + tiny)
    return torch.stack([lo, hi])


def make_reference_point(Y: torch.Tensor, ref_multiplier: float) -> list[float]:
    """Reference point for hypervolume in target-normalized objective space.

    A value of 2.0 means that, by default, twice each target emittance is the
    hypervolume reference. If observed feasible points are worse than that,
    move the reference slightly lower so the current feasible Pareto front is
    still represented.
    """
    if ref_multiplier <= 1.0:
        raise ValueError("ref_multiplier must be > 1.0")

    ref = torch.full((2,), -float(ref_multiplier), dtype=DTYPE, device=DEVICE)
    feasible = feasible_mask(Y)
    if feasible.any():
        observed_min = Y[feasible, :2].min(dim=0).values
        ref = torch.minimum(ref, observed_min - 0.05)
    return ref.tolist()


def propose_qnehvi(
    model: ModelListGP,
    train_X: torch.Tensor,
    train_Y: torch.Tensor,
    global_bounds: torch.Tensor,
    trust_center: torch.Tensor,
    trust_radius_t: float,
    ref_multiplier: float,
    mc_samples: int,
    num_restarts: int,
    raw_samples: int,
) -> torch.Tensor:
    """Propose one point using constrained hypervolume improvement."""
    sampler = SobolQMCNormalSampler(sample_shape=torch.Size([mc_samples]))
    objective = IdentityMCMultiOutputObjective(outcomes=[0, 1])

    # Constraint convention: negative (or zero) means feasible.
    constraints = [lambda samples: samples[..., 2]]
    ref_point = make_reference_point(train_Y, ref_multiplier)

    acq = QNEHVI(
        model=model,
        ref_point=ref_point,
        X_baseline=train_X,
        sampler=sampler,
        prune_baseline=True,
        objective=objective,
        constraints=constraints,
    )

    local_bounds = trust_region_bounds(
        global_bounds=global_bounds,
        center_x=trust_center,
        radius_t=trust_radius_t,
    )
    candidate, _ = optimize_acqf(
        acq_function=acq,
        bounds=local_bounds,
        q=1,
        num_restarts=num_restarts,
        raw_samples=raw_samples,
        options={"batch_limit": 5, "maxiter": 250},
        sequential=True,
    )
    return candidate.detach()
def append_csv(path: Path, iteration: int, phase: str, x: np.ndarray, m: Metrics) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    obj4d = -m.emitt4d / TARGET_EMITT4D
    objz = -m.emittz / TARGET_EMITTZ
    constraint = MIN_TRANSMISSION - m.transmission
    with path.open("a", newline="") as f:
        writer = csv.writer(f)
        if not exists:
            writer.writerow(
                [
                    "iteration", "phase", *PARAMETER_NAMES,
                    "emitt4d_um", "emittz_eV_ms", "transmission",
                    "objective_4d", "objective_z", "transmission_constraint",
                    "feasible", "runtime_s", "failed", "error",
                ]
            )
        writer.writerow(
            [
                iteration, phase, *map(float, x),
                m.emitt4d, m.emittz, m.transmission,
                obj4d, objz, constraint, constraint <= 0.0,
                m.runtime_s, m.failed, m.error,
            ]
        )


def save_checkpoint(path: Path, X: torch.Tensor, Y: torch.Tensor, metrics: list[Metrics]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "X": X.cpu(),
            "Y": Y.cpu(),
            "Y_columns": [
                "-emitt4d/TARGET_EMITT4D",
                "-emittz/TARGET_EMITTZ",
                "MIN_TRANSMISSION-transmission",
            ],
            "metrics": [m.__dict__ for m in metrics],
            "parameter_names": PARAMETER_NAMES,
            "start_raw": START_RAW,
            "start_feasible": START,
            "targets": {
                "emitt4d_um": TARGET_EMITT4D,
                "emittz_eV_ms": TARGET_EMITTZ,
                "min_transmission": MIN_TRANSMISSION,
            },
        },
        path,
    )


def print_result(tag: str, i: int, m: Metrics) -> None:
    status = "FAILED" if m.failed else "ok"
    feasible = (not m.failed) and (m.transmission >= MIN_TRANSMISSION)
    print(
        f"[{tag} {i:03d}] {status:6s} feasible={str(feasible):5s}  "
        f"e4d={m.emitt4d:.6g} um  ez={m.emittz:.6g} eV*ms  "
        f"trans={m.transmission:.6g}  time={m.runtime_s:.1f}s"
    )
    if m.failed:
        print(f"    {m.error}")


def print_pareto_summary(Y: torch.Tensor, metrics: list[Metrics]) -> None:
    pidx = pareto_indices(Y)
    if pidx.numel() == 0:
        best_constraint_idx = int(torch.argmin(Y[:, 2]).item())
        m = metrics[best_constraint_idx]
        print(
            "    no feasible point yet; best transmission so far: "
            f"{m.transmission:.6g} (required >= {MIN_TRANSMISSION})"
        )
        return

    # A single representative is useful for console monitoring only. The full
    # Pareto set remains the actual optimization result.
    pareto_obj = Y[pidx, :2]
    balanced_local = torch.argmin(torch.linalg.vector_norm(pareto_obj, dim=1))
    balanced_idx = int(pidx[balanced_local].item())
    m = metrics[balanced_idx]
    print(
        f"    feasible Pareto points={pidx.numel()}; "
        f"balanced representative: e4d={m.emitt4d:.6g}, "
        f"ez={m.emittz:.6g}, trans={m.transmission:.6g}"
    )


def run_bo(args: argparse.Namespace) -> Tuple[torch.Tensor, torch.Tensor, list[Metrics]]:
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    input_matrix = load_input_matrix(Path(args.beam_file))
    global_bounds = make_search_bounds(args.local_half_width)

    if np.any(START_RAW != START):
        changed = [
            f"{name}: {raw:.6g} -> {feasible:.6g} T"
            for name, raw, feasible in zip(PARAMETER_NAMES, START_RAW, START)
            if raw != feasible
        ]
        print("Hard |B| <= 12 T constraint clips the supplied start: " + ", ".join(changed))

    print(f"BO dimension: {len(PARAMETER_NAMES)}")
    print("Objectives: minimize emitt4d AND minimize emittz (equal target-normalized scale)")
    print(f"Outcome constraint: transmission >= {MIN_TRANSMISSION}")
    print(f"Acquisition: {ACQUISITION_NAME}")
    print(f"Initial evaluations: {args.n_initial}; BO iterations: {args.n_iter}")
    print(f"Trust radius={args.trust_radius} T; search half-width={args.local_half_width} T")
    print(f"Hypervolume reference multiplier={args.ref_multiplier}; MC samples={args.mc_samples}")

    csv_path = Path(args.output_csv)
    ckpt_path = Path(args.checkpoint)
    if args.fresh:
        csv_path.unlink(missing_ok=True)
        ckpt_path.unlink(missing_ok=True)

    X = initial_design(global_bounds, args.n_initial, args.sobol_radius, args.seed)
    Y_rows = []
    all_metrics: list[Metrics] = []

    for i in range(X.shape[0]):
        x_np = X[i].cpu().numpy()
        m = safe_evaluate(x_np, input_matrix)
        Y_rows.append(metrics_to_model_outputs(m))
        all_metrics.append(m)
        append_csv(csv_path, i, "initial", x_np, m)
        print_result("init", i, m)

    Y = torch.stack(Y_rows, dim=0)
    save_checkpoint(ckpt_path, X, Y, all_metrics)
    print_pareto_summary(Y, all_metrics)

    for it in range(args.n_iter):
        model = build_model(X, Y, global_bounds)
        center_x = choose_trust_center(X, Y)

        next_X = propose_qnehvi(
            model=model,
            train_X=X,
            train_Y=Y,
            global_bounds=global_bounds,
            trust_center=center_x,
            trust_radius_t=args.trust_radius,
            ref_multiplier=args.ref_multiplier,
            mc_samples=args.mc_samples,
            num_restarts=args.num_restarts,
            raw_samples=args.raw_samples,
        )
        next_np = next_X.squeeze(0).cpu().numpy()
        m = safe_evaluate(next_np, input_matrix)
        next_Y = metrics_to_model_outputs(m).unsqueeze(0)

        X = torch.cat([X, next_X], dim=0)
        Y = torch.cat([Y, next_Y], dim=0)
        all_metrics.append(m)
        append_csv(csv_path, it, "bo", next_np, m)
        save_checkpoint(ckpt_path, X, Y, all_metrics)
        print_result("bo", it, m)
        print_pareto_summary(Y, all_metrics)

    return X, Y, all_metrics


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--beam-file", type=str, default=str(DEFAULT_BEAM_FILE))
    p.add_argument(
        "--n-initial", type=int, default=24,
        help="Total initial evaluations including the supplied start.",
    )
    p.add_argument(
        "--n-iter", type=int, default=60,
        help="Sequential multi-objective BO evaluations after initialization.",
    )
    p.add_argument(
        "--local-half-width", type=float, default=4.0,
        help="Optimization domain is START +/- this many tesla, clipped to +/-12 T.",
    )
    p.add_argument(
        "--sobol-radius", type=float, default=1.0,
        help="Initial Sobol samples use START +/- this many tesla.",
    )
    p.add_argument(
        "--trust-radius", type=float, default=1.5,
        help="Each proposed candidate stays within this many tesla of the trust center.",
    )
    p.add_argument(
        "--ref-multiplier", type=float, default=2.0,
        help="Hypervolume reference corresponds initially to this multiple of each emittance target.",
    )
    p.add_argument(
        "--mc-samples", type=int, default=128,
        help="Sobol QMC posterior samples used by qLogNEHVI/qNEHVI.",
    )
    p.add_argument("--num-restarts", type=int, default=20)
    p.add_argument("--raw-samples", type=int, default=1024)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--output-csv", type=str, default="bo_3cells_mobo_history.csv")
    p.add_argument("--checkpoint", type=str, default="bo_3cells_mobo_checkpoint.pt")
    p.add_argument(
        "--fresh", action="store_true",
        help="Delete previous CSV/checkpoint before starting.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    X, Y, metrics = run_bo(args)

    pidx = pareto_indices(Y)
    print("\n=== FINAL FEASIBLE PARETO FRONT ===")
    if pidx.numel() == 0:
        print("No feasible configuration satisfying the transmission constraint was observed.")
        return

    for rank, idx_t in enumerate(pidx, start=1):
        idx = int(idx_t.item())
        m = metrics[idx]
        print(
            f"\nPareto {rank}: evaluation={idx}, "
            f"emitt4d={m.emitt4d:.8g} um, "
            f"emittz={m.emittz:.8g} eV*ms, "
            f"transmission={m.transmission:.8g}"
        )
        for name, value in zip(PARAMETER_NAMES, X[idx].cpu().numpy()):
            print(f"  {name:8s} = {value:+.9f} T")

    pareto_obj = Y[pidx, :2]
    balanced_local = torch.argmin(torch.linalg.vector_norm(pareto_obj, dim=1))
    balanced_idx = int(pidx[balanced_local].item())
    balanced_m = metrics[balanced_idx]
    print("\n=== EQUAL-SCALE COMPROMISE FROM THE PARETO FRONT ===")
    print(
        "This is only a convenient representative; the Pareto front above is "
        "the actual multi-objective result."
    )
    print(json.dumps(balanced_m.__dict__, indent=2))


if __name__ == "__main__":
    main()
