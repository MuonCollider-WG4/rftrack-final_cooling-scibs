%% SHORT SIMULATION %%%
%% Date: 06/03

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Bunch creation
%%%%%%%%%%%%%%%%%%%%%%%
RF_Track
mass = RF_Track.muonmass;
Q = 1; % to be checked
charge = 5.3e12; % It is 5.3e12 muons
M = load('../050625BeamInput_10000_4.04124T.txt');
x = M(:, 1); % x in mm
y = M(:, 2); % y in mm
t = M(:, 7) * RF_Track.ns; % z in ns
Px = M(:, 4); % Px in MeV/c
Py = M(:, 5); % Py in MeV/c
Pz = M(:, 6); % Pz in MeV/c
xp = Px ./ Pz * 1e3; % mrad
yp = Py ./ Pz * 1e3; % mrad
P = hypot(Px, Py, Pz);
B0 = Bunch6d(mass, charge, Q, [x xp y yp t P]);
B0.set_lifetime (RF_Track.muonlifetime);
L = Lattice();
V1 = Volume(); 

NX = 32; % Mesh cells in x
NY = 32; % Mesh cells in y
NZ = 32; % Mesh cells in z
Ncol = 200; % Maximum number of collisions per kick qith IBS (if more needed, uses a formula)
SC = SpaceCharge_PIC_FreeSpace(NX, NY, NZ); % Space charge definition
SC.set_smooth(0.5); % SC - set smooth factor
SC.set_mirror(0.0); % set position of cathode
RF_Track.SC_engine = SC; % SC configuration in RF-Track


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Function definition
%%%%%%%%%%%%%%%%%%%%%%%
function peak_field = get_peak_field(cur, r_out, r_in, L)
    % cur: current density in A/m^2
    % r_out, r_in: outer and inner radius in m
    % L: solenoid length in m
    
    mu0 = 4*pi*1e-7;

    % Use center of solenoid approximation (on-axis B)
    term_out = r_out;
    term_in  = r_in;
    
    peak_field = mu0 * cur * (term_out - term_in);  
end


% ---- Solenoid: L_0 ----
L_0 = Solenoid(1, get_peak_field(17000000.0, 0.3, 0.1, 1), 0.1, 0.3, 5);
L_0.set_aperture(0.1);
V1.add(L_0, 0.0, 0.0, 0, reference='center');
% ---- Solenoid: MU_0 ----
MU_0 = Solenoid(0.2, get_peak_field(34079073.34, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MU_0.set_aperture(0.1);
V1.add(MU_0, 0.0, 0.0, 1, reference='center');
% ---- Solenoid: MU_1 ----
MU_1 = Solenoid(0.2, get_peak_field(3321975.541, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MU_1.set_aperture(0.1);
V1.add(MU_1, 0.0, 0.0, 1.6, reference='center');
% ---- Solenoid: MU_2 ----
MU_2 = Solenoid(0.2, get_peak_field(51652979.870000005, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MU_2.set_aperture(0.1);
V1.add(MU_2, 0.0, 0.0, 1.9, reference='center');
% ---- Solenoid: HU_0 ----
HU_0 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5);
HU_0.set_aperture(0.03);
V1.add(HU_0, 0.0, 0.0, 2.106, reference='center');
% ---- Solenoid: HU_1 ----
HU_1 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5);
HU_1.set_aperture(0.03);
V1.add(HU_1, 0.0, 0.0, 2.118, reference='center');
% ---- Solenoid: HU_2 ----
HU_2 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5);
HU_2.set_aperture(0.03);
V1.add(HU_2, 0.0, 0.0, 2.13, reference='center');
% ---- Solenoid: H_0 ----
H_0 = Solenoid(1.51736269, get_peak_field(535299999.99999994, 0.09, 0.03, 1.51736269), 0.03, 0.09, 5);
H_0.set_aperture(0.03);
V1.add(H_0, 0.0, 0.0, 2.894681345, reference='center');
% ---- Absorber: abs_0 ----
abs_0 = Absorber(1.27736269, 'liquid_hydrogen');
abs_0.set_aperture(0.025);
V1.add(abs_0, 0.0, 0.0, 2.894681345, reference='center');
%V1.add(Drift(1.27736269), 0.0, 0.0, 2.894681345, reference='center');

HD_2 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5);
HD_2.set_aperture(0.03);
V1.add(HD_2, 0.0, 0.0, 3.65936269, reference='center');

% ---- Solenoid: HD_1 ----
HD_1 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5);
HD_1.set_aperture(0.03);
V1.add(HD_1, 0.0, 0.0, 3.67136269, reference='center');
% ---- Solenoid: HD_0 ----
HD_0 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5);
HD_0.set_aperture(0.03);
V1.add(HD_0, 0.0, 0.0, 3.68336269, reference='center');
% ---- Solenoid: MD_0 ----
MD_0 = Solenoid(0.2, get_peak_field(29111657.75, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_0.set_aperture(0.1);
V1.add(MD_0, 0.0, 0.0, 3.88936269, reference='center');

rf1_1 = Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1);
rf1_1.set_phid(-90.0);
rf1_1.set_aperture(0.16);
rf1_1.set_t0(1789.179061);
V1.add(rf1_1, 0.0, 0.0, 4.16436269, reference='center');

% ---- Solenoid: MD_1 ----
MD_1 = Solenoid(0.2, get_peak_field(7770668.698, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_1.set_aperture(0.1);
V1.add(MD_1, 0.0, 0.0, 4.43936269, reference='center');
% ---- Pillbox_Cavity: rf1_2 ----
rf1_2 = Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1);
rf1_2.set_phid(-90.0);
rf1_2.set_aperture(0.16);
rf1_2.set_t0(2454.045840);
V1.add(rf1_2, 0.0, 0.0, 4.71436269, reference='center');

% ---- Solenoid: MD_2 ----
MD_2 = Solenoid(0.2, get_peak_field(2078073.5140000002, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_2.set_aperture(0.1);
V1.add(MD_2, 0.0, 0.0, 4.98936269, reference='center');
% ---- Pillbox_Cavity: rf1_3 ----
rf1_3 = Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1);
rf1_3.set_phid(-90.0);
rf1_3.set_aperture(0.16);
rf1_3.set_t0(120.955517);
V1.add(rf1_3, 0.0, 0.0, 5.26436269, reference='center');

% ---- Solenoid: MD_3 ----
MD_3 = Solenoid(0.2, get_peak_field(10499009.54, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_3.set_aperture(0.1);
V1.add(MD_3, 0.0, 0.0, 5.53936269, reference='center');
% ---- Pillbox_Cavity: rf1_4 ----
rf1_4 = Pillbox_Cavity(35000000.0, 100000000.0, 0.25, 1);
rf1_4.set_phid(-90.0);
rf1_4.set_aperture(0.16);
rf1_4.set_t0(785.765173);
V1.add(rf1_4, 0.0, 0.0, 5.81436269, reference='center');

% ---- Solenoid: MD_4 ----
MD_4 = Solenoid(0.2, get_peak_field(-7495575.152, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_4.set_aperture(0.1);
V1.add(MD_4, 0.0, 0.0, 6.08936269, reference='center');
% ---- Pillbox_Cavity: rf1_5 ----
rf1_5 = Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1);
rf1_5.set_phid(-10.0);
rf1_5.set_aperture(0.16);
rf1_5.set_t0(7445.973727);
V1.add(rf1_5, 0.0, 0.0, 6.36436269, reference='center');

MD_5 = Solenoid(0.2, get_peak_field(7953739.713, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_5.set_aperture(0.1);
V1.add(MD_5, 0.0, 0.0, 6.63936269, reference='center');
% ---- Pillbox_Cavity: rf1_6 ----
rf1_6 = Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1);
rf1_6.set_phid(-10.0);
rf1_6.set_aperture(0.16);
rf1_6.set_t0(8098.525725);
V1.add(rf1_6, 0.0, 0.0, 6.91436269, reference='center');

% ---- Solenoid: MD_6 ----
MD_6 = Solenoid(0.2, get_peak_field(-11816292.66, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_6.set_aperture(0.1);
V1.add(MD_6, 0.0, 0.0, 7.18936269, reference='center');
% ---- Pillbox_Cavity: rf1_7 ----
rf1_7 = Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1);
rf1_7.set_phid(-10.0);
rf1_7.set_aperture(0.16);
rf1_7.set_t0(8742.313534);
V1.add(rf1_7, 0.0, 0.0, 7.46436269, reference='center');

% ---- Solenoid: MD_7 ----
MD_7 = Solenoid(0.2, get_peak_field(-2474372.499, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_7.set_aperture(0.1);
V1.add(MD_7, 0.0, 0.0, 7.73936269, reference='center');
% ---- Pillbox_Cavity: rf1_8 ----
rf1_8 = Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1);
rf1_8.set_phid(-10.0);
rf1_8.set_aperture(0.16);
rf1_8.set_t0(9295.232106);
V1.add(rf1_8, 0.0, 0.0, 8.01436269, reference='center');

MD_8 = Solenoid(0.2, get_peak_field(-9754001.591, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_8.set_aperture(0.1);
V1.add(MD_8, 0.0, 0.0, 8.28936269, reference='center');
% ---- Pillbox_Cavity: rf1_9 ----
rf1_9 = Pillbox_Cavity(33000000.0, 1000000.0, 0.25, 1);
rf1_9.set_phid(-10.0);
rf1_9.set_aperture(0.16);
V1.add(rf1_9, 0.0, 0.0, 8.56436269, reference='center');

% ---- Solenoid: MD_9 ----
MD_9 = Solenoid(0.2, get_peak_field(-38405798.160000004, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_9.set_aperture(0.1);
V1.add(MD_9, 0.0, 0.0, 8.83936269, reference='center');
% ---- Solenoid: HU_3 ----
HU_3 = Solenoid(0.012, get_peak_field(-535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5);
HU_3.set_aperture(0.03);
V1.add(HU_3, 0.0, 0.0, 9.04536269, reference='center');
% ---- Solenoid: HU_4 ----
HU_4 = Solenoid(0.012, get_peak_field(-535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5);
HU_4.set_aperture(0.03);
V1.add(HU_4, 0.0, 0.0, 9.05736269, reference='center');
% ---- Solenoid: HU_5 ----
HU_5 = Solenoid(0.012, get_peak_field(-535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5);
HU_5.set_aperture(0.03);
V1.add(HU_5, 0.0, 0.0, 9.06936269, reference='center');
% ---- Solenoid: H_1 ----
H_1 = Solenoid(1.56518469, get_peak_field(-535299999.99999994, 0.09, 0.03, 1.56518469), 0.03, 0.09, 5);
H_1.set_aperture(0.03);
V1.add(H_1, 0.0, 0.0, 9.857955035, reference='center');

% ---- Absorber: abs_1 ----
abs_1 = Absorber(1.32518469, 'liquid_hydrogen');
abs_1.set_aperture(0.025);
V1.add(abs_1, 0.0, 0.0, 9.857955035, reference='center');

% ---- Solenoid: HD_5 ----
HD_5 = Solenoid(0.012, get_peak_field(-535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5);
HD_5.set_aperture(0.03);
V1.add(HD_5, 0.0, 0.0, 10.64654738, reference='center');

% ---- Solenoid: HD_4 ----
HD_4 = Solenoid(0.012, get_peak_field(-535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5);
HD_4.set_aperture(0.03);
V1.add(HD_4, 0.0, 0.0, 10.65854738, reference='center');

% ---- Solenoid: HD_3 ----
HD_3 = Solenoid(0.012, get_peak_field(-535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5);
HD_3.set_aperture(0.03);
V1.add(HD_3, 0.0, 0.0, 10.67054738, reference='center');

% ---- Solenoid: MD_10 ----
MD_10 = Solenoid(0.2, get_peak_field(-28752919.330000002, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_10.set_aperture(0.1);
V1.add(MD_10, 0.0, 0.0, 10.87654738, reference='center');

% ---- Pillbox_Cavity: rf2_1 ----
rf2_1 = Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1);
rf2_1.set_phid(-90.0);
rf2_1.set_aperture(0.16);
V1.add(rf2_1, 0.0, 0.0, 11.15154738, reference='center');

% ---- Solenoid: MD_11 ----
MD_11 = Solenoid(0.2, get_peak_field(-5101607.766, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_11.set_aperture(0.1);
V1.add(MD_11, 0.0, 0.0, 11.42654738, reference='center');

% ---- Pillbox_Cavity: rf2_2 ----
rf2_2 = Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1);
rf2_2.set_phid(-90.0);
rf2_2.set_aperture(0.16);
V1.add(rf2_2, 0.0, 0.0, 11.70154738, reference='center');

% ---- Solenoid: MD_12 ----
MD_12 = Solenoid(0.2, get_peak_field(-10181668.370000001, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_12.set_aperture(0.1);
V1.add(MD_12, 0.0, 0.0, 11.97654738, reference='center');

% ---- Pillbox_Cavity: rf2_3 ----
rf2_3 = Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1);
rf2_3.set_phid(-90.0);
rf2_3.set_aperture(0.16);
V1.add(rf2_3, 0.0, 0.0, 12.25154738, reference='center');

% ---- Solenoid: MD_13 ----
MD_13 = Solenoid(0.2, get_peak_field(-5363069.454, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_13.set_aperture(0.1);
V1.add(MD_13, 0.0, 0.0, 12.52654738, reference='center');
% ---- Pillbox_Cavity: rf2_4 ----
rf2_4 = Pillbox_Cavity(25000000.0, 50000000.0, 0.25, 1);
rf2_4.set_phid(-90.0);
rf2_4.set_aperture(0.16);
V1.add(rf2_4, 0.0, 0.0, 12.80154738, reference='center');

% ---- Solenoid: MD_14 ----
MD_14 = Solenoid(0.2, get_peak_field(7639264.017, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_14.set_aperture(0.1);
V1.add(MD_14, 0.0, 0.0, 13.07654738, reference='center');
% ---- Pillbox_Cavity: rf2_5 ----
rf2_5 = Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1);
rf2_5.set_phid(-10.0);
rf2_5.set_aperture(0.16);
V1.add(rf2_5, 0.0, 0.0, 13.35154738, reference='center');

% ---- Solenoid: MD_15 ----
MD_15 = Solenoid(0.2, get_peak_field(-7899052.061, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_15.set_aperture(0.1);
V1.add(MD_15, 0.0, 0.0, 13.62654738, reference='center');

% ---- Pillbox_Cavity: rf2_6 ----
rf2_6 = Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1);
rf2_6.set_phid(-10.0);
rf2_6.set_aperture(0.16);
V1.add(rf2_6, 0.0, 0.0, 13.90154738, reference='center');

% ---- Solenoid: MD_16 ----
MD_16 = Solenoid(0.2, get_peak_field(5739264.072, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_16.set_aperture(0.1);
V1.add(MD_16, 0.0, 0.0, 14.17654738, reference='center');
% ---- Pillbox_Cavity: rf2_7 ----
rf2_7 = Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1);
rf2_7.set_phid(-10.0);
rf2_7.set_aperture(0.16);
V1.add(rf2_7, 0.0, 0.0, 14.45154738, reference='center');

% ---- Solenoid: MD_17 ----
MD_17 = Solenoid(0.2, get_peak_field(11250204.07, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_17.set_aperture(0.1);
V1.add(MD_17, 0.0, 0.0, 14.72654738, reference='center');
% ---- Pillbox_Cavity: rf2_8 ----
rf2_8 = Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1);
rf2_8.set_phid(-10.0);
rf2_8.set_aperture(0.16);
V1.add(rf2_8, 0.0, 0.0, 15.00154738, reference='center');

% ---- Solenoid: MD_18 ----
MD_18 = Solenoid(0.2, get_peak_field(5818381.93, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_18.set_aperture(0.1);
V1.add(MD_18, 0.0, 0.0, 15.27654738, reference='center');

% ---- Pillbox_Cavity: rf2_9 ----
rf2_9 = Pillbox_Cavity(18000000.0, 1000000.0, 0.25, 1);
rf2_9.set_phid(-10.0);
rf2_9.set_aperture(0.16);
V1.add(rf2_9, 0.0, 0.0, 15.55154738, reference='center');

% ---- Solenoid: MD_19 ----
MD_19 = Solenoid(0.2, get_peak_field(33814699.760000005, 0.3, 0.1, 0.2), 0.1, 0.3, 5);
MD_19.set_aperture(0.1);
V1.add(MD_19, 0.0, 0.0, 15.82654738, reference='center');

% ---- Solenoid: HU_6 ----
HU_6 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.135, 0.03, 0.012), 0.03, 0.135, 5);
HU_6.set_aperture(0.03);
V1.add(HU_6, 0.0, 0.0, 16.03254738, reference='center');

% ---- Solenoid: HU_7 ----
HU_7 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.125, 0.03, 0.012), 0.03, 0.125, 5);
HU_7.set_aperture(0.03);
V1.add(HU_7, 0.0, 0.0, 16.04454738, reference='center');

% ---- Solenoid: HU_8 ----
HU_8 = Solenoid(0.012, get_peak_field(535299999.99999994, 0.115, 0.03, 0.012), 0.03, 0.115, 5);
HU_8.set_aperture(0.03);
V1.add(HU_8, 0.0, 0.0, 16.05654738, reference='center');

% ---- Solenoid: H_2 ----
H_2 = Solenoid(1.34521782, get_peak_field(535299999.99999994, 0.09, 0.03, 1.34521782), 0.03, 0.09, 5);
H_2.set_aperture(0.03);
V1.add(H_2, 0.0, 0.0, 16.73515629, reference='center');
% ---- Absorber: abs_2 ----
abs_2 = Absorber(1.10521782, 'liquid_hydrogen');
abs_2.set_aperture(0.025);
V1.add(abs_2, 0.0, 0.0, 16.73515629, reference='center');

%% FINAL PART
% --------- VOLUMES ---------

V1.verbosity = 2; 
V1.dt_mm = 1;
V1.cfx_dt_mm = 10; 
V1.sc_dt_mm = 10;
V1.odeint_algorithm = 'rk2';
V1.odeint_epsabs = 1e-6;

C = Drift(0.1);
C.set_aperture(0.0);
V1.add(C, 0, 0, -0.001, 'exit');
V1.set_s0(0.0);
V1.set_s1(17.4); % End of tracking

M = B0.get_phase_space();
% Lets get the reference particle from B1
B_refp = Bunch6d(mass, charge, Q, mean(M));
B1_refp = V1.autophase(B_refp);


V1.set_tt_nsteps(250);
L.append(V1);


% --------- TRACKING ---------

% This is done for accelerating the autophase(). It does it for only the reference particle.
disp('Tracking...')

%% TRACKING
B1 = L.track(B0);
M_lost = L.get_lost_particles();
next_bunch = B1.get_phase_space();

S = L{1}.get_bunch_at_tt_screens();

n = numel(S);

mkdir("Phase_space_Bunch6d_SC");

for i = 1:n
    B_i = S{1, i};                 % extract Bunch6d object
    M_i = B_i.get_phase_space();   % get phase space matrix
    S_i = B_i.get_info().S;
    N_i = B_i.get_info().transmission;
    emitt_4d_i = B_i.get_info().emitt_4d;
    emitt_z_i = B_i.get_info().emitt_z;
    mean_K_i = B_i.get_info().mean_K;
    sigma_t_i = B_i.get_info().sigma_t;
    emitt_6d_i = B_i.get_info().emitt_6d;
    % Create the additional info row
    other_info = [i, S_i, N_i, emitt_4d_i, emitt_z_i, mean_K_i, sigma_t_i, emitt_6d_i];
    filename_info = sprintf("Phase_space_Bunch6d_SC/other_info_%03d.txt", i);
    dlmwrite(filename_info, other_info, "delimiter", "\t");

    filename = sprintf("Phase_space_Bunch6d_SC/M_%03d.txt", i);
    dlmwrite(filename, M_i, "delimiter", "\t");
end

% Transport table
T = L.get_transport_table("%S %mean_K %sigma_x %sigma_y %sigma_t %emitt_x %emitt_y %emitt_4d %emitt_z %N %beta_x %beta_y %alpha_x %alpha_y");
save('-text', 'transport_table_3cells_SC.txt', 'T');
%save('-text', 'lost.txt', 'M_lost');
