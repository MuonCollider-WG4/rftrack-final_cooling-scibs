clear; clc;
RF_Track;

mass = RF_Track.muonmass;

folder_0 = "Phase_space_Bunch6d";
folder_SC = "Phase_space_Bunch6d_SC";
folder_SC_IBS = "Phase_space_Bunch6d_SC_IBS";

cut = 3;

% ----------------------------
% Preallocate
% ----------------------------
mean_t0         = NaN(249,1);
mean_t_SC       = NaN(249,1);
mean_t_SC_IBS   = NaN(249,1);

emitt4d         = NaN(249,1);
emitt4d_SC      = NaN(249,1);
emitt4d_SC_IBS  = NaN(249,1);

emittz         = NaN(249,1);
emittz_SC      = NaN(249,1);
emittz_SC_IBS  = NaN(249,1);

emitt6d         = NaN(249,1);
emitt6d_SC      = NaN(249,1);
emitt6d_SC_IBS  = NaN(249,1);

ES_arr      = NaN(249,1);
ES_SC_arr   = NaN(249,1);
ES_SC_IBS_arr = NaN(249,1);

length_arr      = NaN(249,1);
length_SC_arr   = NaN(249,1);
length_SC_IBS_arr = NaN(249,1);

mean_K_arr      = NaN(249,1);
mean_K_SC_arr   = NaN(249,1);
mean_K_SC_IBS_arr = NaN(249,1);

sigma_t_arr      = NaN(249,1);
sigma_t_SC_arr   = NaN(249,1);
sigma_t_SC_IBS_arr = NaN(249,1);

N_arr            = NaN(249,1);
N_SC_arr         = NaN(249,1);
N_SC_IBS_arr     = NaN(249,1);

struct_M = cell(1,249);
struct_M_SC = cell(1,249);
struct_M_SC_IBS = cell(1,249);

allX = [];
allY = [];

% ----------------------------
% Main loop
% ----------------------------
for i = 1:249

    filename        = sprintf("%s/M_%03d.txt", folder_0, i);
    filename_SC     = sprintf("%s/M_%03d.txt", folder_SC, i);
    filename_SC_IBS = sprintf("%s/M_%03d.txt", folder_SC_IBS, i);

    other_info        = sprintf("%s/other_info_%03d.txt", folder_0, i);
    other_info_SC     = sprintf("%s/other_info_%03d.txt", folder_SC, i);
    other_info_SC_IBS = sprintf("%s/other_info_%03d.txt", folder_SC_IBS, i);

    % ---- file check ----
    if ~exist(filename, "file") || ...
       ~exist(filename_SC, "file") || ...
       ~exist(filename_SC_IBS, "file") || ...
       ~exist(other_info, "file") || ...
       ~exist(other_info_SC, "file") || ...
       ~exist(other_info_SC_IBS, "file")

        fprintf("Skipping i=%d (missing file)\n", i);
        continue;
    end

    info        = load(other_info);
    info_SC     = load(other_info_SC);
    info_SC_IBS = load(other_info_SC_IBS);

    if isempty(info) || isempty(info_SC) || isempty(info_SC_IBS)
        fprintf("Skipping i=%d (empty info)\n", i);
        continue;
    end

    length_arr(i)        = info(1,2);
    N_arr(i)             = info(1,3);
    mean_K_arr(i)        = info(1,6);
    sigma_t_arr(i)       = info(1,7);

    length_SC_arr(i)     = info_SC(1,2);
    N_SC_arr(i)          = info_SC(1,3);
    mean_K_SC_arr(i)     = info_SC(1,6);
    sigma_t_SC_arr(i)    = info_SC(1,7);

    length_SC_IBS_arr(i) = info_SC_IBS(1,2);
    N_SC_IBS_arr(i)      = info_SC_IBS(1,3);
    mean_K_SC_IBS_arr(i) = info_SC_IBS(1,6);
    sigma_t_SC_IBS_arr(i)= info_SC_IBS(1,7);

    N_arr_ini        = 5.3e12; % initial bunch charge
    N_SC_arr_ini     = 5.3e12;
    N_SC_IBS_arr_ini = 5.3e12;

    % ----------------------------
    % Load phase space
    % ----------------------------
    M        = load(filename);
    M_SC     = load(filename_SC);
    M_SC_IBS = load(filename_SC_IBS);

    x  = M(:,1);   xp = M(:,2);
    y  = M(:,3);   yp = M(:,4);
    t  = M(:,5) * 3.335641e-09;   P  = M(:,6);
    E = sqrt(P.^2 + mass^2) * 1e6; 

    x_SC  = M_SC(:,1);   xp_SC = M_SC(:,2);
    y_SC  = M_SC(:,3);   yp_SC = M_SC(:,4);
    t_SC  = M_SC(:,5) * 3.335641e-09;   P_SC  = M_SC(:,6);
    E_SC = sqrt(P_SC.^2 + mass^2) * 1e6; 

    x_IBS  = M_SC_IBS(:,1); xp_IBS = M_SC_IBS(:,2);
    y_IBS  = M_SC_IBS(:,3); yp_IBS = M_SC_IBS(:,4);
    t_IBS  = M_SC_IBS(:,5) * 3.335641e-09; P_IBS  = M_SC_IBS(:,6);
    E_IBS = sqrt(P_IBS.^2 + mass^2) * 1e6; 

    struct_M{i} = M;
    struct_M_SC{i} = M_SC;
    struct_M_SC_IBS{i} = M_SC_IBS;
    allX = [allX; struct_M{i}(:,1); struct_M_SC{i}(:,1); struct_M_SC_IBS{i}(:,1)];
    allY = [allY; struct_M{i}(:,2); struct_M_SC{i}(:,2); struct_M_SC_IBS{i}(:,2)];

    % ----------------------------
    % Statistics
    % ----------------------------
    mx = mean(x);   sx = std(x);
    mxp = mean(xp); sxp = std(xp);
    my = mean(y);   sy = std(y);
    myp = mean(yp); syp = std(yp);
    mt = mean(t);   st = std(t);
    mE = mean(E);   sE = std(E);

    mx_SC = mean(x_SC); sx_SC = std(x_SC);
    mxp_SC = mean(xp_SC); sxp_SC = std(xp_SC);
    my_SC = mean(y_SC); sy_SC = std(y_SC);
    myp_SC = mean(yp_SC); syp_SC = std(yp_SC);
    mt_SC = mean(t_SC); st_SC = std(t_SC);
    mE_SC = mean(E_SC); sE_SC = std(E_SC);

    mx_IBS = mean(x_IBS); sx_IBS = std(x_IBS);
    mxp_IBS = mean(xp_IBS); sxp_IBS = std(xp_IBS);
    my_IBS = mean(y_IBS); sy_IBS = std(y_IBS);
    myp_IBS = mean(yp_IBS); syp_IBS = std(yp_IBS);
    mt_IBS = mean(t_IBS); st_IBS = std(t_IBS);
    mE_IBS = mean(E_IBS); sE_IBS = std(E_IBS);

    ES_arr(i)      = sE;
    ES_SC_arr(i)   = sE_SC;
    ES_SC_IBS_arr(i) = sE_IBS;
    % ----------------------------
    % Cut selection
    % ----------------------------
    mask = (abs(x  - mx)  < cut*sx)  & ...
           (abs(xp - mxp) < cut*sxp) & ...
           (abs(y  - my)  < cut*sy)  & ...
           (abs(yp - myp) < cut*syp) & ...
           (abs(t  - mt)  < cut*st)  & ...
           (abs(E  - mE)  < cut*sE);

    mask_SC = (abs(x_SC  - mx_SC)  < cut*sx_SC)  & ...
            (abs(xp_SC - mxp_SC) < cut*sxp_SC) & ...
            (abs(y_SC  - my_SC)  < cut*sy_SC)  & ...
            (abs(yp_SC - myp_SC) < cut*syp_SC) & ...
            (abs(t_SC  - mt_SC)  < cut*st_SC)  & ...
            (abs(E_SC  - mE_SC)  < cut*sE_SC);

    mask_IBS = (abs(x_IBS  - mx_IBS)  < cut*sx_IBS)  & ...
            (abs(xp_IBS - mxp_IBS) < cut*sxp_IBS) & ...
            (abs(y_IBS  - my_IBS)  < cut*sy_IBS)  & ...
            (abs(yp_IBS - myp_IBS) < cut*syp_IBS) & ...
            (abs(t_IBS  - mt_IBS)  < cut*st_IBS)  & ...
            (abs(E_IBS  - mE_IBS)  < cut*sE_IBS);

    X     = [x(mask) xp(mask) y(mask) yp(mask)];
    X_SC  = [x_SC(mask_SC) xp_SC(mask_SC) y_SC(mask_SC) yp_SC(mask_SC)];
    X_IBS = [x_IBS(mask_IBS) xp_IBS(mask_IBS) y_IBS(mask_IBS) yp_IBS(mask_IBS)];

    if size(X,1) < 10 || size(X_SC,1) < 10 || size(X_IBS,1) < 10
        continue;
    end

    Z     = [t(mask)     E(mask)];
    Z_SC  = [t_SC(mask_SC) E_SC(mask_SC)];
    Z_IBS = [t_IBS(mask_IBS) E_IBS(mask_IBS)];

    if size(Z,1) < 10 || size(Z_SC,1) < 10 || size(Z_IBS,1) < 10
        continue;
    end

    % Covariance determinants
    detZ     = det(cov(Z));
    detZ_SC  = det(cov(Z_SC));
    detZ_IBS = det(cov(Z_IBS));
    % Longitudinal emittance
    emittz(i)        = sqrt(abs(detZ));
    emittz_SC(i)     = sqrt(abs(detZ_SC));
    emittz_SC_IBS(i) = sqrt(abs(detZ_IBS));
    % ----------------------------
    % Emittance
    % ----------------------------
    emitt4d(i)        = abs(det(cov(X)))^(1/4)       * (mean(P(mask)) / mass);
    emitt4d_SC(i)     = abs(det(cov(X_SC)))^(1/4)    * (mean(P_SC(mask_SC)) / mass);
    emitt4d_SC_IBS(i) = abs(det(cov(X_IBS)))^(1/4)   * (mean(P_IBS(mask_IBS)) / mass);

end


xmin = min(allX);
xmax = max(allX);
ymin = min(allY);
ymax = max(allY);

% Plot
% ----------------------------
figure(1); clf;
plot(length_arr, emitt4d, 'k', "linewidth", 1.5); hold on;
plot(length_SC_arr, emitt4d_SC, 'b', "linewidth", 1.5);
plot(length_SC_IBS_arr, emitt4d_SC_IBS, 'r', "linewidth", 1.5);
xlabel("z [m]", 'FontSize', 20);
xlim([0 17.4])
ylabel("4D normalized emittance", 'FontSize', 20);
legend("No SC, No IBS", "SC", "SC + IBS", 'FontSize', 20);
grid on;
title("Three first cells of the Muon Final Cooling", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(2); clf;
plot(length_arr, N_arr / N_arr_ini * 100, 'k', "linewidth", 1.5); hold on;
plot(length_SC_arr, N_SC_arr / N_SC_arr_ini * 100, 'b', "linewidth", 1.5);
plot(length_SC_IBS_arr, N_SC_IBS_arr / N_SC_IBS_arr_ini * 100, 'r', "linewidth", 1.5);
xlabel("z [m]", 'FontSize', 20);
ylabel("Transmission (%)", 'FontSize', 20);
xlim([0 17.4])
legend("No SC, No IBS", "SC", "SC + IBS", 'FontSize', 20);
grid on;
title("Three first cells of the Muon Final Cooling", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(3); clf;
plot(length_arr, emittz, 'k', "linewidth", 1.5); hold on;
plot(length_SC_arr, emittz_SC, 'b', "linewidth", 1.5);
plot(length_SC_IBS_arr, emittz_SC_IBS, 'r', "linewidth", 1.5);
xlabel("z [m]", 'FontSize', 20);
xlim([0 17.4])
ylabel("Normalized longitudinal emittance [eV * ms]", 'FontSize', 20);
legend("No SC, No IBS", "SC", "SC + IBS", 'FontSize', 20, 'Location', 'northwest');
grid on;
title("Three first cells of the Muon Final Cooling", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(4); clf;
plot(length_arr, mean_K_arr, 'k', "linewidth", 1.5); hold on;
plot(length_SC_arr, mean_K_SC_arr, 'b', "linewidth", 1.5);
plot(length_SC_IBS_arr, mean_K_SC_IBS_arr, 'r', "linewidth", 1.5);
xlabel("z [m]", 'FontSize', 20);
xlim([0 17.4])
ylabel("Kinetic energy [MeV]", 'FontSize', 20);
legend("No SC, No IBS", "SC", "SC + IBS", 'FontSize', 20, 'Location', 'northwest');
grid on;
title("Three first cells of the Muon Final Cooling", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(5); clf;
plot(length_arr, sigma_t_arr, 'k', "linewidth", 1.5); hold on;
plot(length_SC_arr, sigma_t_SC_arr, 'b', "linewidth", 1.5);
plot(length_SC_IBS_arr, sigma_t_SC_IBS_arr, 'r', "linewidth", 1.5);
xlabel("z [m]", 'FontSize', 20);
xlim([0 17.4])
ylabel("Bunch length [mm/c]", 'FontSize', 20);
legend("No SC, No IBS", "SC", "SC + IBS", 'FontSize', 20, 'Location', 'northwest');
grid on;
title("Three first cells of the Muon Final Cooling", 'FontSize', 22);
set(gca, 'FontSize', 20);

%{
%%%%%%%%%%%%%%%%%%%% PHASE SPACE %%%%%%%%%%%%%%%%%%%%%%%
%% X AND XP initial
figure(6)
scatter(M_ini_SC_IBS(:,1), M_ini_SC_IBS(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.9);
hold on
scatter(M_ini_SC(:,1), M_ini_SC(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.6);
hold on
scatter(M_ini(:,1), M_ini(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel("x [mm]", 'FontSize', 20);
ylabel("xp [mm * mrad]", 'FontSize', 20);
legend("SC + IBS", "SC", "No SC, No IBS",'FontSize', 20);
grid on;
title("Initial horizontal phase space", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(7)
scatter(M_end_SC_IBS(:,1), M_end_SC_IBS(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.9);
hold on
scatter(M_end_SC(:,1), M_end_SC(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.6);
hold on
scatter(M_end(:,1), M_end(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel("x [mm]", 'FontSize', 20);
ylabel("xp [mm * mrad]", 'FontSize', 20);
legend("SC + IBS", "SC", "No SC, No IBS", 'FontSize', 20);
grid on;
title("Final horizontal phase space", 'FontSize', 22);
set(gca, 'FontSize', 20);

%%%%%%%%%%%%%%%%%%%% PHASE SPACE %%%%%%%%%%%%%%%%%%%%%%%
figure(8)
scatter(M_ini_SC_IBS(:,3), M_ini_SC_IBS(:,4), 10, 'filled', 'MarkerFaceAlpha', 0.9);
hold on
scatter(M_ini_SC(:,3), M_ini_SC(:,4), 10, 'filled', 'MarkerFaceAlpha', 0.6);
hold on
scatter(M_ini(:,3), M_ini(:,4), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel("y [mm]", 'FontSize', 20);
ylabel("yp [mm * mrad]", 'FontSize', 20);
legend("SC + IBS", "SC", "No SC, No IBS", 'FontSize', 20);
grid on;
title("Initial vertical phase space", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(9)
scatter(M_end_SC_IBS(:,3), M_end_SC_IBS(:,4), 10, 'filled', 'MarkerFaceAlpha', 0.9);
hold on
scatter(M_end_SC(:,3), M_end_SC(:,4), 10, 'filled', 'MarkerFaceAlpha', 0.6);
hold on
scatter(M_end(:,3), M_end(:,4), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel("y [mm]", 'FontSize', 20);
ylabel("yp [mm * mrad]", 'FontSize', 20);
legend("SC + IBS", "SC", "No SC, No IBS", 'FontSize', 20);
grid on;
title("Final vertical phase space", 'FontSize', 22);
set(gca, 'FontSize', 20);

%%%%%%%%%%%%%%%%%%%% PHASE SPACE - Not done like this anymore %%%%%%%%%%%%%%%%%%%%%%%
figure(10)
scatter(M_ini_SC_IBS(:,5), M_ini_SC_IBS(:,6), 10, 'filled', 'MarkerFaceAlpha', 0.9);
hold on
scatter(M_ini_SC(:,5), M_ini_SC(:,6), 10, 'filled', 'MarkerFaceAlpha', 0.6);
hold on
scatter(M_ini(:,5), M_ini(:,6), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel("t [mm/c]", 'FontSize', 20);
ylabel("MeV [P]", 'FontSize', 20);
legend("SC + IBS", "SC", "No SC, No IBS", 'FontSize', 20);
grid on;
title("Initial longitudinal phase space", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(11)
scatter(M_end_SC_IBS(:,5), M_end_SC_IBS(:,6), 10, 'filled', 'MarkerFaceAlpha', 0.9);
hold on
scatter(M_end_SC(:,5), M_end_SC(:,6), 10, 'filled', 'MarkerFaceAlpha', 0.6);
hold on
scatter(M_end(:,5), M_end(:,6), 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel("t [mm/c]", 'FontSize', 20);
ylabel("P [MeV]", 'FontSize', 20);
legend("SC + IBS", "SC", "No SC, No IBS", 'FontSize', 20);
grid on;
title("Final longitudinal phase space", 'FontSize', 22);
set(gca, 'FontSize', 20);
%}

figure(8); clf;
plot(length_arr, ES_arr / 1e6, 'k', "linewidth", 1.5); hold on;
plot(length_SC_arr, ES_SC_arr / 1e6, 'b', "linewidth", 1.5);
plot(length_SC_IBS_arr, ES_SC_IBS_arr / 1e6, 'r', "linewidth", 1.5);
xlabel("z [m]", 'FontSize', 20);
ylabel("Energy spread (MeV)", 'FontSize', 20);
xlim([0 17.4])
legend("No SC, No IBS", "SC", "SC + IBS", 'FontSize', 20);
grid on;
title("Three first cells of the Muon Final Cooling", 'FontSize', 22);
set(gca, 'FontSize', 20);

figure(9);
for i = 1:249
    clf;

    M = struct_M{i};
    M_SC = struct_M_SC{i};
    M_IBS = struct_M_SC_IBS{i};

    scatter(M_IBS(:,1), M_IBS(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.9);
    hold on;
    scatter(M_SC(:,1), M_SC(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(M(:,1), M(:,2), 10, 'filled', 'MarkerFaceAlpha', 0.3);

    axis([xmin xmax ymin ymax]);
    axis equal;

    title(sprintf('Frame = %.4f', length_arr(i)));
    drawnow;
    pause(0.1);
end

