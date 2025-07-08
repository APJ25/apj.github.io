% main
clc;
clear all;
close all;
%% 
colors = [1 0 1; % Magenta
          1 0 0; % Red
          0 1 0; % Green
          0 0 1; % Blue
          0 1 1; % Cyan
          0 0 0; % Black
          1 0.5 0; % Orange
          0.5 0 0.5; % Purple
          0.5 0.5 0; % Olive
          0.2 0.6 1; % Light Blue
          1 0.2 0.2; % Coral
          0.3 0.3 0.3; % Dark Gray
          0.8 0.4 0; % Brown
         ];

ku  = 0.78; % Filling factor
Po  = 500; % Output power
Vpp = 80; Vps = 160; % Primary voltage for parallel and series
% dps = 0.644e-3;  
ds = 0.812e-3; % litz Wire diameter
dps=0.4*2e-3;
dpp = dps;
Vs1 = 24; Vs2 = 24;
f = 100e3; k = 0.160e3; alpha = 1.078; Beta = 2.92;
Vd = 0.7*2; D = 0.8; % Voltage drop and duty cycle
Ploss_max = 20;

I24 = Po/2/(Vs2 + 2*Vd);
Isrms = ((Vs1/Vps)*I24 + (Vs2/Vps)*I24) * sqrt(D)
Iprms = ((Vs1/Vpp)*I24 + (Vs2/Vpp)*I24) * sqrt(D)
I01rms = (1/2)*(I24*sqrt(1+D));
Is02rms = I01rms;
%% 

Awps = 2 * pi * (0.4e-3 / 2)^2;
Awpp = 2 * pi * (0.4e-3 / 2)^2;
Aws = 0.518e-6;
rho = 1.724e-8; muc = 1.256629e-6;
%Core parameters for 10 different cores selected for studies
%Cores=[EE10/11, EE13, EE16, EE19, EE25/19, EE30, EE40, EE50, EE60, ETD39,EE35,EE42,EE55]

% Core Parameters
lms  = [38.4 45.4 53.6 66.8 78.08 94.2 126.2 156.2 180.6 149 107.7 144.55 165.3]*1e-3;
Acs  = [5.39 8.82 10 11.96 20.1528 55 66 112.5 123.2 58.9 45.2675 118.69 176.571]*1e-6;
Was  = [46.62 68.6 82.8 111.6 158 151.6 328 524 814 514 93.7*2 191.64*2 331.82*2]*(1e-6)/2;
MLTs = [12 15.4 14 15 19.18 32 34 45 47.4 40.21 29.32 42.45 51.46]*1e-3;
hs   = [4.2*2 4.6*2 5.0*2 5.6*2 6.41*2 8.15*2 10.25*2 12.75*2 14.05*2 14.6*2 9.71*2 15.15*2 18.8*2]*1e-3 - 4e-3;
Vse  = lms .* Acs;

lmp  = [52.2 60.4 69 78.8 97.4 115.4 154.6 191.6 220 184.2 127.12 174.85 202.9]*1e-3;
Acp  = Acs;
Wap  = Was;
MLTp = MLTs;
hp   = hs;
Vsp  = lmp .* Acp;

Ss = Awps + Aws;
Sp = Awpp + Aws;
ns = 3; np = 3;

Tamb = 25; Tmax = 75;
% Rth = [63.8 62.5 76 60 32.5 23 20 13 6.5 15.3 18 10.31 7.674];
Rths = [63.8 62.5 76 60 41.49 21.80 16.87 11.2 9.925 16.40 22.53 11.42 8.572];
Rthp=  [63.8 62.5 76 60 36.82 19.54 15.12 10.10 8.921 14.63 20.63 10.31 7.674];
Pmaxs = (Tmax - Tamb) ./ Rths;
Pmaxp = (Tmax - Tamb) ./ Rthp;
Bsmax = (Pmaxs ./ (k .* f.^alpha .* Vse)).^(1/Beta);
Bpmax = (Pmaxp ./ (k .* f.^alpha .* Vsp)).^(1/Beta);

Npsmax = (ku * Was) ./ (Awps + Aws / ns);
Nppmax = (ku * Wap) ./ (Awpp + Aws / np);
Npsmin = (0.5 .* Vps * D) ./ (Acs .* f * 2 .* Bsmax);
Nppmin = (Vpp * D) ./ (Acp .* f * 2 .* Bpmax);

core_names = {'EE10/11', 'EE13', 'EE16', 'EE19', 'EE25/19', 'EE30', 'EE40', 'EE50', 'EE60', 'ETD39', 'EE35', 'EE42', 'EE55'};
num_cores = numel(core_names);

Turns.Series.Primary = cell(1, num_cores);
Turns.Parallel.Primary = cell(1, num_cores);
Flux.Series = cell(1, num_cores);
Flux.Parallel = cell(1, num_cores);
Pfes = cell(1, num_cores);
Pfep = cell(1, num_cores);
Turns.Series.Secondary = cell(1, num_cores);
Turns.Parallel.Secondary = cell(1, num_cores);
Layers.Series.Primary = cell(1, num_cores);
Layers.Series.Secondary = cell(1, num_cores);
Layers.Parallel.Primary = cell(1, num_cores);
Layers.Parallel.Secondary = cell(1, num_cores);
Lfactor.Series.Primary = cell(1, num_cores);
Lfactor.Series.Secondary = cell(1, num_cores);
Lfactor.Parallel.Primary = cell(1, num_cores);
Lfactor.Parallel.Secondary = cell(1, num_cores);

nus = (sqrt(pi) * dps) / (2 * (dps + 0.051e-3));
nup = (sqrt(pi) * dpp) / (2 * (dpp + 0.051e-3));
Sd = sqrt(rho * (1 + 0.0393 * (Tmax - 25)) / (pi * f * muc));
deltaps = (pi/4)^0.75 * (dps / Sd) * sqrt(nus);
deltapp = (pi/4)^0.75 * (dpp / Sd) * sqrt(nup);
deltas = (pi/4)^0.75 * (ds / Sd) * sqrt(nup);
Skeffs = (sinh(2*deltaps) + sin(2*deltaps)) / (cosh(2*deltaps) - cos(2*deltaps));
Skeffp = (sinh(2*deltapp) + sin(2*deltapp)) / (cosh(2*deltapp) - cos(2*deltapp));
Proxeffs = (sinh(deltaps) - sin(deltaps)) / (cosh(deltaps) + cos(deltaps));
Proxeffp = (sinh(deltapp) - sin(deltapp)) / (cosh(deltapp) + cos(deltapp));

Fr.Series.Primary = cell(1, num_cores);
Fr.Series.Secondary = cell(1, num_cores);
Fr.Parallel.Primary = cell(1, num_cores);
Fr.Parallel.Secondary = cell(1, num_cores);
Resistance.Series.Primary = cell(1, num_cores);
Resistance.Series.Secondary = cell(1, num_cores);
Resistance.Parallel.Primary = cell(1, num_cores);
Resistance.Parallel.Secondary = cell(1, num_cores);
CopperLoss.Series = cell(1, num_cores);
CopperLoss.Parallel = cell(1, num_cores);
PowerLoss.Series = cell(1, num_cores);
PowerLoss.Parallel = cell(1, num_cores);

for i = 1:num_cores
    % Calculate Np ranges
Nps_range = (Npsmin(i)):(Npsmax(i));
Npp_range = (Nppmin(i)):(Nppmax(i));

% Calculate B for both series and parallel
B_series_full = ((0.5 * Vps * D) ./ (2 * Acs(i) * f .* Nps_range));
B_parallel_full = (Vpp * D) ./ (2 * Acp(i) * f .* Npp_range);

% Apply flux density limits
valid_idx_series = B_series_full <= Bsmax(i);
valid_idx_parallel = B_parallel_full <= Bpmax(i);

% Filter based on limits
B_series = B_series_full(valid_idx_series);
B_parallel = B_parallel_full(valid_idx_parallel);

Nps_valid = Nps_range(valid_idx_series);
Npp_valid = Npp_range(valid_idx_parallel);

% Save to structures
Turns.Series.Primary{i} = Nps_valid;
Turns.Parallel.Primary{i} = Npp_valid;
Flux.Series{i} = B_series;
Flux.Parallel{i} = B_parallel;

Turns.Series.Secondary{i} = Nps_valid / ns;
Turns.Parallel.Secondary{i} = Npp_valid / np;


    Pfes{i} = k .* f.^alpha .* B_series.^Beta .* Acs(i) .* lms(i);
    Pfep{i} = k .* f.^alpha .* B_parallel.^Beta .* Acp(i) .* lmp(i);

    Layers.Series.Primary{i} = max((0.5 .* Turns.Series.Primary{i} .* dps) ./ hs(i), 1);
    Layers.Series.Secondary{i} = max((0.5 .* Turns.Series.Secondary{i} .* ds) ./ hs(i), 1);
    Layers.Parallel.Primary{i} = max((0.5 .* Turns.Parallel.Primary{i} .* dpp) ./ hp(i), 1);
    Layers.Parallel.Secondary{i} = max((0.5 .* Turns.Parallel.Secondary{i} .* ds) ./ hp(i), 1);

    Lfactor.Series.Primary{i} = (2 .* (Layers.Series.Primary{i}.^2) - 1) / 3;
    Lfactor.Series.Secondary{i} = (2 .* (Layers.Series.Secondary{i}.^2) - 1) / 3;
    Lfactor.Parallel.Primary{i} = (2 .* (Layers.Parallel.Primary{i}.^2) - 1) / 3;
    Lfactor.Parallel.Secondary{i} = (2 .* (Layers.Parallel.Secondary{i}.^2) - 1) / 3;

    Fr.Series.Primary{i} = deltaps .* (Skeffs + Lfactor.Series.Primary{i} .* Proxeffs);
    Fr.Series.Secondary{i} = deltas  .* (Skeffs + Lfactor.Series.Secondary{i} .* Proxeffs);
    Fr.Parallel.Primary{i} = deltapp .* (Skeffp + Lfactor.Parallel.Primary{i} .* Proxeffp);
    Fr.Parallel.Secondary{i} = deltas .* (Skeffp + Lfactor.Parallel.Secondary{i} .* Proxeffp);

    Resistance.Series.Primary{i} = Fr.Series.Primary{i} .* Turns.Series.Primary{i} .* (rho / Awps) .* MLTs(i);
    Resistance.Series.Secondary{i} = Fr.Series.Secondary{i} .* Turns.Series.Secondary{i} .* (rho / Aws) .* MLTs(i);
    Resistance.Parallel.Primary{i} = Fr.Parallel.Primary{i} .* Turns.Parallel.Primary{i} .* (rho / Awpp) .* MLTp(i);
    Resistance.Parallel.Secondary{i} = Fr.Parallel.Secondary{i} .* Turns.Parallel.Secondary{i} .* (rho / Aws) .* MLTp(i);

    % Copper Losses 
    CopperLoss.Series{i} = Resistance.Series.Primary{i} .* Isrms.^2 + Resistance.Series.Secondary{i} .* I01rms.^2;
    CopperLoss.Parallel{i} = 0.5 .* Resistance.Parallel.Primary{i} .* (0.5 .* Iprms).^2 + Resistance.Parallel.Secondary{i} .* I01rms.^2;

    % Total Power Loss (Copper + Core)
    PowerLoss.Series{i} = CopperLoss.Series{i} + Pfes{i};
    PowerLoss.Parallel{i} = CopperLoss.Parallel{i} + Pfep{i};

end

%% Fill Factor and Ku Calculation (Series Only)

t_ins = 0.051e-3;  % Interlayer insulation thickness (m)

% AWG 22 parameters
d_barep = 0.4e-3;        % Bare wire diameter (m)
d_insulatedp = 0.442e-3;   % Outer diameter with heavy build insulation (m)
A_copperp = 2*pi * (d_barep / 2)^2;
A_total_wirep = pi * (d_insulatedp / 2)^2;
d_bares = 0.812e-3;        % Bare wire diameter (m)
d_insulateds = 1.082e-3;   % Outer diameter with heavy build insulation (m)
A_coppers = pi * (d_bares / 2)^2;
A_total_wires = pi * (d_insulateds / 2)^2;
A_total_wire=A_total_wirep+A_total_wires;
S1 = (A_copperp+A_coppers) / A_total_wire;  % Constant for AWG 22

for i = 1:num_cores
    FillFactor.Series{i} = zeros(1, length(Turns.Series.Primary{i}));
    Ku.Series{i} = zeros(1, length(Turns.Series.Primary{i}));

    % Parameters for S3 (constant per core)
    h = hs(i);
    margin = 1.5e-3;%bobbin
    top_clearance = 1e-3;%bobbin
    bottom_clearance = 1e-3;%bobbin

    usable_width = Was(i) / h;
    usable_width_eff = usable_width - margin;%only one side of window the bobbin is present
    usable_height_eff = h - top_clearance - bottom_clearance;
    A_usable_window = usable_width_eff * usable_height_eff;

    S3 = (usable_width_eff * usable_height_eff) / (usable_width * h);

    for j = 1:length(Turns.Series.Primary{i})
        Np = Turns.Series.Primary{i}(j);
        Ns = Turns.Series.Secondary{i}(j);
        Nlp = Layers.Series.Primary{i}(j);
        Nls = Layers.Series.Secondary{i}(j);

        % === S2 Calculation ===layer
        Awp_total = Np * Awps;
        Aws_total = Ns * Aws;
        A_copper_total = Awp_total + Aws_total;

        
        A_allocated = Was(i)-(Was(i)-A_copper_total);%available window area - residual area, that results from the particular winding
%technique used.
        S2_dynamic = A_copper_total / A_allocated;
       layer_penalty = 1 - 0.02 * (Nlp + Nls - 1);  % 2% loss per layer 
       S2 = max(min(S2_dynamic * layer_penalty, 1), 0);


        % === S4 Calculation ===
        A_ins = (Nlp + Nls - 1) * t_ins * h;
        S4 = A_usable_window / (A_usable_window + A_ins);

        % === Ku and Fill Factor ===
        ku_dynamic = S1 * S2 * S3 * S4;
        Ku.Series{i}(j) = ku_dynamic;

        
        FillFactor.Series{i}(j) = ku_dynamic;
    end
end
figure(17);
axis square
hold on; grid on;

for i = 1:num_cores
    if ~isempty(Ku.Series{i})
        color_idx = mod(i-1, size(colors,1)) + 1;
        plot(Turns.Series.Primary{i}, Ku.Series{i}, ...
            'LineWidth', 2, ...
            'Color', colors(color_idx, :), ...
            'DisplayName', core_names{i});
    end
end

xlabel('Number  Turns ');
ylabel('Window Utilization Factor (Ku)');
% title('Ku vs Number of Primary Turns (Series)');
legend('Location', 'eastoutside');
set(gca, 'FontSize', 12);

%% 
%optimal points calculation
for i = 1:num_cores
    % --- Series Optimal ---
    if ~isempty(PowerLoss.Series{i})
        [opt_series_val, Iser] = min(PowerLoss.Series{i});
        Optimal.Series.Power(i) = opt_series_val;
        Optimal.Series.B(i)     = Flux.Series{i}(Iser);
        Optimal.Series.Np(i)    = Turns.Series.Primary{i}(Iser);
    else
        Optimal.Series.Power(i) = NaN;
        Optimal.Series.B(i)     = NaN;
        Optimal.Series.Np(i)    = NaN;
    end

    % --- Parallel Optimal ---
    if ~isempty(PowerLoss.Parallel{i})
        [opt_parallel_val, Ipar] = min(PowerLoss.Parallel{i});
        Optimal.Parallel.Power(i) = opt_parallel_val;
        Optimal.Parallel.B(i)     = Flux.Parallel{i}(Ipar);
        Optimal.Parallel.Np(i)    = Turns.Parallel.Primary{i}(Ipar);
    else
        Optimal.Parallel.Power(i) = NaN;
        Optimal.Parallel.B(i)     = NaN;
        Optimal.Parallel.Np(i)    = NaN;
    end
end




%% Plots
% Specify core indices for group plots
core_indices = 5:13;

% Plot 1: Core Loss vs Flux Density
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 1, "loss_vs_flux", Fr, CopperLoss, PowerLoss, Optimal, Pmaxs);

% Plot 2: Core Loss Comparison for EE50-Z (index 8)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, 8, 2, "compare_core", Fr, CopperLoss, PowerLoss, Optimal, Pmaxs);

% Plot 3: Core Loss vs Turns
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 3, "loss_vs_turns", Fr, CopperLoss, PowerLoss, Optimal, Pmaxs);

% Plot 4: Flux Density vs Turns
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 4, "flux_vs_turns", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);

% Plot 5: Core Loss Comparison for EE42/19 (index 12)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, 12, 5, "compare_core", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);

% Plot 6: Number of Layers vs Turns (Series & Parallel)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 6, "layers_vs_turns", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);

% Plot 7: Layers Comparison for EE50-Z
%plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, 8, 7, "compare_layers", Fr, CopperLoss, PowerLoss, Optimal, Pmax);

% Plot 8: Skeff & Proximity Effects for EE50 (index 8)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, 8, 8, "skeff proximity", Fr, CopperLoss, PowerLoss, Optimal, Pmaxs);

% Plot 9: Copper Loss vs Turns
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 9, "copper loss vs turns", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);
%% 

% Plot 10: Total Power Loss vs Flux Density (Series  + Opt markers)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 10, "Total power loss vs flux density series", Fr, CopperLoss, PowerLoss, Optimal, Pmaxs);
% Plot 11: Total Power Loss vs Flux Density (Parallel  + Opt markers)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 11, "Total power loss vs flux density parallel", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);
% Plot 12: Total Power Loss vs Turns (Series  + Opt markers)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 12, "Total power loss vs turns series", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);
% Plot 13: Total Power Loss vs Turns (Parallel  + Opt markers)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 13, "Total power loss vs turns parallel", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);
% Plot 14: Total Power Loss vs Flux Density (Parallel  + series)
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 14, "Total power loss vs turns combined", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);
% Plot 15: Total Power Loss vs Flux Density (Parallel  + series)
% plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 15, "optimal turns per core", Fr, CopperLoss, PowerLoss, Optimal, Pmax);
% Plot 16: Total turns,error bar
plotCoreLosses(Flux, Pfes, Pfep, Turns, Layers, core_names, colors, core_indices, 16, "turns range errorbar", Fr, CopperLoss, PowerLoss, Optimal, Pmaxp);

%% 
% Create table of Pmax vs Optimal Loss Differences
Core = core_names(:);
Pmax_Series = Pmaxs(:);
Opt_P_Series = Optimal.Series.Power(:);
Diff_Series = Pmax_Series - Opt_P_Series;

Pmax_Parallel = Pmaxp(:);
Opt_P_Parallel = Optimal.Parallel.Power(:);
Diff_Parallel = Pmax_Parallel - Opt_P_Parallel;

% Construct table
loss_diff_table = table(Core, Pmax_Series, Opt_P_Series, Diff_Series, ...
                        Pmax_Parallel, Opt_P_Parallel, Diff_Parallel);

% Display or export the table
disp(loss_diff_table);
% Optionally, write to file:
% writetable(loss_diff_table, 'loss_difference_summary.csv');


