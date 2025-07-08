clc;
clear all;
close all;
sympref('FloatingPointOutput',true);
format shortEng
syms eg N  % air gap and number of turns

% Steinmetz parameters
a1 = -0.1108;
b1 = 0.2134;
c1 = 4.213;
a2 = 0.215;
b2 = 0.1373;
c2 = 0.03359;
a3 = 1.022e+06;
b3 = -1.348;
c3 = -0.02106;

f = 100e3;
alpha = (a2 * f^b2 + c2);
Beta = (a1 * f^b1 + c1);
k = a3 * f^b3 + c3;

%%%% Input data %%%%%%%
Dc_coarse = 4e-3:0.5e-3:14e-3;  % coarse diameter values
s = 100;                        % number of fine points
Dc = interp1(1:length(Dc_coarse), Dc_coarse, linspace(1, length(Dc_coarse), s), 'pchip');

% Preallocate
valueOfeg = zeros(1,s);
valueOfN = zeros(1,s);
A = zeros(1,s);
L1 = zeros(1,s);
e = zeros(1,s);
e1 = zeros(1,s);

for i = 1:s
    rau = 1.72e-8;
    I = 8;
    uo = 4*pi*1e-7;
    ur = 2000;
    B = 0.2;
    lc = 0e-3;

    A(i) = pi * Dc(i)^2 / 4;
    L = 12.7*1e-6;
     Sp = I / 4.5 * 1e-6;
    
    e1(i) = uo * L * I^2 / A(i) / B^2;

    eq1 = eg - uo * N * I / B == -(lc / ur);
    eq2 = uo * ur * A(i) * N^2 / (lc + ur * eg) == L;
    sol = solve([eq1, eq2], [eg, N]);
    valueOfeg(i) = double(sol.eg);
    valueOfN(i) = double(sol.N);
end

Aw = valueOfN .* (I / 4.5 * 1e-6) / 0.78;
dp = sqrt(4 * Sp / pi);

Le = pi * (sqrt(4 * Aw ./ (0.78 * pi)) + (Dc / 2));
lm = valueOfN .* pi .* Dc;
DCR = rau .* lm ./ Sp;
e = e1;

L1 = uo * ur .* valueOfN.^2 .* A ./ (Le + ur .* valueOfeg);
L2 = uo * ur .* 12.^2 .* A ./ (Le + ur .* valueOfeg);%plotting
Ee = 1e6 * 0.5 * L1 * I^2;
B1 = L1 .* I ./ (valueOfN .* A);
Em = 1e6 * 0.5 .* B1.^2 .* A .* valueOfeg / uo;
Vc = Le .* A * 1e6;
Bac = B1 * 0.2;
Pin = 1000 * k * f^alpha .* Bac.^Beta .* Vc * 1e-6;
Pw = DCR * I^2;
P = Pw + Pin;

% === PLOTS ===
figure(1)
plot(valueOfN,valueOfeg*1000,'k--','LineWidth',3)
xlabel('N'); ylabel('l_g (mm)'); grid on; axis square;
set(gca,'fontsize',20); set(gcf,'color','white')

figure(2)
plot(A*1e6,valueOfeg*1000,'k--','LineWidth',3)
xlabel('A (mm^2)'); ylabel('l_g (mm)'); grid on; axis square;
set(gca,'fontsize',20); set(gcf,'color','white')

figure(3)
plot(1./valueOfN,valueOfeg*1000,'k--','LineWidth',2)
xlabel('1/N'); ylabel('l_g (mm)'); grid on; axis square;
set(gca,'fontsize',20); set(gcf,'color','white')

figure(4)
plot(A*1e6,valueOfN,'k--','LineWidth',3)
xlabel('A (mm^2)'); ylabel('N'); grid on; axis square;
set(gca,'fontsize',20); set(gcf,'color','white');

% === Normalize values ===
P_norm = (P - min(P)) / (max(P) - min(P));
Vc_norm = (Vc - min(Vc)) / (max(Vc) - min(Vc));

% === Find intersection point ===
[~, idx_maxVc] = max(Vc);
[~, idx_minVc] = min(Vc);
x_intersect = P_norm(idx_maxVc);
y_intersect = Vc_norm(idx_minVc);
intersection_point = [x_intersect, y_intersect];

% === Compute closest point ===
curve_points = [P_norm(:), Vc_norm(:)];
diffs = curve_points - intersection_point;
distances = sqrt(sum(diffs.^2, 2));
[~, idx_closest] = min(distances);

% === Extract closest values ===
closest_P = P(idx_closest);
closest_Vc = Vc(idx_closest);
closest_N = valueOfN(idx_closest);

% === Plot P vs Vc and closest point only ===
figure(5)
plot(P, Vc, 'k--', 'LineWidth', 3); hold on
xline(P(idx_maxVc), 'b', 'LineWidth', 2)
yline(Vc(idx_minVc), 'b', 'LineWidth', 2)
plot(closest_P, closest_Vc, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b')
plot([P(idx_maxVc), closest_P], [Vc(idx_minVc), closest_Vc], 'b--', 'LineWidth', 1.5)

text(closest_P, closest_Vc + 0.3, ...
    sprintf('Closest\nP = %.3f W\nVc = %.3f cm^3', closest_P, closest_Vc), ...
    'Color', 'b', 'FontSize', 12, 'HorizontalAlignment', 'left')

xlabel('Total Loss (W)')
ylabel('Core Volume (cm^3)')
% title('Closest Point to Normalized Intersection')
% legend('Loss vs Volume', ...
%        'Vertical from max V_c', ...
%        'Horizontal from min V_c', ...
%        'Closest point', ...
%        'Distance to closest point', ...
%        'Location', 'northeastoutside')
set(gca, 'fontsize', 18)
axis square
grid on
set(gcf, 'color', 'white')

% === Plot Turns vs Core Volume, marking optimum point ===
figure(6)
plot(valueOfN,Vc, 'r-', 'LineWidth', 2); hold on
plot(closest_N,closest_Vc,  'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
text(closest_N + 0.5,closest_Vc,  ...
    sprintf('Optimal\nN = %.2f', closest_N), 'Color', 'r', 'FontSize', 12)
xlabel('Number of Turns (N)')
ylabel('Core Volume (cm^3)')
title('Turns vs Core Volume with Optimal Point')
grid on
axis square
set(gca,'fontsize', 18)
set(gcf,'color','white')
