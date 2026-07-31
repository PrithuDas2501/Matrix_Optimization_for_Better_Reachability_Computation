%% Input Matrix Optimization for Reachable Set Warping -- ADMIRE fighter jet
%  Implements Algorithm 1 from Das & Ornik for the linearized ADMIRE model
%  (Section V-A).
%
%  States   x = [p; q; r]     (roll, pitch, yaw rates, rad/s)
%  Controls u in R^4          (canards, left/right elevons, rudder)


clear; close all; clc;

%% ---------- Problem data ----------
A  = [-0.9967   0        0.6176;
       0       -0.5057   0;
      -0.0939   0       -0.2127];

B0 = [ 0       -4.2423   4.2423   1.4871;
       1.6532  -1.2735  -1.2735   0.0024;
       0       -0.2805   0.2805  -0.8820];

Ubound = 0.1 * [-1 1;                 % +/- 0.1 rad (~5.7 deg) per surface
                -1 1;
                -1 1;
                -1 1];

X0     = [0; 0; 0];
tf     = 2;                           % horizon (paper uses T = 2 s)
gamma  = 0.5;                         % squared Frobenius radius:  ||B - B0||_F^2 <= gamma
d      = [1; 0; 0];                   % direction of interest (roll rate p)

%% ---------- Algorithm 1 ----------
% Step 2: P0 = exp(A' * T) * d
P0 = expm(A' * tf) * d;

% Step 3: enumerate the 2^m vertices of the box U
m = size(Ubound, 1);
N = 2^m;
masks = dec2bin(0:N-1, m) - '0';                        % N-by-m of 0/1
U_vertices = Ubound(:,1) .* (1 - masks).' + Ubound(:,2) .* masks.';   % m-by-N

% Steps 4-6: closed-form solution to
%     max_{B : ||B - B0||_F^2 <= gamma}  P0' * B * u_i
%     =>  B_i = B0 + sqrt(gamma) * (P0 * u_i') / ||P0 * u_i'||_F
B_candidates = cell(N, 1);
obj_vals     = zeros(N, 1);
for i = 1:N
    ui              = U_vertices(:, i);
    G               = P0 * ui.';
    B_candidates{i} = B0 + sqrt(gamma) * G / norm(G, 'fro');
    obj_vals(i)     = P0.' * B_candidates{i} * ui;
end

% Steps 7-8: pick the winning (B_i, u_i) pair
[~, i_star] = max(obj_vals);       % swap to min for shrinkage along d
B_star      = B_candidates{i_star};

fprintf('Winning vertex index: %d\n', i_star);
fprintf('u* = [%s]\n', sprintf('%.2f ', U_vertices(:, i_star)));
fprintf('||B* - B0||_F = %.4f   (bound = %.4f)\n', ...
        norm(B_star - B0, 'fro'), sqrt(gamma));
fprintf('B* =\n');  disp(B_star);

%% ---------- Reachable-set visualization (3-D) ----------
% Sweep P(T) over the unit sphere via a spherical-coord grid.
theta_list = linspace(0, pi,   14);       % polar
phi_list   = linspace(0, 2*pi, 28);       % azimuth
[TH, PH]   = meshgrid(theta_list, phi_list);
Dirs = [ sin(TH(:)).*cos(PH(:)) , ...
         sin(TH(:)).*sin(PH(:)) , ...
         cos(TH(:)) ].';                  % 3-by-M unit vectors
M = size(Dirs, 2);

vertices_opt = zeros(M, 3);
vertices_nom = zeros(M, 3);
G_opt = 0;  G_nom = 0;
tspan    = [0 tf];
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

fprintf('\nSweeping %d boundary directions...\n', M);
for c = 1:M
    P_T = Dirs(:, c);
    y0  = [X0; expm(A' * tf) * P_T];      % P(0) from P(T)

    [~, y] = ode45(@(t,y) boundary_traj(y, A, B_star, Ubound), tspan, y0, ode_opts);
    vertices_opt(c, :) = y(end, 1:3);
    G_opt = G_opt + P_T.' * (y(end, 1:3).' - X0);

    [~, y] = ode45(@(t,y) boundary_traj(y, A, B0,     Ubound), tspan, y0, ode_opts);
    vertices_nom(c, :) = y(end, 1:3);
    G_nom = G_nom + P_T.' * (y(end, 1:3).' - X0);
end

fprintf('Sum_d  d''*(X_d - X0)  (nominal)   = %.4f\n', G_nom);
fprintf('Sum_d  d''*(X_d - X0)  (optimized) = %.4f\n', G_opt);

figure; hold on; grid on; axis equal;
xlabel('p (rad/s)');  ylabel('q (rad/s)');  zlabel('r (rad/s)');
title('Reachable set: nominal (red) vs. optimized (green)');

K = convhull(vertices_nom(:,1), vertices_nom(:,2), vertices_nom(:,3));
patch('Faces', K, 'Vertices', vertices_nom, ...
      'FaceColor', 'r', 'FaceAlpha', 0.40, 'EdgeColor', 'none');

K = convhull(vertices_opt(:,1), vertices_opt(:,2), vertices_opt(:,3));
patch('Faces', K, 'Vertices', vertices_opt, ...
      'FaceColor', 'g', 'FaceAlpha', 0.40, 'EdgeColor', 'none');

% Arrow showing d
ax = 0.7 * max(abs([vertices_opt(:); vertices_nom(:)]));
quiver3(0, 0, 0, ax*d(1), ax*d(2), ax*d(3), ...
        0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.5);

view(40, 25);  camlight;  lighting gouraud;

%% ================== local functions ==================
function ydot = boundary_traj(y, A, B, Ubound)
    % State/costate dynamics with Hamiltonian-maximizing control.
    % Box U  =>  u_i* is bang-bang on sign of (B'P)_i.
    n = numel(y) / 2;
    X = y(1:n);
    P = y(n+1:end);

    pos = (B.' * P) > 0;                         % componentwise
    u   = pos .* Ubound(:, 2) + (~pos) .* Ubound(:, 1);

    ydot = [ A*X + B*u;
            -A.'*P ];
end