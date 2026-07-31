%% Input/State Matrix Optimization for Reachable-Set Warping
%  Two-stage extension of Algorithm 1 (Das & Ornik) to the case where
%  BOTH A and B are uncertain matrix zonotopes:
%
%      A = Ac + sum_{i=1..nA} alpha_i * GA(:,:,i),   alpha_i in [-1, 1]
%      B = Bc + sum_{j=1..nB} beta_j  * GB(:,:,j),   beta_j  in [-1, 1]
%
%  X0 is a polytope (given by its vertices V_X0), U is a box.
%  Goal: pick (A*, x_0*) and (B*, u_0*) that maximize the reachable-set
%  extent at time T along a given direction d.
%
%  Notation:
%      P_0(A) = expm(A' * T) * d           (terminal-cost back-propagation)
%      max_{A, x_0}   d' * expm(A*T) * A * x_0   = max_{A, x_0} P_0(A)' * A * x_0
%      max_{B, u_0}   d' * expm(A*T) * B * u_0   = max_{B, u_0} P_0(A)' * B * u_0
%  Stage 1 finds A* (with P_0 itself depending on A, handled by fixed-point
%  iteration). Stage 2 then finds B* using the final P_0.

clear; clc; % close all;

%% ---------- Problem data ----------
% --- State-matrix zonotope ---
Ac = [ 0   1;
      -2  -0.8];
GA = zeros(2,2,2);
GA(:,:,1) = [0.05  0;     0    0   ];   % uncertainty on (1,1)
GA(:,:,2) = [0     0;     0.1  0   ];   % uncertainty on (2,1)

% --- Input-matrix zonotope ---
Bc = 1*[ 0   1;
       1   0];
GB = zeros(2,2,2);
GB(:,:,1) = [0.1  0;      0    0   ];   % uncertainty on (1,1)
GB(:,:,2) = [0    0;      0    0.1 ];   % uncertainty on (2,2)

% --- Initial set X0 (polytope, given by its vertices) ---
V_X0 = [ 0.9  0.9  1.1  1.1;
         0.9  1.1  0.9  1.1 ];          % small box around [1;1]
% (For a CORA polytope object P, use V_X0 = vertices(P);)

% --- Control box U ---
Ubound = 0.1 * [-1  1;                   % per-channel control bounds
                -1  1];

% --- Direction of interest, horizon ---
theta_d = 0*pi/2;
d  = [cos(theta_d); sin(theta_d)];
tf = 5;

% --- Iteration settings for Stage 1 ---
max_iter = 20;
tol      = 1e-8;

%% ---------- Stage 1: Find A* and x_0* ----------
nA = size(GA, 3);
nX = size(V_X0, 2);

J_best  = -inf;          % best TRUE objective d' * expm(A*T) * A * x_0
A_best  = Ac;
x0_best = V_X0(:,1);

A_curr = Ac;
P0     = expm(A_curr.' * tf) * d;

for iter = 1:max_iter

    % ----- inner: maximize  P0' * A * x_0  over (A, x_0)  with P0 fixed -----
    best_J_iter = -inf;
    A_iter      = Ac;
    x0_iter     = V_X0(:,1);

    for v = 1:nX
        xv = V_X0(:, v);

        % Linear-in-alpha coefficients:  P0' * G_i * xv
        coeffs = zeros(nA, 1);
        for i = 1:nA
            coeffs(i) = P0.' * GA(:,:,i) * xv;
        end
        alpha = sign(coeffs);   alpha(alpha == 0) = 1;

        A_try = Ac;
        for i = 1:nA
            A_try = A_try + alpha(i) * GA(:,:,i);
        end

        J_lin = P0.' * A_try * xv;
        if J_lin > best_J_iter
            best_J_iter = J_lin;
            A_iter      = A_try;
            x0_iter     = xv;
        end
    end

    % ----- evaluate the TRUE non-linear objective at this iterate -----
    J_true = d.' * expm(A_iter * tf) * A_iter * x0_iter;
    if J_true > J_best
        J_best  = J_true;
        A_best  = A_iter;
        x0_best = x0_iter;
    end

    % ----- update P0 from the new A; check fixed-point convergence -----
    P0_new = expm(A_iter.' * tf) * d;
    if norm(P0_new - P0) < tol
        fprintf('Stage 1 converged at iteration %d.\n', iter);
        break
    end
    P0 = P0_new;
end

A_star  = A_best;
x0_star = x0_best;
P0_star = expm(A_star.' * tf) * d;

fprintf('A* =\n');                disp(A_star);
fprintf('x_0* = (%g, %g)\n',      x0_star);
fprintf('Stage 1 best obj = %g\n\n', J_best);

%% ---------- Stage 2: Find B* and u_0* given A* ----------
nB = size(GB, 3);

% Vertices of the box U
U_vertices = [Ubound(1,1) Ubound(1,1) Ubound(1,2) Ubound(1,2);
              Ubound(2,1) Ubound(2,2) Ubound(2,1) Ubound(2,2)];
nU = size(U_vertices, 2);

best_J  = -inf;
B_star  = Bc;
u0_star = U_vertices(:,1);

for v = 1:nU
    uv = U_vertices(:, v);

    coeffs = zeros(nB, 1);
    for j = 1:nB
        coeffs(j) = P0_star.' * GB(:,:,j) * uv;
    end
    beta = sign(coeffs);  beta(beta == 0) = 1;

    B_try = Bc;
    for j = 1:nB
        B_try = B_try + beta(j) * GB(:,:,j);
    end

    J = P0_star.' * B_try * uv;
    if J > best_J
        best_J  = J;
        B_star  = B_try;
        u0_star = uv;
    end
end

fprintf('B* =\n');                disp(B_star);
fprintf('u_0* = (%g, %g)\n',      u0_star);
fprintf('Stage 2 best obj = %g\n\n', best_J);

%% ---------- Reachable-set visualization (nominal vs. optimized) ----------
% Sweep terminal costate direction theta around the unit circle.
% For each theta, pick the X0-vertex that maximizes P(0)' * x_0 and
% propagate forward with the bang-bang Hamiltonian-maximizing control.

theta_list = linspace(0, 2*pi, 80);
M = numel(theta_list);

snap_fracs = linspace(0.05, 1, 20);
t_snaps    = snap_fracs * tf;
S          = numel(t_snaps);

snap_opt = zeros(M, S, 2);
snap_nom = zeros(M, S, 2);

tspan    = [0, t_snaps];
ode_opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

for c = 1:M
    P_T = [cos(theta_list(c)); sin(theta_list(c))];

    % ----- optimized system (A_star, B_star) -----
    P0_opt = expm(A_star.' * tf) * P_T;
    x0_opt = pick_X0_vertex(V_X0, P0_opt);
    [~, y] = ode45(@(t,y) boundary_traj(y, A_star, B_star, Ubound), ...
                   tspan, [x0_opt; P0_opt], ode_opts);
    snap_opt(c, :, :) = y(2:end, 1:2);

    % ----- nominal system (Ac, Bc) -----
    P0_nom = expm(Ac.' * tf) * P_T;
    x0_nom = pick_X0_vertex(V_X0, P0_nom);
    [~, y] = ode45(@(t,y) boundary_traj(y, Ac, Bc, Ubound), ...
                   tspan, [x0_nom; P0_nom], ode_opts);
    snap_nom(c, :, :) = y(2:end, 1:2);
end

%% ---------- Plot ----------
figure(1); clf; hold on; grid on; axis equal;
xlabel('x_1'); ylabel('x_2');
title('Reachable set R(t; X_0): nominal (red) vs. optimized (green)');

alphas_face = linspace(0.06, 0.35, S);
alphas_edge = linspace(0.35, 1.00, S);

for s = 1:S
    Vnom = squeeze(snap_nom(:, s, :));
    Vopt = squeeze(snap_opt(:, s, :));
    Knom = convhull(Vnom(:,1), Vnom(:,2));
    Kopt = convhull(Vopt(:,1), Vopt(:,2));

    if s < S
        fill(Vnom(Knom,1), Vnom(Knom,2), 'r', ...
             'FaceAlpha', alphas_face(s), 'EdgeColor', 'r', ...
             'EdgeAlpha', alphas_edge(s), 'LineStyle','--', ...
             'HandleVisibility','off');
        fill(Vopt(Kopt,1), Vopt(Kopt,2), 'g', ...
             'FaceAlpha', alphas_face(s), 'EdgeColor', 'g', ...
             'EdgeAlpha', alphas_edge(s), 'LineStyle','--', ...
             'HandleVisibility','off');
    else
        fill(Vnom(Knom,1), Vnom(Knom,2), 'r', ...
             'FaceAlpha', alphas_face(s), 'EdgeColor', 'r', ...
             'LineWidth', 1.5, 'DisplayName', 'R_{(A_c,B_c)}(T)');
        fill(Vopt(Kopt,1), Vopt(Kopt,2), 'g', ...
             'FaceAlpha', alphas_face(s), 'EdgeColor', 'g', ...
             'LineWidth', 1.5, 'DisplayName', 'R_{(A^*,B^*)}(T)');
    end
end

% X0 polytope
K_X0 = convhull(V_X0(1,:), V_X0(2,:));
fill(V_X0(1, K_X0), V_X0(2, K_X0), 'b', 'FaceAlpha', 0.4, ...
     'EdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'X_0');

% Direction d
quiver(0, 0, 0.2*d(1), 0.2*d(2), 1.5, 'k', 'LineWidth', 2, ...
       'MaxHeadSize', 0.5, 'DisplayName', 'd');

legend('Location', 'best');

% %% ---------- Plot: state-space tube ----------
% figure(2); clf; hold on; grid on; axis equal;
% xlabel('x_1'); ylabel('x_2');
% title('Reachable set R(t; X_0): nominal (red) vs. optimized (green)');
% 
% % Full time-interval reachable tubes
% plot(R_nom, [1 2], 'FaceColor', [1.0 0.6 0.6], 'FaceAlpha', 0.30, ...
%      'EdgeColor', 'r', 'DisplayName', 'R_{(A_c,B_c)}');
% plot(R_opt, [1 2], 'FaceColor', [0.6 1.0 0.6], 'FaceAlpha', 0.30, ...
%      'EdgeColor', 'g', 'DisplayName', 'R_{(A^*,B^*)}');
% 
% % Initial set X0
% plot(X0_zono, [1 2], 'FaceColor', [0.6 0.6 1.0], 'FaceAlpha', 0.6, ...
%      'EdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'X_0');
% 
% % Direction of interest
% quiver(0, 0, d(1), d(2), 1.5, 'k', 'LineWidth', 2, ...
%        'MaxHeadSize', 0.5, 'DisplayName', 'd');
% 
% legend('Location', 'best');
% 
% %% ---------- Plot: state components over time ----------
% figure(3); clf;
% n = size(Ac, 1);
% for i = 1:n
%     subplot(n, 1, i); hold on; box on; grid on;
%     plotOverTime(R_nom, i, 'FaceColor', [1.0 0.6 0.6], ...
%                  'FaceAlpha', 0.4, 'EdgeColor', 'r', ...
%                  'DisplayName', 'nominal');
%     plotOverTime(R_opt, i, 'FaceColor', [0.6 1.0 0.6], ...
%                  'FaceAlpha', 0.4, 'EdgeColor', 'g', ...
%                  'DisplayName', 'optimized');
%     xlabel('t');
%     ylabel(sprintf('x_%d', i));
%     if i == 1
%         title('Reachable set components over time');
%         legend('Location', 'best');
%     end
% end
% 
% %% ---------- Plot: extent along direction d over time ----------
% % For each time-point reachable set R(t), compute support function
% % h_R(d) = max_{x in R(t)} d' x  and  -h_R(-d) = min_{x in R(t)} d' x.
% % This shows directly how much further the optimized system reaches
% % in the direction d.
% 
% t_grid = R_nom.timePoint.time;          % shared grid for both
% T = numel(t_grid);
% 
% upper_nom = zeros(T,1);   lower_nom = zeros(T,1);
% upper_opt = zeros(T,1);   lower_opt = zeros(T,1);
% 
% for k = 1:T
%     upper_nom(k) =  supportFunc(R_nom.timePoint.set{k},  d, 'upper');
%     lower_nom(k) = -supportFunc(R_nom.timePoint.set{k}, -d, 'upper');
%     upper_opt(k) =  supportFunc(R_opt.timePoint.set{k},  d, 'upper');
%     lower_opt(k) = -supportFunc(R_opt.timePoint.set{k}, -d, 'upper');
% end
% 
% t_arr = cell2mat(t_grid(:));
% 
% figure(4); clf; hold on; grid on; box on;
% fill([t_arr; flipud(t_arr)], [upper_nom; flipud(lower_nom)], ...
%      [1.0 0.6 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'r', ...
%      'DisplayName', 'nominal extent along d');
% fill([t_arr; flipud(t_arr)], [upper_opt; flipud(lower_opt)], ...
%      [0.6 1.0 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'g', ...
%      'DisplayName', 'optimized extent along d');
% xlabel('t'); ylabel('d^\top x');
% title('Reachable extent along direction d');
% legend('Location', 'best');

%% ================== local functions ==================
function ydot = boundary_traj(y, A, B, Ubound)
    % State/costate dynamics with the Hamiltonian-maximizing
    % bang-bang control on a box U.
    n = numel(y) / 2;
    X = y(1:n);
    P = y(n+1:end);

    pos = (B.' * P) > 0;
    u   = pos .* Ubound(:,2) + (~pos) .* Ubound(:,1);

    ydot = [ A*X + B*u;
            -A.'*P ];
end

function xv_best = pick_X0_vertex(V_X0, P0)
    % Return the vertex of X0 maximizing P0' * x_0.
    [~, k] = max(P0.' * V_X0);
    xv_best = V_X0(:, k);
end