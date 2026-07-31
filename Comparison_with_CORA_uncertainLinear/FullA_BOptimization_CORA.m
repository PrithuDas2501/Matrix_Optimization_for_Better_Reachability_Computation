%% Input/State Matrix Optimization for Reachable-Set Warping (CORA viz)
%  Same two-stage algorithm as before:
%      A = Ac + sum_i alpha_i * GA(:,:,i),   alpha_i in [-1, 1]
%      B = Bc + sum_j beta_j  * GB(:,:,j),   beta_j  in [-1, 1]
%  Stage 1: find A* and x_0* by fixed-point iteration on P_0(A).
%  Stage 2: find B* and u_0* given A*.
%
%  Visualization swapped from boundary-trajectory sweep to CORA's
%  reach()/plot() pipeline — once A* and B* are fixed, the dynamics
%  for visualization purposes are linear-deterministic.

clear; clc; % close all;

%% ---------- Problem data ----------
Ac = [ 0   1;
      -2  -0.8];
GA = zeros(2,2,2);
GA(:,:,1) = [0.05  0;     0    0   ];
GA(:,:,2) = [0     0;     0.1  0   ];

A = matZonotope(Ac, GA);

Bc = [ 0   1;
       1   0];
GB = zeros(2,2,2);
GB(:,:,1) = [0.1  0;      0    0   ];
GB(:,:,2) = [0    0;      0    0.1 ];

B = matZonotope(Bc, GB);

V_X0 = [ 0.9  0.9  1.1  1.1;
         0.9  1.1  0.9  1.1 ];

Ubound = 0.1 * [-1  1;
                -1  1];

theta_d = 0*pi/2;
d  = [cos(theta_d); sin(theta_d)];
tf = 5;

max_iter = 20;
tol      = 1e-8;

%% ---------- Stage 1: A* and x_0* ----------
nA = size(GA, 3);
nX = size(V_X0, 2);

J_best  = -inf;
A_best  = Ac;
x0_best = V_X0(:,1);

A_curr = Ac;
P0     = expm(A_curr.' * tf) * d;

for iter = 1:max_iter
    best_J_iter = -inf;
    A_iter      = Ac;
    x0_iter     = V_X0(:,1);

    for v = 1:nX
        xv = V_X0(:, v);
        coeffs = zeros(nA, 1);
        for i = 1:nA
            coeffs(i) = P0.' * GA(:,:,i) * xv;
        end
        alpha = sign(coeffs);  alpha(alpha == 0) = 1;

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

    J_true = d.' * expm(A_iter * tf) * A_iter * x0_iter;
    if J_true > J_best
        J_best  = J_true;
        A_best  = A_iter;
        x0_best = x0_iter;
    end

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

%% ---------- Stage 2: B* and u_0* ----------
nB = size(GB, 3);
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

%% ---------- CORA reachability ----------
% Two deterministic LTI systems: nominal and optimized
% sys_nom = linearSys('nominal',   Ac,     Bc);
% sys_nom = linearSys('nominal', A, B);          % <-- doesn't accept matZonotope
sys_nom_unc = linParamSys(A, B, 'constParam');             % uncertain (matZono)
sys_nom_det = linearSys('nominal_det', Ac, Bc);            % nominal deterministic
sys_opt     = linearSys('optimized',   A_star, B_star);    % optimized deterministic

% Bounding zonotope of X0 (exact when X0 is an axis-aligned box,
% over-approximation otherwise — sufficient for visualization).
c_X0 = (max(V_X0,[],2) + min(V_X0,[],2)) / 2;
r_X0 = (max(V_X0,[],2) - min(V_X0,[],2)) / 2;
X0_zono = zonotope(c_X0, diag(r_X0));

% Box U as a zonotope
Uc = (Ubound(:,1) + Ubound(:,2)) / 2;
Ur = (Ubound(:,2) - Ubound(:,1)) / 2;
U_zono = zonotope(Uc, diag(Ur));

params.tFinal = tf;
params.R0     = X0_zono;
params.U      = U_zono;

% --- options for the deterministic optimized system ---
options.timeStep      = 0.05;
options.taylorTerms   = 4;
options.zonotopeOrder = 100;
options.linAlg        = 'standard';

% --- options for the uncertain nominal system ---
options_param = options;
options_param.intermediateTerms = 4;
options_param = rmfield(options_param, 'linAlg');   % linParamSys doesn't take linAlg


tic; R_nom_unc = reach(sys_nom_unc, params, options_param); t_unc = toc;
tic; R_nom_det = reach(sys_nom_det, params, options);       t_det = toc;
tic; R_opt     = reach(sys_opt,     params, options);       t_opt = toc;

fprintf('\nReach computation times:\n');
fprintf('  (A,B) uncertain   : %.3f s\n', t_unc);
fprintf('  (Ac,Bc) nominal   : %.3f s\n', t_det);
fprintf('  (A*,B*) optimized : %.3f s\n', t_opt);

%% ---------- Plot: state-space tube ----------
figure(1); clf; hold on; grid on; axis equal;
xlabel('x_1'); ylabel('x_2');
title('Reachable set R(t; X_0): three systems');

plot(R_nom_unc, [1 2], 'FaceColor', [1.0 0.6 0.6], 'FaceAlpha', 0.25, ...
     'EdgeColor', 'r', 'DisplayName', 'R_{(A,B) uncertain}');
plot(R_nom_det, [1 2], 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.60, ...
     'EdgeColor', [0.3 0.3 0.3], 'DisplayName', 'R_{(A_c,B_c)}');
plot(R_opt,     [1 2], 'FaceColor', [0.6 1.0 0.6], 'FaceAlpha', 0.40, ...
     'EdgeColor', 'g', 'DisplayName', 'R_{(A^*,B^*)}');

plot(X0_zono, [1 2], 'FaceColor', [0.6 0.6 1.0], 'FaceAlpha', 0.6, ...
     'EdgeColor', 'b', 'LineWidth', 1.5, 'DisplayName', 'X_0');
quiver(0, 0, d(1), d(2), 1.5, 'k', 'LineWidth', 2, ...
       'MaxHeadSize', 0.5, 'DisplayName', 'd');
legend('Location', 'best');

%% ---------- Plot: state components over time ----------
figure(2); clf;
n = size(Ac, 1);
for i = 1:n
    subplot(n, 1, i); hold on; box on; grid on;
    plotOverTime(R_nom_unc, i, 'FaceColor', [1.0 0.6 0.6], ...
                 'FaceAlpha', 0.4, 'EdgeColor', 'r', ...
                 'DisplayName', '(A,B) uncertain');
    plotOverTime(R_nom_det, i, 'FaceColor', [0.7 0.7 0.7], ...
                 'FaceAlpha', 0.4, 'EdgeColor', [0.3 0.3 0.3], ...
                 'DisplayName', '(A_c,B_c)');
    plotOverTime(R_opt, i, 'FaceColor', [0.6 1.0 0.6], ...
                 'FaceAlpha', 0.4, 'EdgeColor', 'g', ...
                 'DisplayName', '(A^*,B^*)');
    xlabel('t');  ylabel(sprintf('x_%d', i));
    if i == 1
        title('Reachable set components over time');
        legend('Location', 'best');
    end
end

%% ---------- Plot: extent along direction d over time ----------
t_grid = R_nom_unc.timePoint.time;
T = numel(t_grid);

upper_unc = zeros(T,1);   lower_unc = zeros(T,1);
upper_det = zeros(T,1);   lower_det = zeros(T,1);
upper_opt = zeros(T,1);   lower_opt = zeros(T,1);

for k = 1:T
    upper_unc(k) =  supportFunc(R_nom_unc.timePoint.set{k},  d, 'upper');
    lower_unc(k) = -supportFunc(R_nom_unc.timePoint.set{k}, -d, 'upper');
    upper_det(k) =  supportFunc(R_nom_det.timePoint.set{k},  d, 'upper');
    lower_det(k) = -supportFunc(R_nom_det.timePoint.set{k}, -d, 'upper');
    upper_opt(k) =  supportFunc(R_opt.timePoint.set{k},      d, 'upper');
    lower_opt(k) = -supportFunc(R_opt.timePoint.set{k},     -d, 'upper');
end
t_arr = cell2mat(t_grid(:));

figure(3); clf; hold on; grid on; box on;
fill([t_arr; flipud(t_arr)], [upper_unc; flipud(lower_unc)], ...
     [1.0 0.6 0.6], 'FaceAlpha', 0.25, 'EdgeColor', 'r', ...
     'DisplayName', '(A,B) uncertain');
fill([t_arr; flipud(t_arr)], [upper_det; flipud(lower_det)], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.30, 'EdgeColor', [0.3 0.3 0.3], ...
     'DisplayName', '(A_c,B_c)');
fill([t_arr; flipud(t_arr)], [upper_opt; flipud(lower_opt)], ...
     [0.6 1.0 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'g', ...
     'DisplayName', '(A^*,B^*)');
xlabel('t'); ylabel('d^\top x');
title('Reachable extent along direction d');
legend('Location', 'best');

% --- Numerical summary at t = T ---
fprintf('\nExtent along d at t=T:\n');
fprintf('  uncertain (A,B):  upper = %+.4f,  lower = %+.4f,  width = %.4f\n', ...
        upper_unc(end), lower_unc(end), upper_unc(end)-lower_unc(end));
fprintf('  nominal (Ac,Bc):  upper = %+.4f,  lower = %+.4f,  width = %.4f\n', ...
        upper_det(end), lower_det(end), upper_det(end)-lower_det(end));
fprintf('  optimized (A*,B*): upper = %+.4f,  lower = %+.4f,  width = %.4f\n', ...
        upper_opt(end), lower_opt(end), upper_opt(end)-lower_opt(end));