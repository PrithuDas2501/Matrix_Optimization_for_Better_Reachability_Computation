%% Input/State Matrix Optimization for Reachable-Set Warping (CORA viz)
%  ADMIRE jet model: 3 states, 4 inputs.
%
%      A = Ac + sum_i alpha_i * GA(:,:,i),   alpha_i in [-1, 1]
%      B = Bc + sum_j beta_j  * GB(:,:,j),   beta_j  in [-1, 1]

clear; clc; %close all;

%% ---------- Problem data: ADMIRE ----------
Ac = [-0.9967   0        0.6176;
       0       -0.5057   0;
      -0.0939   0       -0.2127];

Bc = [ 0       -4.2423   4.2423   1.4871;
       1.6532  -1.2735  -1.2735   0.0024;
       0       -0.2805   0.2805  -0.8820];

Ubound = 0.1 * [-1  1;
                -1  1;
                -1  1;
                -1  1];

n = size(Ac, 1);     % state dim  = 3
m = size(Bc, 2);     % input dim  = 4

% --- Build matrix-zonotope generators (10% relative on each non-zero entry) ---
unc_rel = 0.10;

GA_list = {};
for i = 1:n
    for j = 1:n
        if Ac(i,j) ~= 0
            G = zeros(n,n);
            G(i,j) = unc_rel * abs(Ac(i,j));
            GA_list{end+1} = G; %#ok<SAGROW>
        end
    end
end
GA = cat(3, GA_list{:});

GB_list = {};
for i = 1:n
    for j = 1:m
        if Bc(i,j) ~= 0
            G = zeros(n,m);
            G(i,j) = unc_rel * abs(Bc(i,j));
            GB_list{end+1} = G; %#ok<SAGROW>
        end
    end
end
GB = cat(3, GB_list{:});

A = matZonotope(Ac, GA);
B = matZonotope(Bc, GB);

fprintf('Matrix zonotopes: %d generators in A, %d in B.\n', ...
        size(GA,3), size(GB,3));

% --- Initial set X0: small box around the origin (perturbation from trim) ---
X0_half = 0.1 * ones(n,1);             % half-widths per component
X0_cen  = zeros(n,1);                  % center at origin
V_X0    = box_vertices(X0_cen, X0_half);

% --- Direction of interest (3D) ---
d  = [0; 0; 1];                        % maximize extent along x_1
d  = d / norm(d);

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
fprintf('x_0* = ');               fprintf('%g  ', x0_star); fprintf('\n');
fprintf('Stage 1 best obj = %g\n\n', J_best);

%% ---------- Stage 2: B* and u_0* ----------
nB = size(GB, 3);
U_vertices = box_vertices(zeros(m,1), 0.5*(Ubound(:,2)-Ubound(:,1))) + ...
             0.5*(Ubound(:,1)+Ubound(:,2));      % box vertices, generic m
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
fprintf('u_0* = ');               fprintf('%g  ', u0_star); fprintf('\n');
fprintf('Stage 2 best obj = %g\n\n', best_J);

%% ---------- CORA reachability ----------
sys_nom_unc = linParamSys(A, B, 'constParam');
sys_nom_det = linearSys('nominal_det', Ac, Bc);
sys_opt     = linearSys('optimized',   A_star, B_star);

% Bounding zonotope of X0 (exact for an axis-aligned box)
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

options.timeStep      = 0.05;
options.taylorTerms   = 4;
options.zonotopeOrder = 100;
options.linAlg        = 'standard';

options_param = options;
options_param.intermediateTerms = 4;
options_param = rmfield(options_param, 'linAlg');

tic; R_nom_unc = reach(sys_nom_unc, params, options_param); t_unc = toc;
tic; R_nom_det = reach(sys_nom_det, params, options);       t_det = toc;
tic; R_opt     = reach(sys_opt,     params, options);       t_opt = toc;

fprintf('\nReach computation times:\n');
fprintf('  (A,B) uncertain   : %.3f s\n', t_unc);
fprintf('  (Ac,Bc) nominal   : %.3f s\n', t_det);
fprintf('  (A*,B*) optimized : %.3f s\n', t_opt);

%% ---------- Plot: three 2D state-space projections ----------
figure(1); clf;
projs = {[1 2], [1 3], [2 3]};
for k = 1:numel(projs)
    pdims = projs{k};
    subplot(1,3,k); hold on; grid on; axis equal; box on;

    % plot(R_nom_unc, pdims, 'FaceColor', [1.0 0.6 0.6], 'FaceAlpha', 0.25, ...
    %      'EdgeColor', 'r', 'DisplayName', '(A,B) uncertain');
    plot(R_nom_det, pdims, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.75, ...
         'EdgeColor', [0.3 0.3 0.3], 'DisplayName', '(A_c,B_c)');
    plot(R_opt,     pdims, 'FaceColor', [0.6 1.0 0.6], 'FaceAlpha', 0.40, ...
         'EdgeColor', 'g', 'DisplayName', '(A^*,B^*)');
    plot(X0_zono,   pdims, 'FaceColor', [0.6 0.6 1.0], 'FaceAlpha', 0.6, ...
         'EdgeColor', 'b', 'LineWidth', 1.2, 'DisplayName', 'X_0');

    xlabel(sprintf('x_%d', pdims(1)));
    ylabel(sprintf('x_%d', pdims(2)));
    title(sprintf('R(t; X_0) projection onto x_%d-x_%d', pdims(1), pdims(2)));
    if k == 1, legend('Location','best'); end
end
sgtitle('ADMIRE: reachable set projections');

%% ---------- Plot: state components over time ----------
figure(2); clf;
for i = 1:n
    subplot(n, 1, i); hold on; box on; grid on;
    % plotOverTime(R_nom_unc, i, 'FaceColor', [1.0 0.6 0.6], 'FaceAlpha', 0.4, ...
    %              'EdgeColor', 'r', 'DisplayName', '(A,B) uncertain');
    plotOverTime(R_nom_det, i, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.4, ...
                 'EdgeColor', [0.3 0.3 0.3], 'DisplayName', '(A_c,B_c)');
    plotOverTime(R_opt, i, 'FaceColor', [0.6 1.0 0.6], 'FaceAlpha', 0.4, ...
                 'EdgeColor', 'g', 'DisplayName', '(A^*,B^*)');
    xlabel('t'); ylabel(sprintf('x_%d', i));
    if i == 1
        title('ADMIRE: reachable set components over time');
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
% fill([t_arr; flipud(t_arr)], [upper_unc; flipud(lower_unc)], ...
%      [1.0 0.6 0.6], 'FaceAlpha', 0.25, 'EdgeColor', 'r', ...
%      'DisplayName', '(A,B) uncertain');
fill([t_arr; flipud(t_arr)], [upper_det; flipud(lower_det)], ...
     [0.7 0.7 0.7], 'FaceAlpha', 0.30, 'EdgeColor', [0.3 0.3 0.3], ...
     'DisplayName', '(A_c,B_c)');
fill([t_arr; flipud(t_arr)], [upper_opt; flipud(lower_opt)], ...
     [0.6 1.0 0.6], 'FaceAlpha', 0.35, 'EdgeColor', 'g', ...
     'DisplayName', '(A^*,B^*)');
xlabel('t'); ylabel('d^\top x');
title(sprintf('ADMIRE: reachable extent along d = [%g, %g, %g]^\\top', d(1), d(2), d(3)));
legend('Location', 'best');

% --- Numerical summary at t = T ---
fprintf('\nExtent along d at t=T:\n');
fprintf('  uncertain (A,B)  : upper = %+.4f,  lower = %+.4f,  width = %.4f\n', ...
        upper_unc(end), lower_unc(end), upper_unc(end)-lower_unc(end));
fprintf('  nominal (Ac,Bc)  : upper = %+.4f,  lower = %+.4f,  width = %.4f\n', ...
        upper_det(end), lower_det(end), upper_det(end)-lower_det(end));
fprintf('  optimized (A*,B*): upper = %+.4f,  lower = %+.4f,  width = %.4f\n', ...
        upper_opt(end), lower_opt(end), upper_opt(end)-lower_opt(end));
fprintf('  optimization gain over nominal (upper): %+.4f\n', ...
        upper_opt(end) - upper_det(end));

% %% ---------- Plot: red and green reachable sets in 3D ----------
% figure(4); clf;
% hold on; grid on; box on;
% 
% view(3);
% axis equal;
% axis vis3d;
% 
% xlabel('x_1');
% ylabel('x_2');
% zlabel('x_3');
% title('ADMIRE: 3D Reachable Sets');
% 
% % Number of reachable-set snapshots to display.
% % Intermediate snapshots show the evolution of the sets over time.
% nSnap = 12;
% 
% N_unc = numel(R_nom_unc.timePoint.set);
% N_opt = numel(R_opt.timePoint.set);
% 
% % Use the number of time points common to both results
% N_common = min(N_unc, N_opt);
% snap_idx = unique(round(linspace(1, N_common, nSnap)));
% 
%  %% Faint intermediate reachable sets
% % for k = snap_idx(1:end-1)
% % 
% %     % Red: uncertain system
% %     plot(R_nom_unc.timePoint.set{k}, [1 2 3], ...
% %          'FaceColor', [1.0 0.6 0.6], ...
% %          'FaceAlpha', 0.06, ...
% %          'EdgeColor', 'r', ...
% %          'LineWidth', 0.25, ...
% %          'HandleVisibility', 'off');
% % 
% %     % Green: optimized system
% %     plot(R_opt.timePoint.set{k}, [1 2 3], ...
% %          'FaceColor', [0.6 1.0 0.6], ...
% %          'FaceAlpha', 0.08, ...
% %          'EdgeColor', 'g', ...
% %          'LineWidth', 0.25, ...
% %          'HandleVisibility', 'off');
% % end
% 
% %% Final reachable sets at t = T
% plot(R_nom_unc.timePoint.set{end}, [1 2 3], ...
%      'FaceColor', [1.0 0.35 0.35], ...
%      'FaceAlpha', 0.30, ...
%      'EdgeColor', 'r', ...
%      'LineWidth', 1.5, ...
%      'DisplayName', '(A,B) uncertain at t=T');
% 
% plot(R_opt.timePoint.set{end}, [1 2 3], ...
%      'FaceColor', [0.35 0.90 0.35], ...
%      'FaceAlpha', 0.40, ...
%      'EdgeColor', 'g', ...
%      'LineWidth', 1.5, ...
%      'DisplayName', '(A^*,B^*) optimized at t=T');
% 
% %% Initial set
% plot(X0_zono, [1 2 3], ...
%      'FaceColor', [0.6 0.6 1.0], ...
%      'FaceAlpha', 0.65, ...
%      'EdgeColor', 'b', ...
%      'LineWidth', 1.2, ...
%      'DisplayName', 'X_0');
% 
% %% Direction of interest
% % Scale the arrow according to the final reachable-set extent
% arrow_len = 0.25 * max([ ...
%     abs(upper_unc(end)), ...
%     abs(lower_unc(end)), ...
%     abs(upper_opt(end)), ...
%     abs(lower_opt(end))]);
% 
% % Prevent a zero-length arrow in a degenerate case
% if arrow_len < 1e-6
%     arrow_len = 1;
% end
% 
% quiver3(0, 0, 0, ...
%         arrow_len*d(1), ...
%         arrow_len*d(2), ...
%         arrow_len*d(3), ...
%         0, ...
%         'Color', 'k', ...
%         'LineWidth', 2, ...
%         'MaxHeadSize', 0.5, ...
%         'DisplayName', 'Direction d');
% 
% legend('Location', 'bestoutside');
% 
% % Choose a useful initial camera angle
% view(135, 25);
% 
% % Improve the appearance of transparent surfaces
% camlight headlight;
% lighting gouraud;

% %% ---------- Plot: 3D reachable sets ----------
% figure(4); clf; hold on; grid on; box on;
% view(3); axis equal; axis vis3d;
% xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
% title('ADMIRE: 3D reachable sets');
% 
% % Sparse snapshot times so the plot stays readable
% nSnap    = 12;
% snap_idx = round(linspace(1, T, nSnap));
% 
% % Faint snapshots for t < T
% for s = snap_idx(1:end-1)
%     plot(R_nom_unc.timePoint.set{s}, [1 2 3], ...
%          'FaceColor', [1.0 0.6 0.6], 'FaceAlpha', 0.08, ...
%          'EdgeColor', 'r', 'LineStyle', '-', 'LineWidth', 0.3, ...
%          'HandleVisibility', 'off');
%     plot(R_nom_det.timePoint.set{s}, [1 2 3], ...
%          'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.18, ...
%          'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.3, ...
%          'HandleVisibility', 'off');
%     plot(R_opt.timePoint.set{s}, [1 2 3], ...
%          'FaceColor', [0.6 1.0 0.6], 'FaceAlpha', 0.13, ...
%          'EdgeColor', 'g', 'LineWidth', 0.3, ...
%          'HandleVisibility', 'off');
% end
% 
% % Bold t = T snapshot — these get the legend entries
% plot(R_nom_unc.timePoint.set{end}, [1 2 3], ...
%      'FaceColor', [1.0 0.4 0.4], 'FaceAlpha', 0.25, ...
%      'EdgeColor', 'r', 'LineWidth', 1.0, ...
%      'DisplayName', '(A,B) uncertain @ t=T');
% plot(R_nom_det.timePoint.set{end}, [1 2 3], ...
%      'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.55, ...
%      'EdgeColor', 'k', 'LineWidth', 1.0, ...
%      'DisplayName', '(A_c,B_c) @ t=T');
% plot(R_opt.timePoint.set{end}, [1 2 3], ...
%      'FaceColor', [0.3 0.8 0.3], 'FaceAlpha', 0.40, ...
%      'EdgeColor', 'g', 'LineWidth', 1.0, ...
%      'DisplayName', '(A^*,B^*) @ t=T');
% 
% % Initial set
% plot(X0_zono, [1 2 3], 'FaceColor', [0.6 0.6 1.0], 'FaceAlpha', 0.7, ...
%      'EdgeColor', 'b', 'LineWidth', 1.0, 'DisplayName', 'X_0');
% 
% % Direction of interest, scaled to a visible length
% arrow_len = 1.2 * max(abs([upper_unc; lower_unc]));
% quiver3(0, 0, 0, arrow_len*d(1), arrow_len*d(2), arrow_len*d(3), 0, ...
%         'Color', 'k', 'LineWidth', 2, 'MaxHeadSize', 0.3, ...
%         'DisplayName', 'd');
% 
% legend('Location', 'bestoutside');
% view(135, 25);
% camlight headlight; lighting gouraud;

%% ================== local helpers ==================
function V = box_vertices(c, r)
    % Enumerate vertices of an axis-aligned box {c + diag(r)*s : s in {-1,1}^d}.
    d  = numel(c);
    nv = 2^d;
    V  = zeros(d, nv);
    for k = 0:nv-1
        bits = bitget(k, 1:d);          % 0/1 vector, length d
        s    = 2*bits' - 1;              % -1/+1 vector
        V(:, k+1) = c + r .* s;
    end
end