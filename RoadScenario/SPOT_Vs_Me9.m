function overtake_frenet_animation(roadType)
% OVERTAKE_FRENET_ANIMATION
% Frenet-frame version of the green-vs-red overtake/wait demo. The road is a
% reference curve parameterised by arclength s; a curvilinear coordinate (s,d)
% maps to world as   world = ref(s) + d * N(s),   N(s) = [-sin th(s); cos th(s)].
% Lanes are constant lateral offsets: left d=+laneWidth, middle d=0, right d=-laneWidth.
%
% The reach set is computed in the FLAT Frenet displacement space (Ds_long, Dd_lat)
% by the same shared model (green = directional (A*,B*); red = full linParamSys),
% so a_long acts along the road tangent and a_lat perpendicular to it. Each tube
% slice is then mapped vertex-by-vertex through the Frenet transform, so the tube
% bends along the road in world coordinates.
%
% Two cases from ONE engine (run it twice):
%     overtake_frenet_animation('straight')   % straight reference (as before)
%     overtake_frenet_animation('curved')     % a sweeping left turn
%
% NB1  On a straight reference the Frenet map collapses to the earlier
%      translate-and-place behaviour (built-in sanity check).
% NB2  Mapping a flat (s,d) reach set through the curvilinear transform is the
%      standard road-network approach, but it stretches area by (1-kappa*d) on
%      curves -- exact enough for gentle turns, mildly distorted on sharp ones.
% NB3  Ego motion is scripted (overtake in green, wait in red); the sets are
%      drawn but do not gate the motion.
%
% Requires: CORA on path (zonotope, linearSys, linParamSys, matZonotope, reach,
% project, vertices, reduce).

clc; close all;
if nargin<1 || isempty(roadType), roadType = 'curved'; end

%% ============================ TOGGLES ================================
MAKE_VIDEO = true;  FPS = 20;

%% ====================== SCENARIO (shared) ===========================
cfg.laneWidth   = 3.75;
cfg.dLeft = +cfg.laneWidth;  cfg.dMid = 0;  cfg.dRight = -cfg.laneWidth;
cfg.L = 4.7;  cfg.W = 1.9;                 % vehicle length / width
cfg.dt = 0.05;  cfg.T = 9.0;               % animation step / horizon
cfg.predHorizon = 2.0;                     % live prediction shown each frame [s]
cfg.safeGap = 12;

% initial ARC-LENGTH positions [m] along the reference and speeds [m/s]
ego.s0 = 5;    cfg.vEgo  = 22;             % middle lane, catching up
o1.s0  = 30;   cfg.vObs1 = 12;   o1.dLane = cfg.dLeft;    % LEFT lane (overtaken)
o2.s0  = 26;   cfg.vObs2 = 13;   o2.dLane = cfg.dRight;   % RIGHT lane

%% ====================== REFERENCE CURVE =============================
switch lower(roadType)
    case 'straight', fr = makeStraightRef(320);
    case 'curved',   fr = makeCurvedRef(320);
    otherwise, error('roadType must be ''straight'' or ''curved''.');
end

% camera / figure aspect per road type
if strcmpi(roadType,'curved')
    cfg.figW=900; cfg.figH=760; cfg.camHalfX=78; cfg.lookahead=18;
else
    cfg.figW=1200; cfg.figH=430; cfg.camHalfX=60; cfg.lookahead=15;
end
cfg.camHalfY = cfg.camHalfX * cfg.figH/cfg.figW;

fprintf('Frenet demo | road = %s | reference length %.0f m\n', lower(roadType), fr.len);

%% ================= ONE shared model & both tubes ====================
m = reachModel(cfg);
o1.greenPoly = tubePolysGreen(m, cfg.vObs1);  o1.redPoly = tubePolysRed(m, cfg.vObs1);
o2.greenPoly = tubePolysGreen(m, cfg.vObs2);  o2.redPoly = tubePolysRed(m, cfg.vObs2);

%% ===================== ego trajectories (Frenet) ====================
t        = 0:cfg.dt:cfg.T;
egoGreen = egoOvertakeFrenet(ego, o1, fr, cfg, t);
egoWait  = egoWaitFrenet(    ego, o1, o2, fr, cfg, t);

%% ============================ ANIMATE ==============================
tagPre = lower(roadType);
animateWindow(1, sprintf('%s road  |  GREEN (directional)  -  ego overtakes', upper(roadType)), ...
    'green', egoGreen, o1, o2, fr, cfg, t, MAKE_VIDEO, FPS, [tagPre '_window1_green']);

animateWindow(2, sprintf('%s road  |  RED (standard reach)  -  ego must wait', upper(roadType)), ...
    'red',   egoWait,  o1, o2, fr, cfg, t, MAKE_VIDEO, FPS, [tagPre '_window2_red']);

end % ===================== end main ==================================


%% ====================== REFERENCE BUILDERS ==========================
function fr = makeStraightRef(smax)
ds=0.5; s=(0:ds:smax)';
fr = packRef(s, s, zeros(size(s)), zeros(size(s)));
end

function fr = makeCurvedRef(smax)
% straight -> smooth left turn -> straight. Heading ramps 0..thTot over the turn.
ds=0.5; s=(0:ds:smax)';
sTurnStart=30; sTurnEnd=200; thTot=deg2rad(75);
seg = (s - sTurnStart)./(sTurnEnd - sTurnStart);
th  = thTot .* smoothstep(seg);          % smoothstep clamps -> flat before/after
X = cumtrapz(s, cos(th));  Y = cumtrapz(s, sin(th));
fr = packRef(s, X, Y, th);
end

function fr = packRef(s, X, Y, TH)
fr.s=s(:); fr.X=X(:); fr.Y=Y(:); fr.TH=TH(:); fr.len=s(end);
end

function W = frWorld(fr, sq, d)
% (sq,d) -> world. sq any shape; d scalar or same shape. Returns 2xN.
x  = interp1(fr.s, fr.X , sq, 'linear','extrap');
y  = interp1(fr.s, fr.Y , sq, 'linear','extrap');
th = interp1(fr.s, fr.TH, sq, 'linear','extrap');
W  = [x(:).' + d.*(-sin(th(:).')); y(:).' + d.*( cos(th(:).'))];
end

function th = frHeading(fr, sq)
th = interp1(fr.s, fr.TH, sq, 'linear','extrap');
end


%% ===================== SHARED REACH-SET MODEL =======================
function m = reachModel(cfg)
% ONE model for BOTH tubes -- IDENTICAL matrices, uncertainties, R0 and U.
% Only the downstream COMPUTATION differs (green optimises A*,B*; red = full
% non-directional parametric reach).
c_x0=0.10; c_y0=0.50; b_x0=1.00; b_y0=1.00;
dc_x=0.05; dc_y=0.20; db_x=0.10; db_y=0.15; uscale=3;
m.Ac=[0 0 1 0;0 0 0 1;0 0 -c_x0 0;0 0 0 -c_y0];
m.GA=zeros(4,4,2); m.GA(3,3,1)=-dc_x; m.GA(4,4,2)=-dc_y; m.GA=m.GA*uscale;
m.Bc=[0 0;0 0;b_x0 0;0 b_y0];
m.GB=zeros(4,2,2); m.GB(3,1,1)=db_x; m.GB(4,2,2)=db_y; m.GB=m.GB*uscale;
m.ddir=[1;0;0;0]; m.r0=[0.4;0.4;1.0;0.5];
m.a_long_max=2.0; m.a_lat_max=0.6;     % SAME for both
m.dt=0.20; m.tf=cfg.predHorizon; m.maxit=20; m.tol=1e-8;
end


%% ===================== TUBE BUILDERS (Frenet displacement) ==========
function polys = tubePolysGreen(m, v0)
polys = slicesToPolys(computeGreenTube(m, v0), [1 2]);
end
function polys = tubePolysRed(m, v0)
c0=zeros(4,1); c0(3)=v0;
X0=zonotope(c0,diag(m.r0));  U=zonotope([0;0],diag([m.a_long_max;m.a_lat_max]));
params.tFinal=m.tf; params.R0=X0; params.U=U;
options.timeStep=m.dt; options.taylorTerms=4; options.zonotopeOrder=100; options.intermediateTerms=4;
A=matZonotope(m.Ac,m.GA); B=matZonotope(m.Bc,m.GB);
R=reach(linParamSys(A,B,'constParam'), params, options);
polys = slicesToPolys(R.timeInterval.set,[1 2]);
end
function polys = slicesToPolys(Rslices, dims)
% Each slice -> a smooth closed boundary in (Ds,Dd), ready to map through Frenet.
n=numel(Rslices); polys=cell(n,1);
for k=1:n
    Z2 = reduce(project(Rslices{k},dims),'girard',10);
    V  = vertices(Z2);
    if size(V,2)>=3, K=convhull(V(1,:),V(2,:)); V=V(:,K); end
    if size(V,2)>=2 && isequal(V(:,1),V(:,end)), V(:,end)=[]; end
    polys{k} = resampleClosed(V, 80);      % densify so curved edges render smoothly
end
end
function Q = resampleClosed(P, N)
if size(P,2)<2, Q=repmat(P(:,1),1,N); return; end
V=[P, P(:,1)];
d=[0, cumsum(hypot(diff(V(1,:)),diff(V(2,:))))];
if d(end)==0, Q=repmat(P(:,1),1,N); return; end
q=linspace(0,d(end),N);
Q=[interp1(d,V(1,:),q); interp1(d,V(2,:),q)];
end

function Rslices = computeGreenTube(mdl, v0)
Ac=mdl.Ac; GA=mdl.GA; Bc=mdl.Bc; GB=mdl.GB; d=mdl.ddir; tf=mdl.tf;
c0=zeros(size(Ac,1),1); c0(3)=v0; r0=mdl.r0;
uc=[0;0]; ur=[mdl.a_long_max; mdl.a_lat_max];
VX=boxVertices2(c0,r0); VU=boxVertices2(uc,ur);
nA=size(GA,3); A_best=Ac; J_best=-inf; P0=expm(Ac.'*tf)*d;
for it=1:mdl.maxit
    bJ=-inf; A_it=Ac; x0_it=VX(:,1);
    for v=1:size(VX,2)
        xv=VX(:,v); c=zeros(nA,1);
        for i=1:nA, c(i)=P0.'*GA(:,:,i)*xv; end
        al=sign(c); al(al==0)=1; A_try=Ac;
        for i=1:nA, A_try=A_try+al(i)*GA(:,:,i); end
        J=P0.'*A_try*xv; if J>bJ, bJ=J; A_it=A_try; x0_it=xv; end
    end
    Jt=d.'*expm(A_it*tf)*A_it*x0_it;
    if Jt>J_best, J_best=Jt; A_best=A_it; end
    P0n=expm(A_it.'*tf)*d; if norm(P0n-P0)<mdl.tol, break; end; P0=P0n;
end
A_star=A_best; P0s=expm(A_star.'*tf)*d;
nB=size(GB,3); B_star=Bc; bJ=-inf;
for v=1:size(VU,2)
    uv=VU(:,v); c=zeros(nB,1);
    for j=1:nB, c(j)=P0s.'*GB(:,:,j)*uv; end
    be=sign(c); be(be==0)=1; B_try=Bc;
    for j=1:nB, B_try=B_try+be(j)*GB(:,:,j); end
    J=P0s.'*B_try*uv; if J>bJ, bJ=J; B_star=B_try; end
end
X0=zonotope(c0,diag(r0)); U=zonotope(uc,diag(ur));
params.tFinal=tf; params.R0=X0; params.U=U;
options.timeStep=mdl.dt; options.taylorTerms=4; options.zonotopeOrder=100; options.linAlg='standard';
Rslices = reach(linearSys('green',A_star,B_star),params,options).timeInterval.set;
end
function V = boxVertices2(c,r)
n=numel(c); K=2^n; V=zeros(n,K);
for k=0:K-1, b=bitget(k,1:n); sgn=2*b-1; V(:,k+1)=c(:)+sgn(:).*r(:); end
end


%% ===================== EGO TRAJECTORIES (Frenet) ====================
function tr = egoOvertakeFrenet(ego, o1, fr, cfg, t)
s = ego.s0 + cfg.vEgo.*t;
t_pass  = (o1.s0 - ego.s0)/(cfg.vEgo - cfg.vObs1);
t_merge = t_pass + 0.6; dur = 1.6;
a = smoothstep((t - t_merge)./dur);
d = cfg.dMid + a.*(cfg.dLeft - cfg.dMid);
W = frWorld(fr, s, d);
tr = poseFromWorld(W(1,:), W(2,:), t, s, d);
end
function tr = egoWaitFrenet(ego, o1, o2, fr, cfg, t)
lead = min(o1.s0 + cfg.vObs1.*t, o2.s0 + cfg.vObs2.*t);
free = ego.s0 + cfg.vEgo.*t;
s = cummax(min(free, lead - cfg.safeGap));
d = cfg.dMid + 0*t;
W = frWorld(fr, s, d);
tr = poseFromWorld(W(1,:), W(2,:), t, s, d);
end
function tr = poseFromWorld(x, y, t, s, d)
dt=t(2)-t(1); dx=gradient(x,dt); dy=gradient(y,dt);
th=atan2(dy,dx); th(~isfinite(th))=0;
tr=struct('t',t,'x',x,'y',y,'theta',th,'v',hypot(dx,dy),'s',s,'d',d);
end


%% ========================= ANIMATION ================================
function animateWindow(figNo, ttl, which, ego, o1, o2, fr, cfg, t, makeVideo, fps, tag)
fig=figure(figNo); clf(fig); set(fig,'Color','w','Position',[60 60 cfg.figW cfg.figH]);
ax=axes(fig);
if strcmp(which,'green'), col=[0.10 0.75 0.20]; else, col=[0.90 0.20 0.20]; end

vw=[];
if makeVideo
    try,  vw=VideoWriter(tag,'MPEG-4');
    catch, vw=VideoWriter([tag '.avi'],'Motion JPEG AVI'); end
    vw.FrameRate=fps; open(vw);
end

for k=1:numel(t)
    tk=t(k);
    cla(ax); hold(ax,'on'); set(ax,'DataAspectRatio',[1 1 1]);

    % road, clipped to a window of arclength around the ego
    drawRoadFrenet(ax, fr, cfg.laneWidth, ego.s(k)-60, ego.s(k)+95);

    % obstacle Frenet -> world poses now
    s1=o1.s0+cfg.vObs1*tk; W1=frWorld(fr,s1,o1.dLane); th1=frHeading(fr,s1);
    s2=o2.s0+cfg.vObs2*tk; W2=frWorld(fr,s2,o2.dLane); th2=frHeading(fr,s2);

    % live predicted reach tubes, mapped through the Frenet transform
    if strcmp(which,'green')
        drawTubeFrenet(ax, o1.greenPoly, s1, o1.dLane, fr, col);
        drawTubeFrenet(ax, o2.greenPoly, s2, o2.dLane, fr, col);
    else
        drawTubeFrenet(ax, o1.redPoly,   s1, o1.dLane, fr, col);
        drawTubeFrenet(ax, o2.redPoly,   s2, o2.dLane, fr, col);
    end

    plot(ax, ego.x, ego.y, '--','Color',[0.20 0.40 0.90],'LineWidth',0.8);  % planned path

    drawBox(ax, W1(1),W1(2),th1, cfg.L,cfg.W, [0.15 0.35 0.95]);   % left car blue
    drawBox(ax, W2(1),W2(2),th2, cfg.L,cfg.W, [0.75 0.10 0.55]);   % right car
    drawBox(ax, ego.x(k),ego.y(k),ego.theta(k), cfg.L,cfg.W, [0.05 0.05 0.05]); % ego black

    % camera: centre a little ahead of the ego along its heading
    thc=ego.theta(k); cx=ego.x(k)+cfg.lookahead*cos(thc); cy=ego.y(k)+cfg.lookahead*sin(thc);
    ax.XLim=[cx-cfg.camHalfX, cx+cfg.camHalfX];
    ax.YLim=[cy-cfg.camHalfY, cy+cfg.camHalfY];

    text(ax, ax.XLim(1)+3, ax.YLim(2)-3, egoStatus(which,ego,k), ...
        'FontSize',12,'FontWeight','bold','BackgroundColor','w','Margin',3, ...
        'VerticalAlignment','top');
    title(ax, sprintf('%s   |   t = %.2f s   (live prediction %.1f s)', ttl, tk, cfg.predHorizon), ...
        'Interpreter','none','FontWeight','normal');
    xlabel(ax,'x [m]'); ylabel(ax,'y [m]');
    drawnow;
    if makeVideo, writeVideo(vw, getframe(fig)); end
end
if makeVideo, close(vw); fprintf('saved animation: %s\n', vw.Filename); end
end

function drawTubeFrenet(ax, polys, s_obs, d_lane, fr, col)
n=numel(polys);
for k=1:n
    P = polys{k};                              % 2xN in (Ds, Dd)
    W = frWorld(fr, s_obs + P(1,:), d_lane + P(2,:));
    a = 0.08 + 0.30*(k/n);
    patch(ax,'XData',W(1,:),'YData',W(2,:),'FaceColor',col, ...
        'FaceAlpha',a,'EdgeColor',col,'LineWidth',0.4);
end
end

function str = egoStatus(which, ego, k)
if strcmp(which,'green')
    if abs(ego.d(k)-ego.d(1)) < 0.3,     str='ego: approaching the left car';
    elseif ego.d(k) < ego.d(end)-0.3,    str='ego: merging left  (overtaking)';
    else,                                str='ego: overtake complete - cruising left lane'; end
else
    if ego.v(k) < 0.7*max(ego.v),        str='ego: WAITING - red (standard reach) set blocks the corridor';
    else,                                str='ego: closing the gap ...'; end
end
end


%% ======================= DRAWING HELPERS ============================
function drawRoadFrenet(ax, fr, lw, sLo, sHi)
mask = fr.s>=sLo & fr.s<=sHi;
s = fr.s(mask).';
if numel(s)<2, s = fr.s.'; end
Etop=frWorld(fr,s,+1.5*lw); Ebot=frWorld(fr,s,-1.5*lw);
Dtop=frWorld(fr,s,+0.5*lw); Dbot=frWorld(fr,s,-0.5*lw);
grey=[0.82 0.82 0.82];
patch(ax,'XData',[Etop(1,:) fliplr(Ebot(1,:))], ...
         'YData',[Etop(2,:) fliplr(Ebot(2,:))],'FaceColor',grey,'EdgeColor','none');
plot(ax,Etop(1,:),Etop(2,:),'w-','LineWidth',2);
plot(ax,Ebot(1,:),Ebot(2,:),'w-','LineWidth',2);
plot(ax,Dtop(1,:),Dtop(2,:),'--','Color',[0.95 0.85 0.10],'LineWidth',1.2);
plot(ax,Dbot(1,:),Dbot(2,:),'--','Color',[0.95 0.85 0.10],'LineWidth',1.2);
set(ax,'Color',grey);
end

function drawBox(ax,x,y,th,L,W,col)
hl=L/2; hw=W/2; loc=[hl hw; hl -hw; -hl -hw; -hl hw]';
R=[cos(th) -sin(th); sin(th) cos(th)]; C=R*loc+[x;y];
patch(ax,'XData',C(1,:),'YData',C(2,:),'FaceColor',col,'EdgeColor','k','FaceAlpha',0.95);
end

function y=smoothstep(x), x=max(0,min(1,x)); y=3*x.^2-2*x.^3; end