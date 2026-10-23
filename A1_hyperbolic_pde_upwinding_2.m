%% A1_hyperbolic_pde_upwinding.m
% Demonstrates first-order upwind, Lax-Friedrichs, and three MUSCL variants
% (minmod, van Leer, MC) + limiter function plot
%
% PDE:  ∂u/∂t + a * ∂u/∂x = 0 ,  a > 0
% Domain: periodic [0,1]
% Initial condition: square wave + smooth Gaussian

clear; close all; clc;

%% Parameters
a      = 1.0;          % advection speed
Nx     = 200;          % cells
L      = 1.0;
dx     = L / Nx;
x      = (0.5:Nx-0.5)' * dx;   % cell centers
CFL    = 0.4;
dt     = CFL * dx / a;
Tfinal = 1.0;          % one full period
Nt     = round(Tfinal / dt);
dt     = Tfinal / Nt;  % exact final time

%% Initial condition
u0 = zeros(Nx,1);
u0(x > 0.2 & x < 0.4) = 1.0;                 % square wave
u0 = u0 + 0.5 * exp(-200*(x-0.7).^2);         % Gaussian

%% Allocate solutions
u_upwind   = u0;
u_lf       = u0;
u_minmod   = u0;
u_vanleer  = u0;
u_mc       = u0;

%% Time marching
for n = 1:Nt
    % --- First-order upwind ---
    uL = [u_upwind(end); u_upwind(1:end-1)];
    u_upwind = u_upwind - (dt/dx) * a .* (u_upwind - uL);

    % --- Lax-Friedrichs ---
    uL = [u_lf(end); u_lf(1:end-1)];
    uR = [u_lf(2:end); u_lf(1)];
    u_lf = 0.5*(uL + uR) - 0.5*(dt/dx)*a*(uR - uL);

    % --- MUSCL-minmod ---
    sigma = limited_slope(u_minmod, 'minmod');
    u_face = u_minmod + 0.5 * sigma;
    flux   = a * u_face;
    flux_L = [flux(end); flux(1:end-1)];
    u_minmod = u_minmod - (dt/dx)*(flux - flux_L);

    % --- MUSCL-van Leer ---
    sigma = limited_slope(u_vanleer, 'vanleer');
    u_face = u_vanleer + 0.5 * sigma;
    flux   = a * u_face;
    flux_L = [flux(end); flux(1:end-1)];
    u_vanleer = u_vanleer - (dt/dx)*(flux - flux_L);

    % --- MUSCL-MC ---
    sigma = limited_slope(u_mc, 'mc');
    u_face = u_mc + 0.5 * sigma;
    flux   = a * u_face;
    flux_L = [flux(end); flux(1:end-1)];
    u_mc = u_mc - (dt/dx)*(flux - flux_L);
end

%% Exact solution
shift   = mod(round(a*Tfinal/dx), Nx);
u_exact = circshift(u0, shift);

%% ==================================================================
%  FIGURE 1 – Solution comparison
% ==================================================================
figure('Position',[80 80 1000 520], 'Name','Solutions');
plot(x, u0,      'k--', 'LineWidth',1.2); hold on;
plot(x, u_exact, 'k-',  'LineWidth',1.8);
plot(x, u_upwind,'b-',  'LineWidth',1.0);
plot(x, u_lf,    'r-',  'LineWidth',1.0);
plot(x, u_minmod,'g-',  'LineWidth',1.4);
plot(x, u_vanleer,'m-','LineWidth',1.5);
plot(x, u_mc,    'c-',  'LineWidth',1.5);

legend('Initial','Exact',...
       '1st-order Upwind','Lax-Friedrichs',...
       'MUSCL-minmod','MUSCL-van Leer','MUSCL-MC', ...
       'Location','best');
xlabel('x'); ylabel('u');
title(sprintf('1-D Advection after t=%.1f (CFL=%.2f, Nx=%d)', Tfinal, CFL, Nx));
grid on; box on; xlim([0 1]);
plotedit off; 
zoom on;
%% ==================================================================
%  FIGURE 2 – Limiter functions φ(r)
% ==================================================================
r = linspace(-0.5, 4.5, 600);          % ratio range

% Classical limiter definitions (Sweby form)
phi_minmod   = max(0, min(1, r));
phi_vanleer  = (r + abs(r)) ./ (1 + abs(r));          % = 0 for r<0
phi_mc       = max(0, min([2*r; 0.5*(1+r); 2*ones(size(r))]));
phi_superbee = max(0, max( min(2*r,1), min(r,2) ));   % for reference

% TVD region boundaries
tvd_lower = phi_minmod;
tvd_upper = phi_superbee;

figure('Position',[120 120 780 520], 'Name','Limiter functions');
hold on;

% Shade the TVD region
fill([r, fliplr(r)], [tvd_lower, fliplr(tvd_upper)], ...
     [0.92 0.95 1.0], 'EdgeColor','none', 'FaceAlpha',0.6);
plot(r, tvd_lower, 'k:', 'LineWidth',1.0);
plot(r, tvd_upper, 'k:', 'LineWidth',1.0);

% The three limiters we use + superbee for reference
plot(r, phi_minmod,  'g-',  'LineWidth',2.0);
plot(r, phi_vanleer, 'm-',  'LineWidth',2.0);
plot(r, phi_mc,      'c-',  'LineWidth',2.0);
plot(r, phi_superbee,'k--', 'LineWidth',1.3);

% Reference lines
yline(1, 'k-', 'LineWidth',0.6);
xline(0, 'k-', 'LineWidth',0.6);
xline(1, 'k:', 'LineWidth',0.8);

xlabel('r  (ratio of successive differences)');
ylabel('\phi(r)');
title('Limiter functions  \phi(r)   (shaded = TVD region)');
legend('TVD region',...
       'minmod (lower TVD bound)',...
       'superbee (upper TVD bound)',...
       'minmod',...
       'van Leer',...
       'MC',...
       'superbee (reference)',...
       'Location','northwest');
grid on; box on;
xlim([-0.4 4.2]);
ylim([-0.1 2.3]);
set(gca,'FontSize',11);

fprintf('\n=== Demo complete ===\n');
fprintf('Figure 1 : solution comparison\n');
fprintf('Figure 2 : limiter functions (Sweby diagram)\n');
fprintf('  green  = minmod   (most dissipative)\n');
fprintf('  magenta= van Leer (smooth)\n');
fprintf('  cyan   = MC       (monotonized central)\n');

%% ==================================================================
%  LOCAL FUNCTION (must be at the end of the file)
% ==================================================================
function sigma = limited_slope(u, limiter)
    uL_cell = [u(end); u(1:end-1)];   % u_{i-1}
    uR_cell = [u(2:end); u(1)];       % u_{i+1}
    deltaL  = u - uL_cell;            % backward difference
    deltaR  = uR_cell - u;            % forward difference

    switch lower(limiter)
        case 'minmod'
            sigma = 0.5*(sign(deltaL)+sign(deltaR)) .* min(abs(deltaL), abs(deltaR));

        case 'vanleer'
            prod = deltaL .* deltaR;
            sigma = zeros(size(u));
            mask  = prod > 0;
            sigma(mask) = 2*prod(mask) ./ (deltaL(mask) + deltaR(mask));

        case 'mc'
            central = 0.5*(deltaL + deltaR);
            sigma   = 0.5*(sign(deltaL)+sign(deltaR)) .* ...
                      min( abs(central), min(2*abs(deltaL), 2*abs(deltaR)) );
        otherwise
            error('Unknown limiter');
    end
end