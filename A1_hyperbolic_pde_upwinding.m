%% 01_hyperbolic_pde_upwinding.m
% Demonstrates first-order upwind, Lax-Friedrichs, and MUSCL (limited)
% Relevant to: Staff CFD Modeling Engineer
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
u_upwind = u0;
u_lf     = u0;
u_muscl  = u0;

%% Time marching
for n = 1:Nt
    % --- First-order upwind (a > 0) ---
    % F_{i-1/2} = a * u_{i-1}
    uL = [u_upwind(end); u_upwind(1:end-1)];
    u_upwind = u_upwind - (dt/dx) * a .* (u_upwind - uL);

    % --- Lax-Friedrichs (classic) ---
    uL = [u_lf(end); u_lf(1:end-1)];
    uR = [u_lf(2:end); u_lf(1)];
    u_lf = 0.5*(uL + uR) - 0.5*(dt/dx)*a*(uR - uL);

    % --- MUSCL with minmod limiter (2nd-order TVD) ---
    % Correct limited slope: sigma_i = minmod(u_i - u_{i-1}, u_{i+1} - u_i)
    % Left state at face i+1/2:  u_i + 0.5 * sigma_i
    uL_cell = [u_muscl(end); u_muscl(1:end-1)];   % u_{i-1}
    uR_cell = [u_muscl(2:end); u_muscl(1)];       % u_{i+1}
    deltaL  = u_muscl - uL_cell;                  % backward difference
    deltaR  = uR_cell - u_muscl;                  % forward difference
    % minmod limiter
    sigma   = 0.5*(sign(deltaL)+sign(deltaR)) .* min(abs(deltaL), abs(deltaR));
    % Interface state from the left (used because a > 0)
    u_face  = u_muscl + 0.5 * sigma;              % state at i+1/2 from left
    % Flux at i+1/2
    flux    = a * u_face;
    % Conservative update: (F_{i+1/2} - F_{i-1/2})/dx
    flux_L  = [flux(end); flux(1:end-1)];         % F_{i-1/2}
    u_muscl = u_muscl - (dt/dx) * (flux - flux_L);
end

%% Exact solution (periodic shift)
shift = mod(round(a*Tfinal/dx), Nx);
u_exact = circshift(u0, shift);

%% Plot
figure('Position',[100 100 900 500], 'Name','Hyperbolic PDE Upwinding');
plot(x, u0, 'k--', 'LineWidth',1.2); hold on;
plot(x, u_exact, 'k-', 'LineWidth',1.5);
plot(x, u_upwind, 'b-', 'LineWidth',1.2);
plot(x, u_lf, 'r-', 'LineWidth',1.2);
plot(x, u_muscl, 'g-', 'LineWidth',1.5);
legend('Initial','Exact','1st-order Upwind','Lax-Friedrichs','MUSCL (minmod)', ...
       'Location','best');
xlabel('x'); ylabel('u');
title(sprintf('1-D Advection after t=%.1f (CFL=%.2f, Nx=%d)', Tfinal, CFL, Nx));
grid on; box on;
xlim([0 1]);

fprintf('\n=== Hyperbolic PDE Demo complete ===\n');
fprintf('Favorite practical choice for many industrial CFD codes: limited MUSCL / TVD schemes\n');
fprintf('or high-resolution Godunov-type methods with approximate Riemann solvers.\n');
fprintf('First-order upwind is robust but diffusive; pure high-order without limiting oscillates.\n');
