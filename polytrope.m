function [r, rho, p, N2, gam] = polytrope(n, R, M, dxi)
% Constructs a polytrope of index n, radius R (units of solar radii) and
% mass M (units of solar masses), by solving the Lane-Emden equation. The
% parameter dxi controls the overall size of the grid. Example values for
% this are as follows:
%   dxi = 0.0005 gives ~10000 points for n=3
%   dxi = 0.001 gives ~10000 points for n=4
%   dxi = 0.002 gives ~10000 points for n=4.5
%
% Outputs the following parameters:
%  * r - radial grid values, in units of the stellar radius
%  * rho - mass density, in units of M*/R*^3
%  * p - gas pressure, in units of dynamical pressure GM*^2/R*^4
%  * N2 - squared buoyancy frequency, in units of the dynamical frequency
%      sqrt(GM*/R*^3)
%  * gam - adiabatic index
%
% Created 19th Oct 2020          C. Loi

% Physical constants
G = 6.67e-11;             % Newton's gravitational constant in N m^2 kg^-2
Rsun = 6.955e8;           % solar radius in metres
Msun = 1.99e30;           % solar mass in kilograms

% Constant adiabatic index
gamma = 5/3;

% Initial conditions
xi0 = 1e-8;
X = [1-xi0^2/6, -xi0/3];       % X = [zeta, eta]

% Output arrays;
zeta = [];
eta = [];

% Integrate
xi = xi0;
while X(1) > 0
    X_next = LaneEmden(X, xi, dxi, n);
    zeta = [zeta, X_next(1)];
    eta = [eta, X_next(2)];
    
    xi = xi + dxi;
    X = X_next;
end

% Take one step back to avoid negative point
xi = xi - dxi;
zeta = zeta(1:end-1);
eta = eta(1:end-1);

xiVals = xi0:dxi:xi;
xiVals = xiVals(1:length(zeta));

% Compute auxilliary quantities
xi_star = xi;
eta_star = eta(end-1);   % !!!!!!!
Rstar = R*Rsun;
Mstar = M*Msun;
Cn = (4*pi)^(1/n)/(n+1) * (-xi_star^2*eta_star)^(1/n-1) * xi_star^(1-3/n);
K = G * Cn * Rstar^(3/n-1) * Mstar^(1-1/n);
rho_c = (4*pi*G / K / (n+1) * (Rstar/xi_star)^2)^(n/(1-n));
p_c = K * rho_c^(1+1/n);
fm = cumsum(zeta.^n .* xiVals.^2) / sum(zeta.^n .* xiVals.^2);

% Compute desired output quantities
r = xiVals / xi_star;
rho = rho_c * zeta.^n * Rstar^3 / Mstar;
p = p_c * zeta.^(n+1) * Rstar^4 / G / Mstar^2;
N2 = eta ./ zeta * ((n+1)/gamma-n) * xi_star^3 ./ xiVals.^2 .* fm;
gam = gamma * ones(size(r));


function X_new = LaneEmden(X_old, xi, dxi, n)
% Performs one integration step of the Lane-Emden equation using RK4

% RHS function
f = @(X, s) [X(2), -X(1)^n - 2*X(2)/s];

% One step of RK4
f1 = f(X_old, xi);
f2 = f(X_old + 0.5*dxi*f1, xi+0.5*dxi);
f3 = f(X_old + 0.5*dxi*f2, xi+0.5*dxi);
f4 = f(X_old + dxi*f3, xi+dxi);

X_new = X_old + dxi/6 * (f1 + 2*f2 + 2*f3 + f4);