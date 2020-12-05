function [Psi, dPsidr, lambdas] = Prendergast(r, rho, Rfield, min_lam, max_lam, res, vAcen)
% Constructs the Prendergast magnetic field solution, given the radial 
% density profile rho(r) and desired extent of the field, Rfield. The field
% outside will be set to zero. The inputs 'r' and 'rho' need to be vectors
% of the same size.
%
% The parameters min_lam, max_lam and res specify the search range and
% resolution to find the lambda roots. All roots found will be returned as
% the vector 'lam', and corresponding radial flux function 'Psi' as an Nlam
% x Nr array, where Nlam is the number of lambda roots found and Nr is the 
% length of 'r'. The solution will be scaled such that the central Alfven
% speed is vAcen, in units of the dynamical speed.
%
% Created 30th Sep 2020           C. Loi

rm = r(r < Rfield);
n_rm = length(rm);
rhom = rho(r < Rfield);

lam_trial = min_lam : res : max_lam;
lambdas = lambda_roots(rm, rhom, lam_trial);

if isempty(lambdas)
    disp('** WARNING ** No roots found. Search range for lambda may need adjustment.')
end

Psi = zeros(length(lambdas), length(r));
dPsidr = zeros(length(lambdas), length(r));

j1 = @(x) sin(x) ./ x.^2 - cos(x) ./ x;
dj1dx = @(x) 2*cos(x) ./ x.^2 - 2*sin(x) ./ x.^3 + sin(x) ./ x;
f = @(x1, x2) cos(x2 - x1) .* (1 ./ (x1.^2 .* x2) - 1 ./ (x1 .* x2.^2)) ...
    - sin(x2 - x1) .* (1 ./ (x1.^2 .* x2.^2) + 1 ./ (x1 .* x2));
dfdr1 = @(x1, x2, lam) lam * cos(x2 - x1) .* (2 ./ (x1.^2 .* x2.^2) ...
    - 2 ./ (x1.^3 .* x2) + 1 ./ (x1 .* x2)) + lam * sin(x2 - x1) .* ...
    (2 ./ (x1.^3 .* x2.^2) + 2 ./ (x1.^2 .* x2) - 1 ./ (x1 .* x2.^2));

intgnd1 = @(lam) rhom .* rm.^3 .* j1(lam*rm);
intgnd2 = @(lam) rhom .* rm.^3 .* f(lam*rm, lam*Rfield);

for i = 1:length(lambdas)
    lam = lambdas(i);
    
    I1 = cumtrapz(intgnd1(lam));
    I2 = fliplr(cumtrapz(fliplr(intgnd2(lam))));
    
    Psi(i, 1:n_rm) = lam * rm / j1(lam*Rfield) .* ...
        (f(lam*rm, lam*Rfield) .* I1 + j1(lam*rm) .* I2);
    dPsidr(i, 1:n_rm) = lam / j1(lam*Rfield) .* ( (f(lam*rm, lam*Rfield) ...
        + rm .* dfdr1(lam*rm, lam*Rfield, lam)) .* I1 + ...
        (j1(lam*rm) + lam*rm .* dj1dx(lam*rm)) .* I2 );
    
    % Scale to desired strength
    r0 = min(r);
    Psi0 = Psi(i, r==min(r));
    rho0 = rho(r==min(r));
    scale = 0.5 * r0^2 * sqrt(rho0) / Psi0;
    Psi(i,:) = scale * Psi(i,:) * vAcen;
    dPsidr(i,:) = scale * dPsidr(i,:) * vAcen;
end


function lambdas = lambda_roots(xi, rhox, trial)
% Solve for the lambda roots given the trial vector 'trial', radial
% vector 'xi' and density profile 'rhox', for a field occupying the whole
% volume within radius max(xi).

max_iters = 100;     % maximum number of iterations before giving up
frac_tol = 1e-4;     % convergence tolerance for fractional change in lambda

intgnd = @(lam) rhox .* (xi / lam^2 .* sin(lam*xi) - xi.^2 / lam .* cos(lam*xi));

intVals = zeros(size(trial));
for i = 1:length(trial)
    intVals(i) = sum(intgnd(trial(i)));
end

rootLocs = intVals(1:end-1) .* intVals(2:end) < 0;
nRoots = sum(rootLocs);
lam_low = trial([rootLocs false]);
lam_high = trial([false rootLocs]);
int_low = intVals([rootLocs false]);
int_high = intVals([false rootLocs]);

converged = false(1, nRoots);
lambdas = zeros(1, nRoots);

for i = 1:nRoots
    lam1 = lam_low(i);
    lam2 = lam_high(i);
    int1 = int_low(i);
    int2 = int_high(i);
    
    old = NaN;
    n_iters = 1;
    
    while ~converged(i)
        lam3 = (lam1 * int2 - lam2 * int1) / (int2 - int1);
        int3 = sum(intgnd(lam3));
        
        if int3 * int1 < 0
            lam2 = lam3;
            int2 = int3;
        else
            lam1 = lam3;
            int1 = int3;
        end
        
        if abs((lam3 - old) / lam3) < frac_tol, converged(i) = true; end
        old = lam3;
        
        if n_iters > max_iters, break; end
        n_iters = n_iters + 1;
    end
    if n_iters > max_iters, continue; end
    
    lambdas(i) = lam3;
end
