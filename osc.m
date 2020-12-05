%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Solves the eigenvalue problem of linear adiabatic stellar oscillations under
%  the Cowling approximation using the shooting method. Integration is performed
%  using RK4 from each end of the grid, applying appropriate boundary conditions
%  and matched at an interior point via N-R iteration on the determinant of the
%  matching matrix.
%
%  The background stellar model needs to be supplied, in the form of radial
%  profiles of:
%    - mass density, rho
%    - pressure, p
%    - adiabatic index, gam
%    - squared buoyancy frequency, N2
%
%  Requires the use of the following functions/routines:
%    - polytrope.m
%    - get_inout_solns.m
%    - EckartClass.m
%
%  Created 29th Sep 2020               C. Loi
%
%%%%%

clear

%% USER-DEFINED PARAMETERS

l = 1;                 % SH degree (currently available: l=0,1,2,3)

% Stellar model parameters
modelType = 'polytrope';   % on offer: 'polytrope', 'mesa'
mesaInFile = 'MESA_2Msun_profile26';
npoly = 4.2;             % polytropic index
Rstar = 6;             % polytrope radius (units of solar radii)
Mstar = 2;             % polytrope mass (units of solar masses)
dxi = 0.001;          % polytrope grid spacing (e.g. for n=3, dxi=0.0005 gives ~10000 points)

% Mode searching parameters
min_om = 0.5;         % lower bound of search
max_om = 1;         % upper bound
dom = 0.01;            % resolution
max_iters = 100;       % give up after exceeding this number of iterations
frac_tol = 1e-4;       % consider converged if normalised determinant falls below this

% Turn on plotting?
wantPlots = true;

% Save results?
saveResults = false;
if strcmp(modelType, 'polytrope')
    outFile = ['modes_n' num2str(npoly) 'poly_l' num2str(l) '.mat'];
elseif strcmp(modelType, 'mesa')
    outFile = ['modes_' mesaInFile '_l' num2str(l) '.mat'];
end


%% STELLAR MODEL CONSTRUCTION

if strcmp(modelType, 'polytrope')
    [r, rho, p, N2, gam] = polytrope(npoly, Rstar, Mstar, dxi);
    
    % Throw away innermost point (avoid potential issues with divergent terms)
    r = r(2:end);
    rho = rho(2:end);
    p = p(2:end);
    gam = gam(2:end);
    N2 = N2(2:end);

elseif strcmp(modelType, 'mesa')
    load([mesaInFile '.mat'])
end

dr = r(2) - r(1);
nPts_th = 100;
dth = pi/nPts_th;
th = 0.5*dth : dth : pi-0.5*dth;
[~, Mth] = meshgrid(r, th);

% Whatever the above is, by the end of this section need to have rho(r), 
% p(r), gam(r) and N2(r) defined


%% COMPUTE AUXILLIARY QUANTITIES

if ~exist('drho', 'var'), drho = gradient(rho, dr); end
if ~exist('d2rho', 'var'), d2rho = gradient(drho, dr); end
if ~exist('dp', 'var'), dp = gradient(p, dr); end
if ~exist('d2p', 'var'), d2p = gradient(dp, dr); end
if ~exist('d3p', 'var'), d3p = gradient(d2p, dr); end
if ~exist('dN2', 'var'), dN2 = gradient(N2, dr); end
if ~exist('d2N2', 'var'), d2N2 = gradient(dN2, dr); end


%% SCAN FREQUENCY RANGE FOR MODES

r_fit = 0.5;
i_fit = find(r >= r_fit, 1);
if mod(i_fit, 2) == 0, i_fit = i_fit+1; end      % needs to be odd

om_trial = min_om : dom : max_om;
determinants = zeros(size(om_trial));

disp('Scanning for modes...')
msg_no = 1;
for i = 1:length(om_trial)
    om = om_trial(i);
    
    get_inout_solns;
    
    determinants(i) = soln1(end,2) * soln2(1,1) - soln2(1,2) * soln1(end,1);
    
    if i/length(om_trial) > msg_no * 0.1
        fprintf('...%i%%', msg_no*10)
        msg_no = msg_no + 1;
    end
end
fprintf('\n')

rootLocs = determinants(2:end) .* determinants(1:end-1) < 0;
nRoots = sum(rootLocs);

disp(['Result: ' num2str(nRoots) ' modes found between om = ' num2str(min_om) ' and ' num2str(max_om)])


%% SOLVE FOR MODES

disp('  Order     Frequency')
disp('----------------------------')

om_low = om_trial([rootLocs false]);
om_high = om_trial([false rootLocs]);
det_low = determinants([rootLocs false]);
det_high = determinants([false rootLocs]);

% Output arrays
orders = nan(1, nRoots);
frequencies = nan(1, nRoots);
if mod(length(r), 2) == 0
    R_fns = zeros(length(r)/2, nRoots);
    H_fns = zeros(length(r)/2, nRoots);
else
    R_fns = zeros((length(r)+1)/2, nRoots);
    H_fns = zeros((length(r)+1)/2, nRoots);
end
converged = false(1, nRoots);

for i = 1:nRoots
    % Newton-Raphson iteration
    om1 = om_low(i);
    om2 = om_high(i);
    det1 = det_low(i);
    det2 = det_high(i);
    
    n_iters = 1;
    while ~converged(i)
        om = om1 - det1 * (om2 - om1) / (det2 - det1);
        
        get_inout_solns;
        
        det_i = soln1(end,2) * soln2(1,1) - soln2(1,2) * soln1(end,1);
        norm_i = soln1(end,2) * soln2(1,1) + soln2(1,2) * soln1(end,1);
        
        if det_i * det2 < 0
            om1 = om; 
            det1 = det_i;
        else
            om2 = om; 
            det2 = det_i;
        end

        if abs(det_i / norm_i) < frac_tol, converged(i) = true; end
        
        if n_iters > max_iters, break; end
        n_iters = n_iters + 1;
    end
    
    % Skip non-convergent modes
    if n_iters > max_iters
        disp('  (non-converged)')
        continue; 
    end
    
    % Stitch together solution
    soln1 = soln1 / soln1(end,1) * soln2(1,1);
    R_i = [soln1(:,1); soln2(2:end,1)];
    H_i = [soln1(:,2); soln2(2:end,2)];
    
    R_fns(:,i) = R_i;
    H_fns(:,i) = H_i;
    orders(i) = EckartClass(r(1:2:end), R_i, H_i);
    frequencies(i) = om;
    
    disp(['  ' num2str(orders(i)) '     ' num2str(om)])
end


%% VISUALISE RESULTS

if wantPlots
    figure
    plot(orders(converged), frequencies(converged), 'k*')
    xlabel('Radial order (Eckart scheme)'), ylabel('\omega')
    
    for i = 1:nRoots
        if mod(i-1,15)+1 == 1
            figure; 
            set(gcf, 'Position', [1860 200 1500 900])
        end
        
        if converged(i)
            subplot(3, 5, mod(i-1,15)+1)
            plot(r(1:2:end), R_fns(:,i), 'k', r(1:2:end), H_fns(:,i), 'r')
            ax_scale = max(median(abs(R_fns(:,i))), median(abs(H_fns(:,i))));
            axis([0 1 -10*ax_scale 10*ax_scale])
            
            xlabel('r/R_*'), ylabel('Displacement')
            title(['n = ' num2str(orders(i)) 10 ...
                '\omega = ' num2str(frequencies(i))])
        end
    end
end


%% SAVE RESULTS

if saveResults
    r = r(1:2:end);
    rho = rho(1:2:end);
    save(outFile, 'orders', 'frequencies', 'r', 'rho', 'R_fns', 'H_fns')
end
