%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Code to load, process and save MESA output profiles into a MATLAB
%  variable file. Processing includes both smoothing and regridding to a
%  uniform grid. Assumes the following quantities are present:
%   * mass - enclosed mass in units of solar masses
%   * logR - log10 of the radius in units of solar radii
%   * logRho - log10 of the density in g/cm^3
%   * logP - log10 of the pressure in Ba
%   * gamma1 - adiabatic index
%   * log_brunt_N2_dimensionless - log10 of the squared buoyancy frequency
%       in units of 3GM*/R*^3
%
%  The non-dimensional forms that will be saved in the final MATLAB
%  variable file are done according to the following scheme:
%   * r - units of R*
%   * rho - units of M*/R*^3
%   * p - units of GM*^2/R*^4
%   * N2 - units of GM*/R*^3
%
%  Besides the above, derivatives of rho, p and N2 will also be calculated
%  and saved (these are needed by basisfree_osc.m to compute the Lorentz
%  terms).
%
%  Also, the direction of the arrays will be flipped to be monotonically
%  increasing in r, unlike the output of MESA which is decreasing.
%
%  Created 23rd Oct 2020              C. Loi
%
%%%%%

clear

%% USER-SPECIFIED PARAMETERS

inDir = '/store/ASTRO/stl36/mesa-r11701/2Msun/LOGS/';
profileNo = 31;

saveResults = true;
outFile = ['MESA_2Msun_profile' num2str(profileNo) '.mat'];

finalPts = 1000000;    % desired #points in final interpolated grid

boxcar = 13;          % #points in smoothing window (must be odd)

doPlots = false;


%% PHYSICAL CONSTANTS

G = 6.674e-11;           % Newton's gravitational constant (N m^2 / kg^2)
Rsun = 6.955e8;          % solar radius (m)
Msun = 1.99e30;          % solar mass (kg)


%% IMPORT DATA

profile = importdata([inDir 'profile' num2str(profileNo) '.data'], ' ', 6);

for i = 1:length(profile.colheaders)
    eval([profile.colheaders{i} ' = profile.data(:,' num2str(i) ');'])
end


%% DIMENSIONALISE

Rstar = 10^max(logR) * Rsun;      % m
Mstar = max(mass) * Msun;         % kg

r_SI = 10.^flipud(logR) * Rsun;        % m
rho_SI = 10.^flipud(logRho) * 1e3;     % kg/m^3
p_SI = 10.^flipud(logP) * 0.1;         % Pa
N2_SI = 10.^flipud(log_brunt_N2_dimensionless) ...
    * 3*G*Mstar/Rstar^3;               % (rad/s)^2

if doPlots
    maxPlotX = 1e8;
    
    figure, set(gcf, 'Position', [340 350 900 740])
    
    subplot(2,2,1)
    plot(r_SI, rho_SI, 'k.-')
    xlabel('r (m)'), ylabel('\rho (kg/m^3)')
    axis([0 maxPlotX ylim])
    title('Unprocessed, in SI')
    
    subplot(2,2,2)
    plot(r_SI, p_SI, 'k.-')
    xlabel('r (m)'), ylabel('p (Pa)')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,3)
    plot(r_SI, flipud(gamma1), 'k.-')
    xlabel('r (m)'), ylabel('\gamma')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,4)
    plot(r_SI, N2_SI, 'k.-')
    xlabel('r (m)'), ylabel('N^2 (rad/s)^2')
    axis([0 maxPlotX ylim])
end


%% NON-DIMENSIONALISE

r_nonunif = r_SI / Rstar;
rho_nonunif = rho_SI * Rstar^3 / Mstar;
p_nonunif = p_SI * Rstar^4 / G / Mstar^2;
N2_nonunif = N2_SI * Rstar^3 / G / Mstar;
gam_nonunif = flipud(gamma1);

if doPlots
    maxPlotX = 0.04;
    
    figure, set(gcf, 'Position', [340 350 900 740])
    
    subplot(2,2,1)
    plot(r_nonunif, rho_nonunif, 'k.-')
    xlabel('r/R_*'), ylabel('\rho')
    axis([0 maxPlotX ylim])
    title('Unprocessed, non-dimensionalised')
    
    subplot(2,2,2)
    plot(r_nonunif, p_nonunif, 'k.-')
    xlabel('r/R_*'), ylabel('p')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,3)
    plot(r_nonunif, flipud(gamma1), 'k.-')
    xlabel('r/R_*'), ylabel('\gamma')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,4)
    plot(r_nonunif, N2_nonunif, 'k.-')
    xlabel('r/R_*'), ylabel('N^2')
    axis([0 maxPlotX ylim])
end


%% SMOOTHING

% Generalised boxcar method (see handwritten notes dated 24th June 2019)
nPts = length(r_nonunif);

rho_sm = zeros(nPts, 1);
p_sm = zeros(nPts, 1);
N2_sm = zeros(nPts, 1);
gam_sm = zeros(nPts, 1);

genboxcar = @(xs, ys, x) sum((xs-mean(xs)) .* (ys-mean(ys))) ...
    ./ sum((xs-mean(xs)).^2) * (x - mean(xs)) + mean(ys);

for i = 1:nPts
    if i < nPts/2
        jcen = max(ceil(boxcar/2), i);
    else
        jcen = min(nPts-floor(boxcar/2), i);
    end
    
    rs = r_nonunif(jcen-floor(boxcar/2) : jcen+floor(boxcar/2));
    rhos = rho_nonunif(jcen-floor(boxcar/2) : jcen+floor(boxcar/2));
    ps = p_nonunif(jcen-floor(boxcar/2) : jcen+floor(boxcar/2));
    N2s = N2_nonunif(jcen-floor(boxcar/2) : jcen+floor(boxcar/2));
    gams = gam_nonunif(jcen-floor(boxcar/2) : jcen+floor(boxcar/2));
    
    rho_sm(i) = genboxcar(rs, rhos, r_nonunif(i));
    p_sm(i) = genboxcar(rs, ps, r_nonunif(i));
    N2_sm(i) = genboxcar(rs, N2s, r_nonunif(i));
    gam_sm(i) = genboxcar(rs, gams, r_nonunif(i));
end

if doPlots
    maxPlotX = 0.04;
    
    figure, set(gcf, 'Position', [340 350 900 740])
    
    subplot(2,2,1)
    plot(r_nonunif, rho_sm, 'k.-')
    xlabel('r/R_*'), ylabel('\rho')
    axis([0 maxPlotX ylim])
    title('Smoothed')
    
    subplot(2,2,2)
    plot(r_nonunif, p_sm, 'k.-')
    xlabel('r/R_*'), ylabel('p')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,3)
    plot(r_nonunif, gam_sm, 'k.-')
    xlabel('r/R_*'), ylabel('\gamma')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,4)
    plot(r_nonunif, N2_sm, 'k.-')
    xlabel('r/R_*'), ylabel('N^2')
    axis([0 maxPlotX ylim])
end

%% INTERPOLATION

% Linearly interpolate to uniform grid with same number of points
r_unif = linspace(min(r_nonunif), max(r_nonunif), nPts);
rho_unif = interp1(r_nonunif, rho_sm, r_unif);
p_unif = interp1(r_nonunif, p_sm, r_unif);
N2_unif = interp1(r_nonunif, N2_sm, r_unif);
gam_unif = interp1(r_nonunif, gam_sm, r_unif);

if doPlots
    maxPlotX = 0.04;
    
    figure, set(gcf, 'Position', [340 350 900 740])
    
    subplot(2,2,1)
    plot(r_unif, rho_unif, 'k.-')
    xlabel('r/R_*'), ylabel('\rho')
    axis([0 maxPlotX ylim])
    title('Smoothed, uniformly interpolated')
    
    subplot(2,2,2)
    plot(r_unif, p_unif, 'k.-')
    xlabel('r/R_*'), ylabel('p')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,3)
    plot(r_unif, gam_unif, 'k.-')
    xlabel('r/R_*'), ylabel('\gamma')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,4)
    plot(r_unif, N2_unif, 'k.-')
    xlabel('r/R_*'), ylabel('N^2')
    axis([0 maxPlotX ylim])
end

% Downsample, spline interpolate and resample to desired grid size
r = linspace(min(r_unif), max(r_unif), finalPts);

dsfac = 5;         % downsampling factor

rho = interp1(r_unif(1:dsfac:end), rho_unif(1:dsfac:end), r, 'spline');
p = interp1(r_unif(1:dsfac:end), p_unif(1:dsfac:end), r, 'spline');
N2 = interp1(r_unif(1:dsfac:end), N2_unif(1:dsfac:end), r, 'spline');
gam = interp1(r_unif(1:dsfac:end), gam_unif(1:dsfac:end), r, 'spline');

if doPlots
    maxPlotX = 0.04;
    
    figure, set(gcf, 'Position', [340 350 900 740])
    
    subplot(2,2,1)
    plot(r, rho, 'k.-')
    xlabel('r/R_*'), ylabel('\rho')
    axis([0 maxPlotX ylim])
    title('Resampled')
    
    subplot(2,2,2)
    plot(r, p, 'k.-')
    xlabel('r/R_*'), ylabel('p')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,3)
    plot(r, gam, 'k.-')
    xlabel('r/R_*'), ylabel('\gamma')
    axis([0 maxPlotX ylim])
    
    subplot(2,2,4)
    plot(r, N2, 'k.-')
    xlabel('r/R_*'), ylabel('N^2')
    axis([0 maxPlotX ylim])
end


%% CALCULATE RADIAL DERIVATIVES

dr = r(2) - r(1);

drho = gradient(rho, dr);
d2rho = gradient(drho, dr);
dp = gradient(p, dr);
d2p = gradient(dp, dr);
d3p = gradient(d2p, dr);
dN2 = gradient(N2, dr);
d2N2 = gradient(dN2, dr);

%% SAVE RESULTS

if saveResults
    save(outFile, 'r', 'rho', 'drho', 'd2rho', ...
        'p', 'dp', 'd2p', 'd3p', 'gam', 'N2', 'dN2', 'd2N2');
end