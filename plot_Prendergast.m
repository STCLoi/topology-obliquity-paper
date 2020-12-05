function plot_Prendergast(hfig, r, Psi, lambda, nc, Rfield, cmap, dsamp)
% Generates a plot of the meridional section of the Prendergast solution,
% given input radial vector r, radial flux function Psi and the eigenvalue 
% lambda. The field components are determined from these according to e.g. 
% Eq (11) of Loi 2020 (MNRAS, 493, 5726). Plots the projection of the 
% poloidal component of the field lines as 'nc' contours and the toroidal 
% component as underlying colour onto the figure with number 'hfig', using
% colour map 'cmap'.
%
% Note that r and Psi are assumed to be 1D vectors, and lambda is
% assumed to be a scalar. The resolution will be lowered by factor dsamp.
%
% Created 24th Nov 2020              C. Loi

% Restrict to field region and downsample
rm = r(r < 1.1*Rfield);
rm = rm(1:dsamp:end);

% Set up Cartesian grid
x = linspace(0, max(rm), length(rm));
y = linspace(-max(rm), max(rm), 2*length(rm));
[Mx, My] = meshgrid(x, y);
Mr = sqrt(Mx.^2 + My.^2);
Msth = Mx ./ Mr;

Psi_fun = @(xi) interp1(r, Psi, xi);
psi = Psi_fun(Mr) .* Msth.^2;
Bphi = -lambda * Psi_fun(Mr) ./ Mr .* Msth;

figure(hfig)
imagesc(x, y, Bphi)
hold on
contour(Mx, My, psi, nc, 'k')
xlabel('x (R_*)'), ylabel('z (R_*)'), title(['\lambda = ' num2str(lambda)])
colormap(cmap)
caxis([-1 1]*max(abs(Bphi(:))))
axis equal, axis tight
set(gca, 'YDir', 'normal')
colorbar
