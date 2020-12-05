%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Perform shooting of two solutions from the endpoints of a grid to solve
%  for stellar oscillation modes. Intended for use together with osc.m,
%  the main code. Requires that the stellar backgrounds are all present in
%  the workspace, and assumes the relevant frequency is stored as the
%  variable 'om'. The outputs desired from this routine are the 1D arrays
%  'soln1' and 'soln2', representing the solutions obtained by shooting
%  from the inner and outer boundaries, respectively.
%
%  The relevant expressions for the problem can be found in Eqs (7.33) to 
%  (7.36) in my handwritten notes.
%
%  Requires the use of the following functions:
%  * shoot.m
%
%  Created 29th Sep 2020            C. Loi
%
%%%%%

global cRH cRR cHR cHH

cHR = (1 - N2 / om^2) ./ r;
cHH = -rho .* N2 ./ dp - 1 ./ r;
cRH = l*(l+1) ./ r - rho .* r ./ gam ./ p * om^2;
cRR = -2 ./ r - dp ./ gam ./ p;

% Boundary conditions (assumes non-dim scheme where R_* = 1)
BCin = [min(r)^(l-1), min(r)^(l-1)/l];
BCout = [1, 1/om^2];

% Shoot from each end
if mod(length(r), 2) == 0 
    i_end = length(r) - 1;   % number of radial grid points must be odd 
else
    i_end = length(r);
end
soln1 = shoot(BCin, 1, i_fit, dr, 'RHSfunc');
soln2 = shoot(BCout, i_end, i_fit, dr, 'RHSfunc');
