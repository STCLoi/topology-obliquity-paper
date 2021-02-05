function soln = shoot(BC, i_start, i_end, dr, RHSfunc)
% Integrate vector of dependent variables from i_start to i_end using the
% shooting method. To be used with get_inout_solns.m. The input 
% vector 'BC' needs to be a row vector containing starting boundary 
% conditions for each dependent variable. Returns the Nr x Nv array xGrid 
% where Nr is the number of radial gridpoints (excluding halfway points) 
% and Nv is the number of dependent variables.
%
% Created 30th Sep 2020                 C. Loi

global cRH cRR cHR cHH

soln = zeros(abs(i_end - i_start)/2 + 1, length(BC));

% Work out if we are integrating forwards or backwards
backwd = i_end < i_start;

% Set up indices and increments accordingly
j = size(soln,1)^backwd;
soln(j,:) = BC;
pm1 = (-1)^backwd;

for i = i_start : 2*pm1 : i_end-2*pm1
    soln(j+pm1,:) = RK4(soln(j,:), pm1, i, pm1*dr, RHSfunc);
    j = j + pm1;
end


function x_new = RK4(x_old, pm1, i, dr, RHSfunc)
% Performs one step of RK4 to evaluate the new dependent variables at the
% next radial grid point. The flag 'pm1' indicates whether the integration
% is forwards or backwards in r, and dr specifies half the radial distance
% between start and end points. Note that the background quantities are
% defined at not just these points but also the halfway points. The RHS
% function name should be passed as the string 'RHSfunc'.

global cRH cRR cHR cHH

j = i + pm1;       % index of halfway point
k = i + 2*pm1;     % index of endpoint

f1 = feval(RHSfunc, x_old, cRH(i), cRR(i), cHR(i), cHH(i));
f2 = feval(RHSfunc, x_old + dr*f1, cRH(j), cRR(j), cHR(j), cHH(j));
f3 = feval(RHSfunc, x_old + dr*f2, cRH(j), cRR(j), cHR(j), cHH(j));
f4 = feval(RHSfunc, x_old + 2*dr*f3, cRH(k), cRR(k), cHR(k), cHH(k));

x_new = x_old + dr * (f1 + 2*f2 + 2*f3 + f4) / 3;


function f = RHSfunc(x, c1, c2, c3, c4)
% Evaluates the RHS function in Eq 7.36 of my notes

R = x(1);
H = x(2);

f = zeros(size(x));

f(1) = c1 * H + c2 * R;
f(2) = c3 * R + c4 * H;
