function [Yms, dYms, d2Yms] = construct_spherical_harmonics(l)
% Syntax: [Yms, dYms, d2Yms] = construct_spherical_harmonics(l)
% where l is the spherical harmonic degree (currently available: l=0,1,2,3
% only). Returns three structs each containing 2l+1 anonymous functions, of 
% theta and phi, for each value of m in increasing order. These correspond 
% to the spherical harmonic itself (Yms), its first theta-derivative (dYms)
% and second theta-derivative (d2Yms).
%
% E.g. if l=2, then Yms{1} corresponds to m=-2, Yms{2} to m=-1, etc.
%
% Created 29th Sep 2020            C. Loi

if l==0
    Yms = {@(th,ph) 0.5 * sqrt(1/pi)};
    dYms = {@(th,ph) 0.0};
    d2Yms = {@(th,ph) 0.0};
elseif l==1
    Yms = {@(th,ph) 0.5 * sqrt(3/(2*pi)) * sin(th) .* exp(-1i*ph), ...
        @(th,ph) 0.5 * sqrt(3/pi) * cos(th), ...
        @(th,ph) -0.5 * sqrt(3/(2*pi)) * sin(th) .* exp(1i*ph)};
    dYms = {@(th,ph) 0.5 * sqrt(3/(2*pi)) .* cos(th) .* exp(-1i*ph), ...
        @(th,ph) -0.5 * sqrt(3/pi) * sin(th), ...
        @(th,ph) -0.5 * sqrt(3/(2*pi)) * cos(th) .* exp(1i*ph)};
    d2Yms = {@(th,ph) -0.5 * sqrt(3/(2*pi)) .* sin(th) .* exp(-1i*ph), ...
        @(th,ph) -0.5 * sqrt(3/pi) * cos(th), ...
        @(th,ph) 0.5 * sqrt(3/(2*pi)) * sin(th) .* exp(1i*ph)};
elseif l==2
    Yms = {@(th,ph) 0.25 * sqrt(15/(2*pi)) * sin(th).^2 .* exp(-2i*ph), ...
        @(th,ph) 0.5 * sqrt(15/(2*pi)) * sin(th) .* cos(th) .* exp(-1i*ph), ...
        @(th,ph) 0.25 * sqrt(5/pi) * (3*cos(th).^2 - 1), ...
        @(th,ph) -0.5 * sqrt(15/(2*pi)) * sin(th) .* cos(th) .* exp(1i*ph), ...
        @(th,ph) 0.25 * sqrt(15/(2*pi)) * sin(th).^2 .* exp(2i*ph)};
    dYms = {@(th,ph) 0.25 * sqrt(15/(2*pi)) * sin(2*th) .* exp(-2i*ph), ...
        @(th,ph) 0.5 * sqrt(15/(2*pi)) * cos(2*th) .* exp(-1i*ph), ...
        @(th,ph) -0.75 * sqrt(5/pi) * sin(2*th), ...
        @(th,ph) -0.5 * sqrt(15/(2*pi)) * cos(2*th) .* exp(1i*ph), ...
        @(th,ph) 0.25 * sqrt(15/(2*pi)) * sin(2*th) .* exp(2i*ph)};
    d2Yms = {@(th,ph) 0.5 * sqrt(15/(2*pi)) * cos(2*th) .* exp(-2i*ph), ...
        @(th,ph) -sqrt(15/(2*pi)) * sin(2*th) .* exp(-1i*ph), ...
        @(th,ph) -1.5 * sqrt(5/pi) * cos(2*th), ...
        @(th,ph) sqrt(15/(2*pi)) * sin(2*th) .* exp(1i*ph), ...
        @(th,ph) 0.5 * sqrt(15/(2*pi)) * cos(2*th) .* exp(2i*ph)};
elseif l==3
    Yms = {@(th,ph) 0.125 * sqrt(35/pi) * sin(th).^3 .* exp(-3i*ph), ...
        @(th,ph) 0.25 * sqrt(105/(2*pi)) * sin(th).^2 .* cos(th) .* exp(-2i*ph), ...
        @(th,ph) 0.125 * sqrt(21/pi) .* sin(th) .* (5*cos(th).^2 - 1) .* exp(-1i*ph), ...
        @(th,ph) 0.25 * sqrt(7/pi) * cos(th) .* (5*cos(th).^2 - 3), ...
        @(th,ph) -0.125 * sqrt(21/pi) * sin(th) .* (5*cos(th).^2 - 1) .* exp(1i*ph), ...
        @(th,ph) 0.25 * sqrt(105/(2*pi)) * sin(th).^2 .* cos(th) .* exp(2i*ph), ...
        @(th,ph) -0.125 * sqrt(35/pi) * sin(th).^3 .* exp(3i*ph)};
    dYms = {@(th,ph) 0.375 * sqrt(35/pi) * sin(th).^2 .* cos(th) .* exp(-3i*ph), ...
        @(th,ph) 0.25 * sqrt(105/(2*pi)) * sin(th) .* (2*cos(th).^2 - sin(th).^2) .* exp(-2i*ph), ...
        @(th,ph) 0.125 * sqrt(21/pi) * cos(th) .* (5*cos(th).^2 - 10*sin(th).^2 - 1) .* exp(-1i*ph), ...
        @(th,ph) 0.25 * sqrt(7/pi) * sin(th) .* (3 - 15*cos(th).^2), ...
        @(th,ph) -0.125 * sqrt(21/pi) * cos(th) .* (5*cos(th).^2 - 10*sin(th).^2 - 1) .* exp(1i*ph), ...
        @(th,ph) 0.25 * sqrt(105/(2*pi)) * sin(th) .* (2*cos(th).^2 - sin(th).^2) .* exp(2i*ph), ...
        @(th,ph) -0.375 * sqrt(35/pi) * sin(th).^2 .* cos(th) .* exp(3i*ph)};
    d2Yms = {@(th,ph) 0.375 * sqrt(35/pi) * sin(th) .* (2*cos(th).^2 - sin(th).^2) .* exp(-3i*ph), ...
        @(th,ph) 0.25 * sqrt(105/(2*pi)) * cos(th) .* (2*cos(th).^2 - 7*sin(th).^2) .* exp(-2i*ph), ...
        @(th,ph) 0.125 * sqrt(21/pi) * sin(th) .* (10*sin(th).^2 + 1 - 35*cos(th).^2) .* exp(-1i*ph), ...
        @(th,ph) 0.25 * sqrt(7/pi) * cos(th) .* (3 - 15*cos(th).^2 + 30*sin(th).^2), ...
        @(th,ph) -0.125 * sqrt(21/pi) * sin(th) .* (10*sin(th).^2 + 1 - 35*cos(th).^2) .* exp(1i*ph), ...
        @(th,ph) 0.25 * sqrt(105/(2*pi)) * cos(th) .* (2*cos(th).^2 - 7*sin(th).^2) .* exp(2i*ph), ...
        @(th,ph) -0.375 * sqrt(35/pi) * sin(th) .* (2*cos(th).^2 - sin(th).^2) .* exp(3i*ph)};
elseif l==4    
    Yms = {@(th,ph) 0.1875 * sqrt(35/(2*pi)) * sin(th).^4 .* exp(-4i*ph), ...
        @(th,ph) 0.375 * sqrt(35/pi) * sin(th).^3 .* cos(th) .* exp(-3i*ph), ...
        @(th,ph) 0.375 * sqrt(5/(2*pi)) * sin(th).^2 .* (7*cos(th).^2 - 1) .* exp(-2i*ph), ...
        @(th,ph) 0.375 * sqrt(5/pi) * sin(th) .* (7*cos(th).^2 - 3) .* cos(th) .* exp(-1i*ph), ...
        @(th,ph) 0.1875 / sqrt(pi) * (35*cos(th).^4 - 30*cos(th).^2 + 3), ...
        @(th,ph) -0.375 * sqrt(5/pi) * sin(th) .* (7*cos(th).^2 - 3) .* cos(th) .* exp(1i*ph), ...
        @(th,ph) 0.375 * sqrt(5/(2*pi)) * sin(th).^2 .* (7*cos(th).^2 - 1) .* exp(2i*ph), ...
        @(th,ph) -0.375 * sqrt(35/pi) .* sin(th).^3 .* cos(th) .* exp(3i*ph), ...
        @(th,ph) 0.1875 * sqrt(35/(2*pi)) .* sin(th).^4 .* exp(4i*ph)};
    dYms = {@(th,ph) 0.75 * sqrt(35/(2*pi)) * sin(th).^3 .* cos(th) .* exp(-4i*ph), ...
        @(th,ph) 0.375 * sqrt(35/pi) * sin(th).^2 * (4*cos(th).^2 - 1) .* exp(-3i*ph), ...
        @(th,ph) 1.5 * sqrt(5/(2*pi)) * sin(th) .* cos(th) .* (3*cos(th).^2 - 4*sin(th).^2) .* exp(-2i*ph), ...
        @(th,ph) 0.375 * sqrt(5/pi) * (cos(th).^2 .* (7*cos(th).^2 - 3) - sin(th).^2 .* (21*cos(th).^2 - 3)) .* exp(-1i*ph), ...
        @(th,ph) 3.75 / sqrt(pi) * (3 - 7*cos(th).^2) .* cos(th) .* sin(th), ...
        @(th,ph) -0.375 * sqrt(5/pi) * (cos(th).^2 .* (7*cos(th).^2 - 3) - sin(th).^2 .* (21*cos(th).^2 - 3)) .* exp(1i*ph), ...
        @(th,ph) 1.5 * sqrt(5/(2*pi)) * sin(th) .* cos(th) .* (3*cos(th).^2 - 4*sin(th).^2) .* exp(2i*ph), ...
        @(th,ph) -0.375 * sqrt(35/pi) * sin(th).^2 * (4*cos(th).^2 - 1) .* exp(3i*ph), ...
        @(th,ph) 0.75 * sqrt(35/(2*pi)) * sin(th).^3 .* cos(th) .* exp(4i*ph)};
    d2Yms = {@(th,ph) 0.75 * sqrt(35/(2*pi)) * sin(th).^2 .* (3*cos(th).^2 - sin(th).^2) .* exp(-4i*ph), ...
        @(th,ph) 0.75 * sqrt(35/pi) * (3*cos(th).^2 - 5*sin(th).^2) .* sin(th) .* cos(th) .* exp(-3i*ph), ...
        @(th,ph) 1.5 * sqrt(5/(2*pi)) * (3*cos(th).^4 - 21*cos(th).^2 .* sin(th).^2 + 4*sin(th).^4) .* exp(-2i*ph), ...
        @(th,ph) -0.75 * sqrt(5/pi) * sin(th) .* cos(th) .* (35*cos(th).^2 - 6 - 21*sin(th).^2) .* exp(-1i*ph), ...
        @(th,ph) 3.75 / sqrt(pi) * (21*sin(th).^2 .* cos(th).^2 + 3*(cos(th).^2 - sin(th).^2) - 7*cos(th).^4), ...
        @(th,ph) 0.75 * sqrt(5/pi) * sin(th) .* cos(th) .* (35*cos(th).^2 - 6 - 21*sin(th).^2) .* exp(1i*ph), ...
        @(th,ph) 1.5 * sqrt(5/(2*pi)) * (3*cos(th).^4 - 21*cos(th).^2 .* sin(th).^2 + 4*sin(th).^4) .* exp(2i*ph), ...
        @(th,ph) -0.75 * sqrt(35/pi) * (3*cos(th).^2 - 5*sin(th).^2) .* sin(th) .* cos(th) .* exp(3i*ph), ...
        @(th,ph) 0.75 * sqrt(35/(2*pi)) * sin(th).^2 .* (3*cos(th).^2 - sin(th).^2) .* exp(4i*ph)};
end