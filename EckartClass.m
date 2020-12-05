function n = EckartClass(r, xi, pp)
% Computes the order of a mode according to the Eckart (1960) approach,
% with inputs being the radial grid 'r', fluid displacements 'xi', and
% pressure perturbations 'pp'. These are expected to be column vectors. The
% output value n is an integer representing the order of the mode. As per
% this classification scheme, p-modes have positive order, g-modes have
% negative order, and f-modes have order zero.
%
% Note that the Eckart classification scheme is only valid for modes
% computed under the Cowling approximation.
%
% Created 11th Nov 2015           C. Loi

% Identify nodes
nodeLocs = xi(1:end-1) .* xi(2:end) < 0;

% Identify the xi and pp elements just inward (in r) of the node
if r(1) < r(end)
    xi_in = xi([nodeLocs; false]);
    pp_in = pp([nodeLocs; false]);
else
    xi_in = xi([false; nodeLocs]);
    pp_in = pp([false; nodeLocs]);
end

n = 0;
for i = 1:sum(nodeLocs)
    if pp_in(i) > 0
        if xi_in(i) > 0, n = n+1; else n = n-1; end
    else
        if xi_in(i) < 0, n = n+1; else n = n-1; end
    end
end