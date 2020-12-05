%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Script to solve the matrix eigenvalue problem for stellar oscillations
%  perturbed by both rotation and magnetism. Each influence is treated as
%  axisymmetric, but they are in general allowed to be misaligned, 
%  controlled through the angle beta. Performs the calculation for one
%  spherical harmonic degree at a time (currently available: l=0,1,2,3),
%  outputting the frequency shifts for all possible m values, in the
%  inertial frame. The perturbative effects of rotation and magnetism are
%  included only up to first order, and distortional effects are not 
%  treated. Computes both the frequency shift and the expansion
%  coefficients of the mode of interest in terms of all other available
%  modes in the loaded file.
%
%  Requires the use of the following functions:
%  * construct_spherical_harmonics.m 
%  * Prendergast.m
%
%  Created 24th Aug 2020                 C. Loi
%
%%%%%

clear

%% USER-DEFINED PARAMETERS
   
n = -31;               % Eckart radial order of the mode we wish to consider
                       % (N.B. up to the user to ensure this is present in 
                       % the loaded file containing the eigenfunctions)
l = 1;                 % SH degree (currently available: l=0,1,2,3)
beta = pi/4;           % obliquity angle
Omega = 0.002;          % rotation frequency (fraction of dynamical freq)

% File containing the unperturbed eigenfunctions
modesInFile = ['modes_n4.2poly_l' num2str(l) '.mat'];

% Magnetic field parameters
vAcen = 1e-4;          % central Alfven speed (fraction of dynamical speed)
Rfield = 0.1;          % radial extent of the field region
min_lam = 0;           % lower search bound for lambda
max_lam = 1000;         % upper search bound
dlam = 2;            % search resolution
lamOrder = 1;          % complexity order of the field solution desired

% Turn on plotting?
wantPlots = true;

%% LOAD UNPERTURBED EIGENFUNCTION(S)

load(modesInFile)

om_n = frequencies(orders==n);
R_1d_n = R_fns(:, orders==n)';
H_1d_n = H_fns(:, orders==n)';

% Construct angular grid
nPts_r = length(r);
nPts_th = 100;
nPts_ph = 1;
dr = r(2) - r(1);
dth = pi/nPts_th;
th = 0.5*dth : dth : pi-0.5*dth;
[Mr, Mth] = meshgrid(r, th);

dR_1d_n = gradient(R_1d_n, dr);
dH_1d_n = gradient(H_1d_n, dr);
d2R_1d_n = gradient(dR_1d_n, dr);
d2H_1d_n = gradient(dH_1d_n, dr);

R_n = repmat(R_1d_n, nPts_th, 1);
H_n = repmat(H_1d_n, nPts_th, 1);
dR_n = repmat(dR_1d_n, nPts_th, 1);
dH_n = repmat(dH_1d_n, nPts_th, 1);
d2R_n = repmat(d2R_1d_n, nPts_th, 1);
d2H_n = repmat(d2H_1d_n, nPts_th, 1);

%% CHOOSE MARKER STYLES AND COLOURS (USED IF PLOTTING)

if l == 0
    mmark = 'o';           % marker styles for different m, m'
    mcol = 'k';            % colour
    mcolM = [0 0 0];       % corresponding RGB component matrix
elseif l == 1
    mmark = 'vo^';         % increasing order from -l:l
    mcol = 'bkr';
    mcolM = [0 0 1;
        0 0 0;
        1 0 0];
elseif l == 2
    mmark = 'svo^d';
    mcol = 'ybkrg';
    mcolM = [1 1 0;
        0 0 1;
        0 0 0;
        1 0 0;
        0 1 0];
elseif l == 3
    mmark = 'psvo^d*';
    mcol = 'cybkrgm';
    mcolM = [0 1 1;
        1 1 0;
        0 0 1;
        0 0 0;
        1 0 0;
        0 1 0;
        1 0 1];
end

%% COMPUTE AUXILLIARY QUANTITIES

[Yms, dYms, d2Yms] = construct_spherical_harmonics(l);

% Construct Wigner d-matrix
dmm = zeros(2*l+1);
for i = 1 : 2*l+1
    for j = 1 : 2*l+1
        m1 = i - 1 - l;
        m2 = j - 1 - l;
        
        dmm(i,j) = sqrt(factorial(l+m1) * factorial(l-m1) * ...
            factorial(l+m2) * factorial(l-m2));
        
        sfact = 0;
        for s = 0 : 2*l
            if l+m2-s < 0, continue, end
            if m1-m2+s < 0, continue, end
            if l-m1-s < 0, continue, end
            
            sfact = sfact + (-1)^(m1-m2+s) * ...
                cos(beta/2)^(2*l+m2-m1-2*s) * sin(beta/2)^(m1-m2+2*s) / ...
                factorial(l+m2-s) / factorial(s) / factorial(m1-m2+s) / ...
                factorial(l-m1-s);
        end
        
        dmm(i,j) = dmm(i,j) * sfact;
    end
end


%% MAGNETIC FIELD CONSTRUCTION

[Psis, dPsis, lambdas] = Prendergast(r, rho, Rfield, min_lam, max_lam, dlam, vAcen);
Psi_1d = Psis(lamOrder, :);
dPsi_1d = dPsis(lamOrder, :);
d2Psi_1d = gradient(dPsi_1d, dr);
d3Psi_1d = gradient(d2Psi_1d, dr);
lam = lambdas(lamOrder);

Psi = repmat(Psi_1d, nPts_th, 1);
dPsi = repmat(dPsi_1d, nPts_th, 1);
d2Psi = repmat(d2Psi_1d, nPts_th, 1);
d3Psi = repmat(d3Psi_1d, nPts_th, 1);

Br = 2*Psi .* cos(Mth) ./ Mr.^2;
Bth = -dPsi .* sin(Mth) ./ Mr;
Bph = -lam * Psi .* sin(Mth) ./ Mr;

if wantPlots
    % Visualise inputs
    figure; set(gcf, 'Position', [1750 80 1080 800])
    subplot(2,3,1), plot(r, R_1d_n, 'k'), xlabel('r'), ylabel('R')
    subplot(2,3,2), plot(r, H_1d_n, 'k'), xlabel('r'), ylabel('H')
    subplot(2,3,3), plot(r, Psi_1d, 'k'), xlabel('r'), ylabel('\Psi')
    subplot(2,3,4), imagesc(r, th, Br), xlabel('r'), ylabel('\theta'), title('B_r'), colorbar
    subplot(2,3,5), imagesc(r, th, Bth), xlabel('r'), ylabel('\theta'), title('B_\theta'), colorbar
    subplot(2,3,6), imagesc(r, th, Bph), xlabel('r'), ylabel('\theta'), title('B_\phi'), colorbar
    colormap jet
end

Ar = @(m, R, H) -(1i*m*lam * Psi .* R + 2*cos(Mth) .* (R .* dPsi - l*(l+1) * ...
    H .* Psi ./ Mr)) ./ Mr.^2 .* Yms{m+l+1}(Mth,0) + sin(Mth) ./ Mr.^2 .* (2*Psi .* H ./ Mr ...
    - dPsi .* R) .* dYms{m+l+1}(Mth,0);
Ath = @(m, R, H, dR, dH) sin(Mth) ./ Mr .* (dPsi .* dR - l*(l+1) * dPsi .* H ./ Mr + ...
    R .* d2Psi) .* Yms{m+l+1}(Mth,0) + (cos(Mth) .* (2*Psi .* dH + H .* dPsi - 2*H .* ...
    Psi ./ Mr) - 1i*m*lam * Psi .* H) ./ Mr.^2 .* dYms{m+l+1}(Mth,0) - dPsi .* H .* ...
    sin(Mth) ./ Mr.^2 .* d2Yms{m+l+1}(Mth,0);
Aph = @(m, R, H, dR, dH) (2i*m * cot(Mth) ./ Mr .* (dH .* Psi + dPsi .* H - H .* Psi ...
    ./ Mr) + m^2*lam * Psi .* H ./ Mr ./ sin(Mth) + lam * sin(Mth) .* ...
    (Psi .* dR - l*(l+1) * Psi .* H ./ Mr + R .* dPsi)) ./ Mr .* Yms{m+l+1}(Mth,0) - ...
    1i*m * dPsi .* H ./ Mr.^2 .* dYms{m+l+1}(Mth,0);

Jr = @(m, R, H, dR, dH) (Aph(m, R, H, dR, dH) .* cot(Mth) - 1i*m .* ...
    Ath(m, R, H, dR, dH) ./ sin(Mth) + (2i*m ./ (Mr .* ...
    sin(Mth)).^2 .* (H .* Psi ./ Mr - H .* dPsi - Psi .* dH) + lam*Psi ...
    .* cos(Mth) ./ Mr .* (dR - H ./ Mr .* (m^2 ./ sin(Mth).^2 + l*(l+1))) ...
    + lam*R .* dPsi .* cos(Mth) ./ Mr) .* Yms{m+l+1}(Mth,0) + (2i*m * cot(Mth) ./ Mr.^2 ...
    .* (Psi .* dH - Psi .* H ./ Mr + H .* dPsi) - lam*Psi .* H ./ Mr.^2 ...
    .* (l*(l+1) * sin(Mth) - m^2 ./ sin(Mth)) + lam*sin(Mth) ./ Mr .* ...
    (Psi .* dR + R .* dPsi)) .* dYms{m+l+1}(Mth,0) - 1i*m * H .* dPsi ...
    ./ Mr.^2 .* d2Yms{m+l+1}(Mth,0)) ./ Mr;
Jth = @(m, R, H, dR, dH, d2R, d2H) 1i*m ./ Mr ./ sin(Mth) .* Ar(m, R, H) - ...
    (2i*m * cot(Mth) ./ Mr.^2 .* (2*(dH - H ./ Mr) .* (dPsi - Psi ./ Mr) + ...
    Psi .* d2H + H .* d2Psi) - lam*sin(Mth) ./ Mr.^2 .* (H .* (dPsi - Psi ./ Mr) + Psi .* dH) ...
    .* (l*(l+1) - m^2 ./ sin(Mth).^2) + lam*sin(Mth) ./ Mr .* (Psi .* d2R ...
    + R .* d2Psi + 2*dPsi .* dR)) .* Yms{m+l+1}(Mth,0) + 1i*m ./ Mr.^2 .* (d2Psi ...
    .* H - dPsi .* (H ./ Mr - dH)) .* dYms{m+l+1}(Mth,0);
Jph = @(m, R, H, dR, dH, d2R, d2H) sin(Mth) ./ Mr .* (dPsi .* d2R + 2*d2Psi ...
    .* dR - R .* (2*dPsi ./ Mr.^2 - d3Psi) - l*(l+1) ./ Mr .* (H .* ...
    (d2Psi - dPsi ./ Mr - 2*Psi ./ Mr.^2) + dH .* dPsi)) .* Yms{m+l+1}(Mth,0) ...
    + (cos(Mth) ./ Mr.^2 .* (3*dPsi .* (dH + (R - H) ./ Mr) - 4*Psi .* dH ./ Mr + ...
    2*(1 - l*(l+1)) * Psi .* H ./ Mr.^2 + 2*Psi .* d2H + H .* d2Psi) ...
    - 1i*m*lam ./ Mr.^2 .* (H .* (dPsi - Psi ./ Mr) + Psi .* (dH - R ./ Mr) ...
    )) .* dYms{m+l+1}(Mth,0) - sin(Mth) ./ Mr.^2 .* (H .* (d2Psi - dPsi ./ Mr + 2*Psi ...
    ./ Mr.^2) + dPsi .* (dH - R ./ Mr)) .* d2Yms{m+l+1}(Mth,0);


%% CALCULATE PERTURBATIONS

disp('Computing frequency perturbations and expansion coefficients...')

% Construct matrices for rotation and magnetism (Eq 7.8)
Mmag = zeros(2*l+1); Mmag_j = zeros(2*l+1);
Mrot = zeros(2*l+1); Mrot_j = zeros(2*l+1);
dV = 2*pi * Mr.^2 .* sin(Mth) * dr * dth;

% Storage arrays for m superposition coeffs and frequency shifts
a_all = zeros(2*l+1, 2*l+1, length(orders));
om1_all = zeros(2*l+1, 2*l+1, length(orders));

% Eigenfunction expansion coefficients (Eq 7.5)
cmag = zeros(2*l+1, length(orders));
crot = zeros(2*l+1, length(orders));

for j = 1:length(orders)
    disp(['n = ' num2str(orders(j))])
    
    for i = 1:2*l+1
        m = i-l-1;
        
        xir_j = repmat(R_fns(:,j)', nPts_th, 1) .* Yms{i}(Mth,0);
        xith_j = repmat(H_fns(:,j)', nPts_th, 1) .* dYms{i}(Mth,0);
        xiph_j = repmat(H_fns(:,j)', nPts_th, 1) .* 1i*m ./ sin(Mth) .* Yms{i}(Mth,0);
        xidotxi_j = sum(rho .* r.^2 .* (R_fns(:,j).^2 + l*(l+1) * H_fns(:,j).^2)' * dr);
        
        xir_n = R_n .* Yms{i}(Mth,0);
        xith_n = H_n .* dYms{i}(Mth,0);

        R_1d_j = R_fns(:,j)';
        H_1d_j = H_fns(:,j)';
        dR_1d_j = gradient(R_1d_j, dr);
        dH_1d_j = gradient(H_1d_j, dr);
        d2R_1d_j = gradient(dR_1d_j, dr);
        d2H_1d_j = gradient(dH_1d_j, dr);

        R_j = repmat(R_1d_j, nPts_th, 1);
        H_j = repmat(H_1d_j, nPts_th, 1);
        dR_j = repmat(dR_1d_j, nPts_th, 1);
        dH_j = repmat(dH_1d_j, nPts_th, 1);
        d2R_j = repmat(d2R_1d_j, nPts_th, 1);
        d2H_j = repmat(d2H_1d_j, nPts_th, 1);
        
        Jr_j = Jr(m, R_j, H_j, dR_j, dH_j);
        Jth_j = Jth(m, R_j, H_j, dR_j, dH_j, d2R_j, d2H_j);
        Jph_j = Jph(m, R_j, H_j, dR_j, dH_j, d2R_j, d2H_j);
        
        Jr_n = Jr(m, R_n, H_n, dR_n, dH_n);
        Jth_n = Jth(m, R_n, H_n, dR_n, dH_n, d2R_n, d2H_n);
        Jph_n = Jph(m, R_n, H_n, dR_n, dH_n, d2R_n, d2H_n);
        
        % Eq 7.11, self-coupling version
        intgnd_mag_jj = (Bth .* conj(Jph_j) - Bph .* conj(Jth_j)) .* xir_j ...
            + (Bph .* conj(Jr_j) - Br .* conj(Jph_j)) .* xith_j ...
            + (Br .* conj(Jth_j) - Bth .* conj(Jr_j)) .* xiph_j;
        
        % Eq 7.12, self-coupling version
        intgnd_rot_jj = repmat(rho, nPts_th, 1) .* Omega .* ...
            imag(xiph_j .* (conj(xith_j) .* cos(Mth) + conj(xir_j) .* sin(Mth)));
        
        Mmag_j(i,i) = 0.5/om_n/xidotxi_j * sum(sum(intgnd_mag_jj .* dV));
        Mrot_j(i,i) = 2/xidotxi_j * sum(sum(intgnd_rot_jj .* dV));
        [a_j, om1_j] = eig(Mrot_j + dmm * Mmag_j * dmm');
        
        % For each eigenvector, alter sign if needed to ensure dominant
        % component is positive
        a_j = a_j .* (max(a_j) == max(abs(a_j))) ...
            - a_j .* (max(a_j) ~= max(abs(a_j)));
        
        % Make a_j is as diagonal as poss (largest entries along diagonal)
        reorder = a_j == max(a_j);
        
        if abs(det(double(reorder))) > eps
            % Only actually reorder if this is unique
            a_j = a_j * reorder;
        
            % Swap around entries of om1_j accordingly
            om1_j = diag(reorder * diag(om1_j));
        end
        
        a_all(:,:,j) = a_j;
        om1_all(:,:,j) = om1_j;
        
        if orders(j) == n
            % Set aside for plotting later
            Mmag(i,i) = Mmag_j(i,i);
            Mrot(i,i) = Mrot_j(i,i);
            a = a_j;
            om1 = om1_j;
        else
            % Eq 7.11, cross-coupling version
            intgnd_mag_nj = (Bth .* conj(Jph_n) - Bph .* conj(Jth_n)) .* xir_j ...
                + (Bph .* conj(Jr_n) - Br .* conj(Jph_n)) .* xith_j ...
                + (Br .* conj(Jth_n) - Bth .* conj(Jr_n)) .* xiph_j;
        
            % Eq 7.12, cross-coupling version
            intgnd_rot_nj = repmat(rho, nPts_th, 1) .* Omega .* ...
                imag(xiph_j .* (conj(xith_n) .* cos(Mth) + conj(xir_n) .* sin(Mth)));
        
            cmag(i,j) = sum(sum(intgnd_mag_nj .* dV)) ...
                / xidotxi_j / (om_n^2 - frequencies(j)^2);
            crot(i,j) = sum(sum(intgnd_rot_nj .* dV)) ...
                / xidotxi_j / (om_n^2 - frequencies(j)^2);
        end
    end
end

% Check size of imag. components (should be small since matrix self-adjoint)
disp(['Largest imaginary component of eigenvector: ' num2str(max(abs(imag(a(:)))))])
disp(['       "          "            eigenvalue: ' num2str(max(abs(imag(om1(:)))))])
disp('-----------------------------------------------------------')

% Inertial frame values (chosen mode)
om1_iner = repmat(real(diag(om1))', 2*l+1, 1) - repmat((-l:l)', 1, 2*l+1) * Omega;
c_iner = crot + dmm * cmag;

% Inertial frame values (all modes)
om1_all_iner = zeros(2*l+1, 2*l+1, length(orders));
for j = 1:length(orders)
    om1_all_iner(:,:,j) = repmat(real(diag(om1_all(:,:,j)))', 2*l+1, 1) ...
        - repmat((-l:l)', 1, 2*l+1) * Omega;
end

% Magnetic frame values
c = cmag + dmm' * crot;

% Reconstruct eigenfunctions
R1 = zeros(nPts_r, 2*l+1);
H1 = zeros(nPts_r, 2*l+1);
Rp = zeros(nPts_r, 2*l+1);
Hp = zeros(nPts_r, 2*l+1);

for i = 1:2*l+1
    R1(:,i) = sum(R_fns * c(i,:)', 2);
    H1(:,i) = sum(H_fns * c(i,:)', 2);
    Rp(:,i) = R_1d_n' + R1(:,i);
    Hp(:,i) = H_1d_n' + H1(:,i);
end

%% DISPLAY RESULTS

disp('Eigenvectors:')
for i = 1:2*l+1
    disp(sprintf('  %11.4g', real(a(i,:))))
end
disp('Frequency shifts (rotating frame):')
disp(sprintf('  %11.4g', real(diag(om1))'))
disp('Frequency shifts (inertial frame):')
for i = 1:2*l+1
    disp(sprintf('  %11.4g', om1_iner(i,:)))
end

if wantPlots
    load RBcmap
    
    % Plot frequency shifts in corotating frame for all modes
    figure, set(gcf, 'Position', [1800 300 1670 450])
    subplot(1,3,1)
    lgd_str = 'legend(';
    for i = 1:2*l+1
        plot(orders, reshape(real(om1_all(i,i,:)), size(orders)), [mcol(i) 'o-'])
        hold on
        lgd_str = [lgd_str '''m \approx ' num2str(i-l-1) ''', '];
    end
    xlabel('n'), ylabel('Frequency shift')
    title('Corotating frame')
    lgd_str = [lgd_str(1:end-2), ')'];
    eval(lgd_str)
    
    subplot(1,3,2)
    for i = 1:2*l+1   % index of corotating frame mode
        for j = 1:2*l+1   % index of advected component
            plot(orders, reshape(real(om1_all_iner(j,i,:)), size(orders)), ...
                [mcol(i) mmark(j) '-'])
            hold on
        end
    end
    xlabel('n'), ylabel('Frequency shift')
    title('Inertial frame')
    
    subplot(1,3,3)
    for i = 1:2*l+1   % index of corotating frame mode
        for j = 1:2*l+1   % index of advected component
            plot(orders, reshape(real(om1_all_iner(j,i,:)), size(orders)), ...
                [mcol(i) '-'])
            hold on
            for k = 1:length(orders)
                plot(orders(k), real(om1_all_iner(j,i,k)), mmark(j), ...
                    'MarkerEdgeColor', mcolM(i,:), ...
                    'MarkerFaceColor', mcolM(i,:) + (1-mcolM(i,:))*abs(1-abs(real(a_all(j,i,k)))))
            end
        end
    end
    xlabel('n'), ylabel('Frequency shift')
    title('Inertial frame')

    % Plot a coefficients for all modes
    figure, set(gcf, 'Position', [1800 100 1670 770])
    for i = 1:2*l+1   % index of corotating frame mode
        subplot(2*l+1, 1, i)
        
        imagesc(orders, -l:l, reshape(real(a_all(:,i,:)), 2*l+1, length(orders)))
        xlabel('n'), ylabel('m'), title(['Modes with m \approx ' num2str(i-l-1)])
        colorbar, colormap(RBcmap), caxis([-1 1])
        set(gca, 'YTick', -l:l, 'YTickLabel', -l:l)
    end
    
    % Visualise frequency perturbations for chosen mode
    figure; set(gcf, 'Position', [1960 -15 960 800])
    subplot(2,2,1)
    imagesc(-l:l, -l:l, real(Mmag)), xlabel('m'''), ylabel('m''')
    colorbar, colormap(RBcmap), caxis([-1 1]*max(abs(Mmag(:))))
    set(gca, 'XTick', -l:l, 'XTickLabel', -l:l, 'YTick', -l:l, 'YTickLabel', -l:l)
    title('M_{mag}')

    subplot(2,2,2)
    imagesc(-l:l, -l:l, real(Mrot)), xlabel('m'), ylabel('m')
    colorbar, caxis([-1 1]*max(abs(Mrot(:))))
    set(gca, 'XTick', -l:l, 'XTickLabel', -l:l, 'YTick', -l:l, 'YTickLabel', -l:l)
    title('M_{rot}')

    subplot(2,2,3)
    imagesc(1:2*l+1, -l:l, real(a)), xlabel('Mode index'), ylabel('m')
    colorbar, set(gca, 'XTick', 1:2*l+1, 'XTickLabel', 1:2*l+1, ...
        'YTick', -l:l, 'YTickLabel', -l:l)
    caxis([-1 1])
    title('a')

    subplot(2,2,4)
    plot([0 1], [0 0], 'k', 'LineWidth', 3)
    hold on
    for i = 1:2*l+1
        plot([3 4], [1 1]*real(om1(i,i)), 'k', 'LineWidth', 3.0)
        plot([1 3], [0 1]*real(om1(i,i)), 'k:', 'LineWidth', 1.0)

        for j = 1:2*l+1
            plot([5+i 6+i], [1 1]*om1_iner(j,i), ...
                'color', abs([1 1 1] - abs(real(a(j,i)))), ...
                'LineWidth', 3.0)
            plot([4 5+i], [real(om1(i,i)) om1_iner(j,i)], 'k:', 'LineWidth', 1.0)
            ylabel('Frequency shift')
            set(gca, 'XTick', [], 'XTickLabel', [])
        end
    end
    title(['v_{A,cen} = ' num2str(vAcen) ', \Omega = ' num2str(Omega) 10 ...
        '\beta = ' num2str(beta)])
    
    % Visualise expansion coefficients of perturbation to eigenfunction    
    lgd_str1 = 'legend(';    % wrt rotation axis (m)
    lgd_str2 = 'legend(';    % wrt magnetic axis (m')
    for i = -l:l
        lgd_str1 = [lgd_str1 '''m = ' num2str(i) ''', '];
        lgd_str2 = [lgd_str2 '''m'''' = ' num2str(i) ''', '];
    end
    lgd_str1 = [lgd_str1(1:end-2) ')'];
    lgd_str2 = [lgd_str2(1:end-2) ')'];
    
    figure; set(gcf, 'Position', [2560 120 960 800])
    subplot(2,2,1)
    for i = 1:2*l+1
        plot(orders, real(crot(i,:)), ['b' mmark(i) '-'])
        hold on
    end
    plot(n, 0, 'b.', 'MarkerSize', 14)
    xlabel('n'), eval(lgd_str1), title('c_{rot}')
    
    subplot(2,2,2)
    for i = 1:2*l+1
        plot(orders, real(cmag(i,:)), ['r' mmark(i) '-'])
        hold on
    end
    plot(n, 0, 'r.', 'MarkerSize', 14)
    xlabel('n'), eval(lgd_str2), title('c_{mag}')
    
    subplot(2,2,3)
    for i = 1:2*l+1
        plot(orders, real(c_iner(i,:)), ['k' mmark(i) '-'])
        hold on
    end
    plot(n, 0, 'k.', 'MarkerSize', 14)
    xlabel('n'), eval(lgd_str1), title('c (inertial frame)')
    
    subplot(2,2,4)
    for i = 1:2*l+1
        plot(orders, real(c(i,:)), ['k' mmark(i) '-'])
        hold on
    end
    plot(n, 0, 'k.', 'MarkerSize', 14)
    xlabel('n'), eval(lgd_str2), title('c (magnetic frame)')
    
end
