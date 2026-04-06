function [deltap_bar, tau_bar] = get_deltap_fast_optimized(E, P, UQ, VQ, Tcond, ...
    LAT, LON, LAT2, LON2, delta_e, csv_file, Plim, dmax)
% get_deltap_fast_optimized
% E: 2D evaporation (kg/m2/s)
% P: 2D precipitation (kg/m2/s)
% UQ: 2D zonal IVT (kg/m/s)
% VQ: 2d meridional IVT (kg/m/s)
% Tcond: 2d mean condensation temperature (K)
% delta_e: 2d delta 18O of evaporative flux
% csv file: name of csv file containing 2 columns: 
%   (1) condensation temperature (in C)
%   (2) effective fractionation factors (from Siler et al., 2021)
% Plim: 2d precipitation cutoff values, below which calculation is skipped
% dmax: scalar maximum value (in km) to integrate along upstream path
% LAT/LON: regular 2D grids of latitude and longitude, with dimensions
%   equal to the those of all other input variables
% LAT2/LON2: regular 2D grids of latitude and longitude points at which
%   deltap_bar and tau_bar will be calculated
% 

%---------------------------------------
% Setup constants and grid extensions
%---------------------------------------
dx  = 15;              % km
Dx  = dx / 111;        % deg latitude increment
Nmax = ceil(dmax / dx);

alpha_eq=load(csv_file); % csv file should contain columns of (1) T (in C) 
% and (2) equilibrium alpha

% Clamp delta_e
delta_e(delta_e < -200) = -200;
delta_e(delta_e > 200)  = 200;

% Remove nans from delta_e and Tcond
if sum(isnan(delta_e(:)))>0
    delta_e=single(inpaint_nans(double(delta_e)));
end
if sum(isnan(Tcond(:)))>0
    Tcond=single(inpaint_nans(double(Tcond)));
end

% Expand grid cyclically
LON = cat(1, LON(end,:) - 360, LON, LON(1,:) + 360);
LAT = cat(1, LAT(end,:), LAT, LAT(1,:));
E   = cat(1, E(end,:), E, E(1,:));
P   = cat(1, P(end,:), P, P(1,:));
UQ  = cat(1, UQ(end,:), UQ, UQ(1,:));
VQ  = cat(1, VQ(end,:), VQ, VQ(1,:));
Tcond = cat(1, Tcond(end,:), Tcond, Tcond(1,:));
Plim  = cat(1, Plim(end,:), Plim, Plim(1,:));
delta_e = cat(1, delta_e(end,:), delta_e, delta_e(1,:));

%---------------------------------------
% Precompute interpolants
%---------------------------------------
Pfit       = griddedInterpolant(LAT', LON', P', 'linear');
vqfit      = griddedInterpolant(LAT', LON', VQ', 'linear');
uqfit      = griddedInterpolant(LAT', LON', UQ', 'linear');
Efit       = griddedInterpolant(LAT', LON', E', 'linear');
Tcondfit   = griddedInterpolant(LAT', LON', Tcond', 'linear');
deltaefit  = griddedInterpolant(LAT', LON', delta_e', 'linear');
Plimfit    = griddedInterpolant(LAT', LON', Plim', 'linear');

%---------------------------------------
% Initialize outputs
%---------------------------------------
[A, B] = size(LAT2);
deltap_bar = nan(A, B);
tau_bar    = nan(A, B);

%---------------------------------------
% Precompute P2 and Plim2 for validity
%---------------------------------------
P2    = Pfit(LAT2', LON2')';
Plim2 = Plimfit(LAT2', LON2')';

%---------------------------------------
% Main grid loops
%---------------------------------------
for i = 1:A
    for j = 1:B

        % Skip invalid points
        if P2(i,j) < Plim2(i,j)
            continue
        end

        %---------------------------------------
        % Initialize streamline integration
        %---------------------------------------
        lat0 = LAT2(i,j);
        lon0 = LON2(i,j);

        tau = zeros(Nmax, 1);
        E0  = zeros(Nmax, 1);
        Tcond0 = zeros(Nmax, 1);
        Fmag0  = zeros(Nmax, 1);
        P0  = zeros(Nmax, 1);
        delta_e0 = zeros(Nmax, 1);

        vq0 = vqfit(lat0, lon0);
        uq0 = uqfit(lat0, lon0);
        E0(1)       = Efit(lat0, lon0);
        Tcond0(1)   = Tcondfit(lat0, lon0);
        Fmag0(1)    = sqrt(uq0.^2 + vq0.^2);
        P0(1)       = Pfit(lat0, lon0);
        delta_e0(1) = deltaefit(lat0, lon0);

        jj = 1;
        while tau(jj) < 10 && jj < Nmax-1
            % Streamline step
            dtheta = -vq0 ./ Fmag0(jj) .* Dx;
            dphi   = -uq0 ./ Fmag0(jj) .* Dx ./ cosd(lat0 + dtheta / 2);

            lat1 = lat0 + dtheta;
            lon1 = lon0 + dphi;

            % Latitude/longitude wrapping
            if lat1 < -90
                lat1 = 180 + lat1; lon1 = lon1 - 180;
            elseif lat1 > 90
                lat1 = 180 - lat1; lon1 = lon1 - 180;
            end
            if lon1 < 0
                lon1 = lon1 + 360;
            elseif lon1 > 360
                lon1 = lon1 - 360;
            end

            % Interpolate fields at new point
            vq1 = vqfit(lat1, lon1);
            uq1 = uqfit(lat1, lon1);
            E0(jj+1)       = Efit(lat1, lon1);
            Tcond0(jj+1)   = Tcondfit(lat1, lon1);
            Fmag0(jj+1)    = sqrt(uq1.^2+vq1.^2);
            P0(jj+1)       = Pfit(lat1, lon1);
            delta_e0(jj+1) = deltaefit(lat1, lon1);

            % Advance
            vq0 = vq1; uq0 = uq1;
            lat0 = lat1; lon0 = lon1;

            jj = jj + 1;
            mu = P0(jj-1) / Fmag0(jj-1);
            tau(jj) = tau(jj-1) + mu * dx * 1000;
        end

        %---------------------------------------
        % Compute along-streamline properties
        %---------------------------------------
        mu = P0 ./ Fmag0;
        tau = cumsum(mu) .* dx * 1000;

        valid_idx = (tau > 0 & delta_e0 < 120 & delta_e0 > -120);
        valid_idx(1) = true;

        alphac0 = interp1(alpha_eq(:,1),alpha_eq(:,2),Tcond0 - 273.15,...
            'linear','extrap');
        epsilonc0 = (alphac0 - 1) * 1000;

        % Weight function
        weight_core = mu(valid_idx);
        decay = exp(-flip(cumsum(flip(weight_core) * dx * 1000)));
        weight = cat(1, weight_core .* decay, zeros(sum(~valid_idx), 1));

        % Mean alpha and epsilon
        alpha_bar0 = cumsum(weight .* alphac0) ./ cumsum(weight);
        alpha_bar0(isnan(alpha_bar0)) = 1.01;
        epsilon_bar0 = (alpha_bar0 - 1) * 1000;

        % delta_p and exponential weighting
        delta_p0 = -tau .* epsilon_bar0 + epsilonc0(1) + delta_e0;
        exp_neg_tau = exp(-tau(valid_idx));
        eweight = E0(valid_idx) .* exp_neg_tau;

        num_deltap = sum(delta_p0(valid_idx) .* eweight);
        num_tau    = sum(tau(valid_idx) .* eweight);
        den        = sum(eweight);

        if den > 0
            deltap_bar(i,j) = num_deltap / den;
            tau_bar(i,j)    = num_tau / den;
        end
    end
end
end

