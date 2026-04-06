% generate_synthetic_testcase.m
% Generate synthetic ERA5-like climate fields for testing get_deltap_fast_optimized.
% Produces two test cases: tiny (6x12) and small (30x60).
% Run in GNU Octave: octave generate_synthetic_testcase.m
%
% IMPORTANT: The Matlab code expects the convention:
%   Dimension 1 (rows) = longitude
%   Dimension 2 (cols) = latitude
% This matches how griddedInterpolant is called with transposed inputs.

function generate_synthetic_testcase()

  % =====================================================================
  % Case 1: Tiny (12x6) — 12 lon x 6 lat, for fast Octave compat checking
  % =====================================================================
  fprintf('Generating tiny test case (12 lon x 6 lat)...\n');
  [E, P, UQ, VQ, Tcond, LAT, LON, delta_e, Plim] = make_fields(6, 12);
  LAT2 = LAT; LON2 = LON;
  dmax = 1500;

  save('-v7', '../test/fixtures/synthetic_tiny_inputs.mat', ...
       'E', 'P', 'UQ', 'VQ', 'Tcond', 'LAT', 'LON', 'LAT2', 'LON2', ...
       'delta_e', 'Plim', 'dmax');
  fprintf('  Saved synthetic_tiny_inputs.mat\n');

  % =====================================================================
  % Case 2: Small (60x30) — 60 lon x 30 lat, main validation case
  % =====================================================================
  fprintf('Generating small test case (60 lon x 30 lat)...\n');
  [E, P, UQ, VQ, Tcond, LAT, LON, delta_e, Plim] = make_fields(30, 60);
  LAT2 = LAT; LON2 = LON;
  dmax = 20000;

  save('-v7', '../test/fixtures/synthetic_small_inputs.mat', ...
       'E', 'P', 'UQ', 'VQ', 'Tcond', 'LAT', 'LON', 'LAT2', 'LON2', ...
       'delta_e', 'Plim', 'dmax');
  fprintf('  Saved synthetic_small_inputs.mat\n');

  fprintf('Done.\n');
end


function [E, P, UQ, VQ, Tcond, LAT, LON, delta_e, Plim] = make_fields(nlat, nlon)
% Create synthetic 2D climate fields on a regular lat/lon grid.
% OUTPUT CONVENTION: size(field) = (nlon, nlat)
%   Dimension 1 (rows) = longitude
%   Dimension 2 (cols) = latitude
% This matches the Matlab code's expected grid convention.

  % --- Grid ---
  dlat = 180 / nlat;
  dlon = 360 / nlon;
  lat1d = linspace(-90 + dlat/2, 90 - dlat/2, nlat);
  lon1d = linspace(dlon/2, 360 - dlon/2, nlon);

  % meshgrid with (lon, lat) so that rows=lon, cols=lat
  [LAT, LON] = meshgrid(lat1d, lon1d);
  % Now LAT is (nlon, nlat), LON is (nlon, nlat)
  % LAT(i,j) = lat1d(j), LON(i,j) = lon1d(i)

  % --- Surface temperature (K) ---
  Tsurf = 300 - 40 * (LAT / 90).^2;

  % --- Condensation temperature (K) ---
  Tcond = Tsurf - 10 + 3 * sin(LON * pi / 180);

  % --- Precipitation (kg/m2/s) ---
  P_itcz = 8e-5 * exp(-LAT.^2 / (2 * 8^2));
  P_midlat = 3e-5 * (exp(-(LAT - 45).^2 / (2 * 12^2)) + ...
                      exp(-(LAT + 45).^2 / (2 * 12^2)));
  P = P_itcz + P_midlat;
  P = P .* (1 + 0.3 * sin(3 * LON * pi / 180));
  P = max(P, 1e-7);

  % --- Evaporation (kg/m2/s) ---
  E_sub = 5e-5 * (exp(-(LAT - 25).^2 / (2 * 15^2)) + ...
                   exp(-(LAT + 25).^2 / (2 * 15^2)));
  E_base = 2e-5 * exp(-LAT.^2 / (2 * 30^2));
  E = E_sub + E_base;
  E = E .* (1 + 0.2 * cos(2 * LON * pi / 180));
  E = max(E, 1e-8);

  % --- Moisture transport (kg/m/s) ---
  UQ_westerly = 200 * sin(LAT * pi / 60);
  UQ_trade = -80 * exp(-LAT.^2 / (2 * 15^2));
  UQ = UQ_westerly + UQ_trade;
  UQ = UQ .* (1 + 0.15 * sin(2 * LON * pi / 180));

  VQ_hadley = 40 * sin(2 * LAT * pi / 180);
  VQ = VQ_hadley .* (1 + 0.2 * cos(3 * LON * pi / 180));

  % Ensure Fmag > 0 everywhere
  Fmag = sqrt(UQ.^2 + VQ.^2);
  too_small = Fmag < 1;
  UQ(too_small) = UQ(too_small) + sign(UQ(too_small) + 0.1) * 2;
  VQ(too_small) = VQ(too_small) + sign(VQ(too_small) + 0.1) * 2;

  % --- delta_e (permil) ---
  delta_e = -8 - 12 * (abs(LAT) / 90).^1.5;
  delta_e = delta_e + 2 * sin(4 * LON * pi / 180);
  % Sprinkle a few NaNs (to test inpaint_nans)
  nan_rows = max(1, min(nlon, round(nlon * [0.2, 0.5, 0.8])));
  nan_cols = max(1, min(nlat, round(nlat * [0.3, 0.6, 0.9])));
  for k = 1:length(nan_rows)
    delta_e(nan_rows(k), nan_cols(k)) = NaN;
  end

  % --- Plim (precipitation threshold) ---
  P_sorted = sort(P(:));
  Plim_scalar = P_sorted(max(1, round(0.1 * numel(P))));
  Plim = Plim_scalar * ones(size(P));

end

% Run if called as script
generate_synthetic_testcase();
