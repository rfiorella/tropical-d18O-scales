% generate_halfdeg_testcase.m
% Generate a 0.5°x0.5° global synthetic test case for profiling.
% Grid: 720 lon x 360 lat = 259,200 points
% Run in GNU Octave: octave generate_halfdeg_testcase.m

function generate_halfdeg_testcase()

  addpath(fileparts(mfilename('fullpath')));  % for inpaint_nans

  nlat = 360;  % 0.5° spacing: 360 lat points
  nlon = 720;  % 0.5° spacing: 720 lon points

  fprintf('Generating 0.5° test case (%d lon x %d lat = %d points)...\n', ...
          nlon, nlat, nlon*nlat);

  [E, P, UQ, VQ, Tcond, LAT, LON, delta_e, Plim] = make_fields(nlat, nlon);
  LAT2 = LAT; LON2 = LON;
  dmax = 20000;

  fprintf('  Saving synthetic_halfdeg_inputs.mat...\n');
  save('-v7', '../test/fixtures/synthetic_halfdeg_inputs.mat', ...
       'E', 'P', 'UQ', 'VQ', 'Tcond', 'LAT', 'LON', 'LAT2', 'LON2', ...
       'delta_e', 'Plim', 'dmax');
  fprintf('  Done.\n');

  % Also generate the reference output
  fprintf('Running get_deltap_fast_optimized_octave...\n');

  % Alpha_eq CSV path
  csv_file = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'alpha_eq.csv');
  alpha_eq = csvread(csv_file, 1, 0);
  fprintf('  alpha_eq: %d rows, T range [%.1f, %.1f]\n', ...
          size(alpha_eq, 1), min(alpha_eq(:,1)), max(alpha_eq(:,1)));

  tic;
  [deltap_bar, tau_bar] = get_deltap_fast_optimized_octave( ...
      E, P, UQ, VQ, Tcond, LAT, LON, LAT2, LON2, ...
      delta_e, csv_file, Plim, dmax);
  elapsed = toc;

  valid = ~isnan(deltap_bar);
  fprintf('  Elapsed: %.1f seconds\n', elapsed);
  fprintf('  Valid points: %d / %d\n', sum(valid(:)), numel(deltap_bar));
  fprintf('  deltap_bar range: [%.2f, %.2f]\n', ...
          min(deltap_bar(valid)), max(deltap_bar(valid)));

  fprintf('  Saving reference_halfdeg.mat...\n');
  save('-v7', '../test/fixtures/reference_halfdeg.mat', ...
       'E', 'P', 'UQ', 'VQ', 'Tcond', 'LAT', 'LON', 'LAT2', 'LON2', ...
       'delta_e', 'alpha_eq', 'Plim', 'dmax', 'deltap_bar', 'tau_bar');
  fprintf('  Done.\n');

end


function [E, P, UQ, VQ, Tcond, LAT, LON, delta_e, Plim] = make_fields(nlat, nlon)
  % Same synthetic field generation as generate_synthetic_testcase.m
  % OUTPUT: size(field) = (nlon, nlat), rows=lon, cols=lat

  dlat = 180 / nlat;
  dlon = 360 / nlon;
  lat1d = linspace(-90 + dlat/2, 90 - dlat/2, nlat);
  lon1d = linspace(dlon/2, 360 - dlon/2, nlon);

  [LAT, LON] = meshgrid(lat1d, lon1d);

  Tsurf = 300 - 40 * (LAT / 90).^2;
  Tcond = Tsurf - 10 + 3 * sin(LON * pi / 180);

  P_itcz = 8e-5 * exp(-LAT.^2 / (2 * 8^2));
  P_midlat = 3e-5 * (exp(-(LAT - 45).^2 / (2 * 12^2)) + ...
                      exp(-(LAT + 45).^2 / (2 * 12^2)));
  P = P_itcz + P_midlat;
  P = P .* (1 + 0.3 * sin(3 * LON * pi / 180));
  P = max(P, 1e-7);

  E_sub = 5e-5 * (exp(-(LAT - 25).^2 / (2 * 15^2)) + ...
                   exp(-(LAT + 25).^2 / (2 * 15^2)));
  E_base = 2e-5 * exp(-LAT.^2 / (2 * 30^2));
  E = E_sub + E_base;
  E = E .* (1 + 0.2 * cos(2 * LON * pi / 180));
  E = max(E, 1e-8);

  UQ_westerly = 200 * sin(LAT * pi / 60);
  UQ_trade = -80 * exp(-LAT.^2 / (2 * 15^2));
  UQ = UQ_westerly + UQ_trade;
  UQ = UQ .* (1 + 0.15 * sin(2 * LON * pi / 180));

  VQ_hadley = 40 * sin(2 * LAT * pi / 180);
  VQ = VQ_hadley .* (1 + 0.2 * cos(3 * LON * pi / 180));

  Fmag = sqrt(UQ.^2 + VQ.^2);
  too_small = Fmag < 1;
  UQ(too_small) = UQ(too_small) + sign(UQ(too_small) + 0.1) * 2;
  VQ(too_small) = VQ(too_small) + sign(VQ(too_small) + 0.1) * 2;

  delta_e = -8 - 12 * (abs(LAT) / 90).^1.5;
  delta_e = delta_e + 2 * sin(4 * LON * pi / 180);
  nan_rows = max(1, min(nlon, round(nlon * [0.2, 0.5, 0.8])));
  nan_cols = max(1, min(nlat, round(nlat * [0.3, 0.6, 0.9])));
  for k = 1:length(nan_rows)
    delta_e(nan_rows(k), nan_cols(k)) = NaN;
  end

  P_sorted = sort(P(:));
  Plim_scalar = P_sorted(max(1, round(0.1 * numel(P))));
  Plim = Plim_scalar * ones(size(P));

end

generate_halfdeg_testcase();
