% generate_reference.m
% Generate reference input/output files for Julia port validation.
% Run in GNU Octave: octave generate_reference.m
%
% Produces .mat files in v7 format (HDF5-compatible, readable by Julia's MAT.jl).

function generate_reference()

  addpath('.');

  csv_file = '../data/alpha_eq.csv';
  alpha_eq = load(csv_file);

  % =====================================================================
  % Case 1: Tiny (12 lon x 6 lat)
  % =====================================================================
  fprintf('=== Case 1: Tiny (12x6) ===\n');
  data = load('../test/fixtures/synthetic_tiny_inputs.mat');

  tic;
  [deltap_bar, tau_bar] = get_deltap_fast_optimized_octave( ...
      data.E, data.P, data.UQ, data.VQ, data.Tcond, ...
      data.LAT, data.LON, data.LAT2, data.LON2, ...
      data.delta_e, csv_file, data.Plim, data.dmax);
  elapsed = toc;

  E = data.E; P = data.P; UQ = data.UQ; VQ = data.VQ;
  Tcond = data.Tcond; LAT = data.LAT; LON = data.LON;
  LAT2 = data.LAT2; LON2 = data.LON2;
  delta_e = data.delta_e; Plim = data.Plim; dmax = data.dmax;

  save('-v7', '../test/fixtures/reference_tiny.mat', ...
       'E', 'P', 'UQ', 'VQ', 'Tcond', 'LAT', 'LON', 'LAT2', 'LON2', ...
       'delta_e', 'Plim', 'dmax', 'alpha_eq', ...
       'deltap_bar', 'tau_bar', 'elapsed');

  valid = ~isnan(deltap_bar(:));
  fprintf('  Elapsed: %.2f seconds\n', elapsed);
  fprintf('  Non-NaN: %d / %d\n', sum(valid), numel(deltap_bar));
  fprintf('  deltap_bar: [%.4f, %.4f]\n', min(deltap_bar(valid)), max(deltap_bar(valid)));
  fprintf('  tau_bar:    [%.6f, %.6f]\n', min(tau_bar(valid)), max(tau_bar(valid)));
  fprintf('  Saved reference_tiny.mat\n\n');

  % =====================================================================
  % Case 2: Small (60 lon x 30 lat)
  % =====================================================================
  fprintf('=== Case 2: Small (60x30) ===\n');
  data = load('../test/fixtures/synthetic_small_inputs.mat');

  tic;
  [deltap_bar, tau_bar] = get_deltap_fast_optimized_octave( ...
      data.E, data.P, data.UQ, data.VQ, data.Tcond, ...
      data.LAT, data.LON, data.LAT2, data.LON2, ...
      data.delta_e, csv_file, data.Plim, data.dmax);
  elapsed = toc;

  E = data.E; P = data.P; UQ = data.UQ; VQ = data.VQ;
  Tcond = data.Tcond; LAT = data.LAT; LON = data.LON;
  LAT2 = data.LAT2; LON2 = data.LON2;
  delta_e = data.delta_e; Plim = data.Plim; dmax = data.dmax;

  save('-v7', '../test/fixtures/reference_small.mat', ...
       'E', 'P', 'UQ', 'VQ', 'Tcond', 'LAT', 'LON', 'LAT2', 'LON2', ...
       'delta_e', 'Plim', 'dmax', 'alpha_eq', ...
       'deltap_bar', 'tau_bar', 'elapsed');

  valid = ~isnan(deltap_bar(:));
  fprintf('  Elapsed: %.2f seconds\n', elapsed);
  fprintf('  Non-NaN: %d / %d\n', sum(valid), numel(deltap_bar));
  fprintf('  deltap_bar: [%.4f, %.4f]\n', min(deltap_bar(valid)), max(deltap_bar(valid)));
  fprintf('  tau_bar:    [%.6f, %.6f]\n', min(tau_bar(valid)), max(tau_bar(valid)));
  fprintf('  Saved reference_small.mat\n\n');

  fprintf('Done. Reference files in test/fixtures/\n');
end

generate_reference();
