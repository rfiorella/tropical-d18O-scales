function A = inpaint_nans(A)
% INPAINT_NANS  Simple iterative Laplacian infill of NaN values.
% Octave-compatible replacement for John D'Errico's inpaint_nans.
% Uses iterative averaging of non-NaN neighbors until convergence.

  nan_mask = isnan(A);
  if ~any(nan_mask(:))
    return;
  end

  % Initialize NaN locations with mean of non-NaN values
  mu = mean(A(~nan_mask));
  A(nan_mask) = mu;

  [M, N] = size(A);
  max_iter = 2000;
  tol = 1e-6;

  for iter = 1:max_iter
    max_change = 0;
    for j = 1:N
      for i = 1:M
        if ~nan_mask(i, j)
          continue;
        end
        n = 0; s = 0;
        if i > 1;  n = n + 1; s = s + A(i-1, j); end
        if i < M;  n = n + 1; s = s + A(i+1, j); end
        if j > 1;  n = n + 1; s = s + A(i, j-1); end
        if j < N;  n = n + 1; s = s + A(i, j+1); end
        new_val = s / n;
        max_change = max(max_change, abs(new_val - A(i, j)));
        A(i, j) = new_val;
      end
    end
    if max_change < tol
      break;
    end
  end
end
