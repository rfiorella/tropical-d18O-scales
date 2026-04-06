"""
    inpaint_nans!(A::Matrix{T}; max_iter=2000, tol=1e-6) where T

Simple iterative Laplacian infill of NaN values in a 2D matrix.
Replaces NaN entries with the average of their non-NaN neighbors,
iterating until convergence. Modifies A in-place.
"""
function inpaint_nans!(A::Matrix{T}; max_iter=2000, tol=1e-6) where T
    nan_mask = isnan.(A)
    any(nan_mask) || return A

    # Initialize NaN locations with mean of non-NaN values
    valid_vals = filter(!isnan, A)
    μ = isempty(valid_vals) ? zero(T) : T(mean(valid_vals))
    A[nan_mask] .= μ

    for iter in 1:max_iter
        max_change = zero(T)
        for j in axes(A, 2), i in axes(A, 1)
            nan_mask[i, j] || continue
            n = 0; s = zero(T)
            if i > 1;            n += 1; s += A[i-1, j]; end
            if i < size(A, 1);   n += 1; s += A[i+1, j]; end
            if j > 1;            n += 1; s += A[i, j-1]; end
            if j < size(A, 2);   n += 1; s += A[i, j+1]; end
            new_val = s / n
            max_change = max(max_change, abs(new_val - A[i, j]))
            A[i, j] = new_val
        end
        max_change < tol && break
    end
    return A
end
