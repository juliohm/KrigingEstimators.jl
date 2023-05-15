# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    kriging_ui(domain, variogram, mean, degree, drifts)

Return the appropriate Kriging estimator for the `domain` based on
end-user inputs such as `variogram`, `mean`, `degree` and `drifts`.
"""
function kriging_ui(domain, variogram, mean, degree, drifts)
  if drifts ≠ nothing
    ExternalDriftKriging(variogram, drifts)
  elseif degree ≠ nothing
    UniversalKriging(variogram, degree, embeddim(domain))
  elseif mean ≠ nothing
    SimpleKriging(variogram, mean)
  else
    OrdinaryKriging(variogram)
  end
end
