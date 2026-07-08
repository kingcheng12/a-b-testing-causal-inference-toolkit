"""
A/B Testing & Causal Inference Toolkit

Assembled from your step-by-step solutions.
"""

import numpy as np

# Step 1 - standard_normal_cdf
import math
import numpy as np

def standard_normal_cdf(z):
    # TODO: return P(Z <= z) for a standard normal Z, supporting float or numpy array input
    erf = np.vectorize(math.erf)
    return 0.5 * (1.0 + erf(z / np.sqrt(2.0)))

# Step 2 - standard_normal_ppf
import math

def standard_normal_ppf(p):
    """Return z such that Phi(z) = p for p in (0, 1)."""
    # TODO: implement a rational approximation to the inverse standard normal CDF
    if p <= 0:
        return -np.inf
    if p >= 1:
        return np.inf
    
    p_low = 0.02425
    p_high = 1 - p_low

    # Coefficients for central region
    a = [
        -3.969683028665376e+01,
         2.209460984245205e+02,
        -2.759285104469687e+02,
         1.383577518672690e+02,
        -3.066479806614716e+01,
         2.506628277459239e+00,
    ]

    b = [
        -5.447609879822406e+01,
         1.615858368580409e+02,
        -1.556989798598866e+02,
         6.680131188771972e+01,
        -1.328068155288572e+01,
    ]

    # Coefficients for tails
    c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
         4.374664141464968e+00,
         2.938163982698783e+00,
    ]

    d = [
         7.784695709041462e-03,
         3.224671290700398e-01,
         2.445134137142996e+00,
         3.754408661907416e+00,
    ]

    if p < p_low:
        # Lower tail
        q = math.sqrt(-2.0 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)

    elif p <= p_high:
        # Central region
        q = p - 0.5
        r = q * q
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)

    else:
        # Upper tail
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)

# Step 3 - pooled_proportion
def pooled_proportion(successes_a, total_a, successes_b, total_b):
    # TODO: Compute the pooled success proportion across two groups for the null of equal rates.
    return float((successes_a + successes_b) / (total_a + total_b))

# Step 4 - pooled_standard_error
import math

def pooled_standard_error(pooled_p, total_a, total_b):
    """Standard error of the difference in two proportions under the pooled null."""
    # TODO: compute sqrt( p*(1-p) * (1/n_a + 1/n_b) ) using the pooled proportion.
    var = pooled_p * (1-pooled_p) * (1/total_a + 1/total_b)
    std = math.sqrt(var)

    return std

# Step 5 - two_proportion_z_statistic
def two_proportion_z_statistic(p_a, p_b, pooled_se):
    # TODO: return the z-statistic for a two-proportion test using p_b - p_a and pooled_se
    z = (p_b - p_a) / pooled_se

    return z

# Step 6 - two_sided_p_value
def two_sided_p_value(z):
    # TODO: convert a z-statistic into a two-sided p-value under the standard normal
    p = 2 * (1-standard_normal_cdf(abs(z)))

    return p

# Step 7 - unpooled_standard_error
import math
def unpooled_standard_error(successes_a, total_a, successes_b, total_b):
    # TODO: return the unpooled SE of the difference between two sample proportions.
    p_a = successes_a/total_a
    p_b = successes_b/total_b

    se = math.sqrt(p_a*(1-p_a)/total_a + p_b*(1-p_b)/total_b)

    return se

# Step 8 - confidence_interval_from_se
def confidence_interval_from_se(point_estimate, standard_error, confidence_level):
    # TODO: build a two-sided normal-approximation CI (lower, upper) from estimate and SE
    p = 1/2 + confidence_level/2
    z = standard_normal_ppf(p)

    return point_estimate - z*standard_error, point_estimate + z*standard_error

# Step 9 - required_sample_size_per_variant (not yet solved)
# TODO: implement

# Step 10 - statistical_power (not yet solved)
# TODO: implement

# Step 11 - chi_square_statistic (not yet solved)
# TODO: implement

# Step 12 - sample_ratio_mismatch_check (not yet solved)
# TODO: implement

# Step 13 - bonferroni_correction (not yet solved)
# TODO: implement

# Step 14 - benjamini_hochberg_correction (not yet solved)
# TODO: implement

# Step 15 - group_mean_change (not yet solved)
# TODO: implement

# Step 16 - difference_in_differences_simple (not yet solved)
# TODO: implement

# Step 17 - build_did_design_matrix (not yet solved)
# TODO: implement

# Step 18 - ols_normal_equations (not yet solved)
# TODO: implement

# Step 19 - did_effect_from_regression (not yet solved)
# TODO: implement

# Step 20 - fit_synthetic_control_weights (not yet solved)
# TODO: implement

# Step 21 - synthetic_control_effect (not yet solved)
# TODO: implement

# Step 22 - ship_decision (not yet solved)
# TODO: implement

