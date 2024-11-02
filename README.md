
# KeplerParameterHolder Class Explanation

The `KeplerParameterHolder` class is designed to manage and maintain consistency among parameters of a Keplerian orbit. 
This class can calculate missing parameters based on given values, ensuring accurate orbital properties. 
Below is a detailed explanation of its structure, methods, and usage.

## 1. Class Overview
The `KeplerParameterHolder` class is primarily responsible for:
- Holding various orbital parameters (semi-major axis, eccentricity, anomalies, etc.).
- Calculating missing parameters based on input.
- Maintaining consistency between values and checking for any parameter errors.

### Constructor Parameters
- **a, b**: Semi-major and semi-minor axes (must be positive).
- **e**: Eccentricity (0 ≤ e < 1).
- **p**: Parameter for orbit calculation (positive).
- **e_angle**: Eccentricity angle (φ), converted based on units.
- **i**: Inclination (0 to π radians).
- **OMEGA, w**: RAAN (Ω) and Argument of Perigee (ω).
- **M**: Mean anomaly.
- **t, t0**: Current time and epoch time.
- **v, E**: True and Eccentric anomalies.
- **unit**: Angle unit (either `radians` or `degrees`).
- **tolerance, strict_tolerance**: Tolerances for calculations and checks.

## 2. Key Methods

### `_check_value()`
This helper function validates a parameter value based on given conditions. If the value does not meet the conditions, 
it logs an error message.

### `_convert_and_check()`
Converts angles to radians if given in degrees and verifies they fall within an expected range.

### `_normalize_angle()`
Normalizes an angle to the range `[0, 2π)` radians.

### `calculate_and_verify_parameters()`
Calculates any missing parameters based on available values, ensuring consistency. This includes calculating:
- **Eccentricity (e)** from **eccentricity angle (φ)** if only `φ` is given.
- **Semi-minor axis (b)** and parameter `p` from **eccentricity angle (φ)** and **semi-major axis (a)**.
- **Parameter (p)** from **semi-major axis (a)** and **eccentricity (e)**.
- **Eccentricity angle (φ)** from eccentricity **e** or both **a** and **b**.

Additionally, it calculates:
- **Mean Anomaly (M)** if `t`, `t0`, and `a` are provided.
- **Eccentric Anomaly (E)** from **Mean Anomaly (M)** and **Eccentricity (e)**.
- **True Anomaly (v)** from **Eccentric Anomaly (E)** and **Eccentricity (e)**.

### `calculate_radius()`
Calculates the radial distance (`r`) from the central body. The method uses parameters such as `a`, `e`, and `E`, or 
alternatively, `p`, `e`, and `v`.

### `calculate_position_vector()`
Calculates the `(x, y)` position vector in the orbital plane based on radius (`r`) and true anomaly (`v`).

### `calculate_mean_anomaly()`
Calculates Mean Anomaly (`M`) using either Eccentric Anomaly (`E`) or mean motion and time difference.

### `calculate_eccentric_anomaly_newton()`
Calculates Eccentric Anomaly (`E`) using the Newton-Raphson method, refining the result based on tolerance.

### `calculate_true_anomaly()`
Calculates True Anomaly (`v`) from **Eccentric Anomaly (E)** and **Eccentricity (e)**.

### `calculate_eccentric_anomaly_brouwer_hybrid()`
Uses a combination of Brouwer-Clemence series expansion and Newton-Raphson iteration for a refined calculation of `E`.

### `compare_eccentric_anomaly_methods()`
Compares Eccentric Anomaly results from Newton-Raphson and Brouwer-Clemence methods to determine any differences.

## 3. Display and Debugging Methods

### `display_parameters()`
Logs and displays all current orbital parameters, showing values in both radians and degrees where applicable.

### `__str__()`
Returns a string representation of the object, useful for debugging and easy display.

## Example Usage

### Example 1
```python
orbit = KeplerParameterHolder(a=7000000, b=234134)
```

Output:
```
--- Current Orbital Parameters ---
a = 7000000
b = 234134
e = 0.9994404686668742
p = 7831.247136571428
e_angle = 1.537342 radians | 88.083230 degrees
i = None
Ω = None
ω = None
M = None
t = None
t0 = None
ν = None
E = None
---------------------------------
```

This example shows the class calculating additional parameters automatically based on the input values for `a` and `b`.


### Example 2
```python
orbit = KeplerParameterHolder(a=7000000, b=234134, t=5, t0=2)
```

Output:
```
--- Current Orbital Parameters ---
a = 7000000
b = 234134
e = 0.9994404686668742
p = 7831.247136571428
e_angle = 1.537342 radians | 88.083230 degrees
i = None
Ω = None
ω = None
M = 0.003234 radians | 0.185296 degrees
t = 5
t0 = 2
ν = 3.283759 radians | 188.145551 degrees
E = 5.821674 radians | 333.557365 degrees
---------------------------------
```

In this example, additional time parameters `t` and `t0` enable the calculation of `M`, `ν`, and `E` based on the time difference.

### Example 3
```python
orbit = KeplerParameterHolder(a=7000000, b=234134, M=85, unit='degrees')
```

Output:
```
--- Current Orbital Parameters ---
a = 7000000
b = 234134
e = 0.9994404686668742
p = 7831.247136571428
e_angle = 1.537342 radians | 88.083230 degrees
i = None
Ω = None
ω = None
M = 1.483530 radians | 85.000000 degrees
t = None
t0 = None
ν = 3.125745 radians | 179.091996 degrees
E = 2.256849 radians | 129.307942 degrees
---------------------------------
```

Here, `M` is provided in degrees, and the class correctly converts and uses this value, calculating `ν` and `E` accordingly.

---
