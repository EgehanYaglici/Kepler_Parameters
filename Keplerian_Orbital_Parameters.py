import logging
import math


class KeplerParameterHolder:
    """
    The main purpose of this class is to hold all parameters of Keplerian orbits.
    It calculates missing parameters based on the given values and maintains consistency.
    """

    def __init__(self, a=None, b=None, e=None, p=None, e_angle=None, i=None, OMEGA=None, w=None, M=None, t=None,
                 t0=None,
                 v=None, E=None,
                 unit="radians", tolerance=1e-6, strict_tolerance=1e-3):
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.tolerance = tolerance
        self.strict_tolerance = strict_tolerance
        self.unit = unit.lower()

        # Initialize primary values and apply initial calculation
        self._a = self._check_value(a, a > 0 if a else True, "Semi-major axis (a) must be positive.", "a")
        self._b = self._check_value(b, b > 0 if b else True, "Semi-minor axis (b) must be positive.", "b")
        self._e = self._check_value(e, 0 <= e < 1 if e else True, "Eccentricity (e) must be 0 ≤ e < 1.", "e")
        self._p = self._check_value(p, p > 0 if p else True, "Parameter p must be positive.", "p")
        self._e_angle = self._convert_and_check(e_angle, 0, math.pi / 2, 0, 90, "Eccentricity angle (φ)")

        # Angular parameters with conversion checks
        self._i = self._convert_and_check(i, 0, math.pi, 0, 180, "Inclination (i)")
        self._OMEGA = self._convert_and_check(OMEGA, 0, 2 * math.pi, 0, 360, "RAAN (Ω)")
        self._w = self._convert_and_check(w, 0, 2 * math.pi, 0, 360, "Argument of Perigee (ω)")

        # Anomaly parameters, allowing for degrees or radians
        self._M = self._convert_and_check(M, 0, 2 * math.pi, 0, 360, "Mean Anomaly (M)")
        self._v = self._convert_and_check(v, 0, 2 * math.pi, 0, 360, "True Anomaly (ν)")
        self._E = self._convert_and_check(E, 0, 2 * math.pi, 0, 360, "Eccentric Anomaly (E)")

        self._t = t  # Current time
        self._t0 = t0  # Reference time (epoch)

        # Initial calculation and consistency check
        self.calculate_and_verify_parameters(critical_check=True)
        self.display_parameters()

    @staticmethod
    def _check_value(value, condition, error_message, param_name):
        if not condition:
            logging.error(f"{param_name} - {error_message} Received: {value}")
            raise ValueError(f"{param_name}: {error_message}")
        if value is not None:
            logging.info(f"{param_name} set to: {value}")
        return value

    def _convert_and_check(self, angle, min_rad, max_rad, min_deg, max_deg, name):
        if angle is None:
            return None
        if self.unit == "degrees":
            if not (min_deg <= angle <= max_deg):
                logging.error(f"{name} must be between {min_deg} and {max_deg} degrees. Received: {angle}")
                raise ValueError(f"{name} must be within [{min_deg}, {max_deg}] degrees.")
            angle = math.radians(angle)
            logging.info(f"{name} converted from degrees to radians: {angle}")
        else:
            if not (min_rad <= angle <= max_rad):
                logging.error(f"{name} must be between {min_rad} and {max_rad} radians. Received: {angle}")
                raise ValueError(f"{name} must be within [{min_rad}, {max_rad}] radians.")
            logging.info(f"{name} set to {angle} radians (no conversion needed).")
        return angle

    def calculate_and_verify_parameters(self, critical_check=False):
        """Calculate missing parameters and verify consistency."""
        significant_errors = []

        # Calculate e from e_angle if provided
        if self._e_angle is not None:
            calculated_e = math.sin(self._e_angle)
            if self._e is None or not critical_check:
                self._e = calculated_e
                logging.info(f"Eccentricity (e) calculated as: {self._e} using eccentricity angle (φ).")
            elif critical_check and not math.isclose(self._e, calculated_e, rel_tol=self.strict_tolerance):
                significant_errors.append(('e', self._e, calculated_e, "e = sin(φ)"))

        # Calculate p and b from e_angle and a if provided
        if self._e_angle is not None and self._a is not None:
            cos_phi = math.cos(self._e_angle)
            calculated_p = self._a * (cos_phi ** 2)
            calculated_b = self._a * cos_phi
            if self._p is None or not critical_check:
                self._p = calculated_p
                logging.info(f"Parameter p calculated as: {self._p} using eccentricity angle (φ) and a.")
            elif critical_check and not math.isclose(self._p, calculated_p, rel_tol=self.strict_tolerance):
                significant_errors.append(('p', self._p, calculated_p, "p = a * cos^2(φ)"))

            if self._b is None or not critical_check:
                self._b = calculated_b
                logging.info(f"Semi-minor axis (b) calculated as: {self._b} using eccentricity angle (φ) and a.")
            elif critical_check and not math.isclose(self._b, calculated_b, rel_tol=self.strict_tolerance):
                significant_errors.append(('b', self._b, calculated_b, "b = a * cos(φ)"))

        # Calculate b and p from a and e if a and e are provided
        if self._a is not None and self._e is not None:
            calculated_b = self._a * math.sqrt(1 - self._e ** 2)
            calculated_p = self._a * (1 - self._e ** 2)
            if self._b is None:
                self._b = calculated_b
                logging.info(f"Semi-minor axis (b) calculated as: {self._b} using a and e.")
            elif critical_check and not math.isclose(self._b, calculated_b, rel_tol=self.strict_tolerance):
                significant_errors.append(('b', self._b, calculated_b, "b = a * sqrt(1 - e^2)"))

            if self._p is None:
                self._p = calculated_p
                logging.info(f"Parameter p calculated as: {self._p} using a and e.")
            elif critical_check and not math.isclose(self._p, calculated_p, rel_tol=self.strict_tolerance):
                significant_errors.append(('p', self._p, calculated_p, "p = a * (1 - e^2)"))

        # Calculate a and b from p and e if p and e are provided
        if self._p is not None and self._e is not None:
            calculated_a = self._p / (1 - self._e ** 2)
            calculated_b = self._p / math.sqrt(1 - self._e ** 2)
            if self._a is None:
                self._a = calculated_a
                logging.info(f"Semi-major axis (a) calculated as: {self._a} using p and e.")
            elif critical_check and not math.isclose(self._a, calculated_a, rel_tol=self.strict_tolerance):
                significant_errors.append(('a', self._a, calculated_a, "a = p / (1 - e^2)"))

            if self._b is None:
                self._b = calculated_b
                logging.info(f"Semi-minor axis (b) calculated as: {self._b} using p and e.")
            elif critical_check and not math.isclose(self._b, calculated_b, rel_tol=self.strict_tolerance):
                significant_errors.append(('b', self._b, calculated_b, "b = p / sqrt(1 - e^2)"))

        # Calculate p from a and b if a and b are provided
        if self._a is not None and self._b is not None:
            calculated_p = (self._b ** 2) / self._a
            if self._p is None or not critical_check:
                self._p = calculated_p
                logging.info(f"Parameter p recalculated as: {self._p} using a and b.")
            elif critical_check and not math.isclose(self._p, calculated_p, rel_tol=self.strict_tolerance):
                significant_errors.append(('p', self._p, calculated_p, "p = b^2 / a"))

        # Calculate e from a and b if needed and ensure a >= b
        if self._a is not None and self._b is not None:
            calculated_e = math.sqrt(1 - (self._b / self._a) ** 2)
            if self._e is None or not critical_check:
                self._e = calculated_e
                logging.info(f"Eccentricity (e) recalculated as: {self._e} using a and b.")
            elif critical_check and not math.isclose(self._e, calculated_e, rel_tol=self.strict_tolerance):
                significant_errors.append(('e', self._e, calculated_e, "e = sqrt(1 - (b / a)^2)"))
        # Calculate e_angle (φ) from e if e is provided
        if self._e is not None:
            calculated_e_angle = math.asin(self._e)
            if self._e_angle is None or not critical_check:
                self._e_angle = calculated_e_angle
                logging.info(f"Eccentricity angle (φ) calculated as: {self._e_angle} using eccentricity (e).")
            elif critical_check and not math.isclose(self._e_angle, calculated_e_angle,
                                                     rel_tol=self.strict_tolerance):
                significant_errors.append(('e_angle', self._e_angle, calculated_e_angle, "φ = arcsin(e)"))

                # Calculate e_angle (φ) from a and b if both are provided
            if self._a is not None and self._b is not None:
                calculated_e_angle = math.acos(self._b / self._a)
                if self._e_angle is None or not critical_check:
                    self._e_angle = calculated_e_angle
                    logging.info(
                        f"Eccentricity angle (φ) calculated as: {self._e_angle} using semi-major axis (a) and "
                        f"semi-minor axis (b).")
                elif critical_check and not math.isclose(self._e_angle, calculated_e_angle,
                                                         rel_tol=self.strict_tolerance):
                    significant_errors.append(('e_angle', self._e_angle, calculated_e_angle, "φ = arccos(b / a)"))

        # Calculate Mean Anomaly (M) if t, t0, and a are provided
        if self._t is not None and self._t0 is not None and self._a is not None and self._M is None:
            self.calculate_mean_anomaly()
            logging.info(f"Mean Anomaly (M) calculated as: {self._M} using mean motion and time difference (t - t0).")

        # Calculate Eccentric Anomaly (E) if M and e are provided
        if self._M is not None and self._e is not None and self._E is None:
            self.calculate_eccentric_anomaly_newton()  # or self.calculate_eccentric_anomaly_brouwer()
            logging.info(f"Eccentric Anomaly (E) calculated as: {self._E} using Mean Anomaly (M) and Eccentricity (e).")

        # Calculate True Anomaly (ν) if E and e are provided
        if self._E is not None and self._e is not None and self._v is None:
            self.calculate_true_anomaly()
            logging.info(f"True Anomaly (ν) calculated as: {self._v} using Eccentric Anomaly (E) and Eccentricity (e).")

        # Log critical errors if any are found during the initial consistency check
        if critical_check and significant_errors:
            for param, given, calculated, formula in significant_errors:
                error_message = (
                    f"CRITICAL ERROR: {param} inconsistency. Given {given}, expected {calculated:.6f} "
                    f"based on {formula}. Check initial values."
                )
                logging.error(error_message)
                print(error_message)

    @staticmethod
    def _normalize_angle(angle):
        """Normalize angle to be within [0, 2π) in radians."""
        return angle % (2 * math.pi)

    def calculate_radius(self):
        """
        Calculate the radial distance (r) from the central body.

        Returns:
            float: The radial distance (r) if calculable; otherwise None.
        """
        if self._a is not None and self._e is not None and self._E is not None:
            r = self._a * (1 - self._e * math.cos(self._E))
            logging.info(f"Radius (r) calculated using eccentric anomaly (E): {r}")
            return r
        elif self._p is not None and self._e is not None and self._v is not None:
            r = self._p / (1 + self._e * math.cos(self._v))
            logging.info(f"Radius (r) calculated using true anomaly (ν): {r}")
            return r
        elif self._a is not None and self._e is not None and self._v is not None:
            p = self._a * (1 - self._e ** 2)
            r = p / (1 + self._e * math.cos(self._v))
            logging.info(
                f"Radius (r) calculated using semi-major axis (a), eccentricity (e), and true anomaly (ν): {r}")
            return r
        logging.warning("Radius (r) cannot be calculated with the available information.")
        return None

    def calculate_position_vector(self):
        """
        Calculate the position vector (x, y) in the orbital plane.

        Returns:
            tuple: (x, y) position vector if calculable; otherwise None.
        """
        r = self.calculate_radius()
        if r is not None and self._v is not None:
            x = r * math.cos(self._v)
            y = r * math.sin(self._v)
            logging.info(f"Position vector (x, y) calculated as: ({x}, {y}) using radius (r) and true anomaly (ν).")
            return x, y
        logging.warning("Position vector (x, y) cannot be calculated with the available information.")
        return None

    def calculate_mean_anomaly(self):
        """Calculate Mean Anomaly (M) based on available data."""
        if self._E is not None and self._e is not None:
            # Calculate M using Eccentric Anomaly (E) and Eccentricity (e)
            self._M = self._normalize_angle(self._E - self._e * math.sin(self._E))
            logging.info(f"Mean Anomaly (M) calculated and normalized as: {self._M}")
            return self._M
        elif self._t is not None and self._t0 is not None and self._a is not None:
            # Mean Motion (n) and time difference (t - t0) used to calculate M
            GM = 398600.4418 * 10 ** 9  # Gravitational parameter for Earth (m^3/s^2)
            n = math.sqrt(GM / self._a ** 3)
            self._M = self._normalize_angle(n * (self._t - self._t0))
            logging.info(f"Mean Anomaly (M) calculated and normalized as: {self._M}")
            return self._M
        logging.warning("Mean Anomaly (M) cannot be calculated with the available information.")
        return None

    def calculate_eccentric_anomaly_newton(self):
        """Calculate Eccentric Anomaly (E) from Mean Anomaly (M) using Newton-Raphson method."""
        if self._M is not None and self._e is not None:
            E = self._M
            for _ in range(10):  # Newton-Raphson iteration
                E_next = E - (E - self._e * math.sin(E) - self._M) / (1 - self._e * math.cos(E))
                if abs(E_next - E) < self.tolerance:
                    break
                E = E_next
            self._E = self._normalize_angle(E)  # Normalize Eccentric Anomaly
            logging.info(f"Eccentric Anomaly (E) calculated and normalized as: {self._E}")
            return self._E
        logging.warning("Eccentric Anomaly (E) cannot be calculated with the available information.")
        return None

    def calculate_true_anomaly(self):
        """Calculate True Anomaly (ν) from Eccentric Anomaly (E) and Eccentricity (e)."""
        if self._E is not None and self._e is not None:
            # Calculate ν using the given formula
            tan_nu_over_2 = math.sqrt((1 + self._e) / (1 - self._e)) * math.tan(self._E / 2)
            self._v = self._normalize_angle(2 * math.atan(tan_nu_over_2))  # Normalize True Anomaly
            logging.info(f"True Anomaly (ν) calculated and normalized as: {self._v}")
            return self._v
        logging.warning("True Anomaly (ν) cannot be calculated with the available information.")
        return None

    def calculate_eccentric_anomaly_brouwer_hybrid(self):
        """Calculate Eccentric Anomaly (E) using an initial Brouwer-Clemence approximation refined by Newton-Raphson
        iterations."""
        if self._M is not None and self._e is not None:
            e = self._e
            M = self._M

            # Initial Brouwer-Clemence series expansion (up to the 10th order)
            E_approx = (M
                        + (e - (1 / 8) * e ** 3 + (1 / 192) * e ** 5 - (1 / 9216) * e ** 7 + (
                            1 / 368640) * e ** 9) * math.sin(M)
                        + ((1 / 2) * e ** 2 - (1 / 6) * e ** 4 + (1 / 48) * e ** 6 - (1 / 384) * e ** 8)
                        * math.sin(2 * M)
                        + ((3 / 8) * e ** 3 - (27 / 128) * e ** 5 + (243 / 5120) * e ** 7 - (
                            1701 / 40960) * e ** 9) * math.sin(3 * M)
                        + ((1 / 3) * e ** 4 - (4 / 15) * e ** 6 + (4 / 105) * e ** 8) * math.sin(4 * M)
                        + ((125 / 384) * e ** 5 - (3125 / 9216) * e ** 7 + (16875 / 286720) * e ** 9) * math.sin(5 * M)
                        + (27 / 80) * e ** 6 * math.sin(6 * M)
                        - (16807 / 46080) * e ** 7 * math.sin(7 * M)
                        + (531441 / 1161216) * e ** 8 * math.sin(8 * M)
                        - (45927 / 40960) * e ** 9 * math.sin(9 * M)
                        + (8192 / 1474560) * e ** 10 * math.sin(10 * M)
                        )

            # Refine using Newton-Raphson for better accuracy
            E = E_approx
            for _ in range(10):  # Limit to a few iterations to refine accuracy without excessive computation
                E_next = E - (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
                if abs(E_next - E) < self.tolerance:
                    break
                E = E_next

            # Normalize and store the refined Eccentric Anomaly
            self._E = self._normalize_angle(E)
            logging.info(
                f"Eccentric Anomaly (E) calculated and refined as: {self._E} using hybrid series and Newton-Raphson.")
            return self._E

        logging.warning("Eccentric Anomaly (E) cannot be calculated with the available information.")
        return None

    def compare_eccentric_anomaly_methods(self):
        """Compare Eccentric Anomaly results from Newton-Raphson and Brouwer-Clemence methods."""
        E_newton = self.calculate_eccentric_anomaly_newton()
        E_brouwer = self.calculate_eccentric_anomaly_brouwer_hybrid()

        print(f"Eccentric Anomaly (Newton-Raphson method): {E_newton}")
        print(f"Eccentric Anomaly (Brouwer-Clemence series, hybrid via Newton-Raphson): {E_brouwer}")
        if E_newton is not None and E_brouwer is not None:
            difference = abs(E_newton - E_brouwer)
            print(f"Difference between methods: {difference}")
        else:
            print("One or both methods failed to compute Eccentric Anomaly.")

    def display_parameters(self):
        """Display all current parameters of the orbit for user verification in radians and degrees."""
        parameters = {
            "a": self._a,
            "b": self._b,
            "e": self._e,
            "p": self._p,
            "e_angle": self._e_angle,
            "i": self._i,
            "Ω": self._OMEGA,
            "ω": self._w,
            "M": self._M,
            "t": self._t,
            "t0": self._t0,
            "ν": self._v,
            "E": self._E,
        }

        print("\n--- Current Orbital Parameters ---")
        for param, value in parameters.items():
            if value is not None:
                if param in ["e_angle", "i", "Ω", "ω", "M", "ν", "E"]:
                    # Convert radians to degrees
                    value_deg = math.degrees(value)
                    print(f"{param} = {value:.6f} radians | {value_deg:.6f} degrees")
                else:
                    print(f"{param} = {value}")
            else:
                print(f"{param} = None")
        print("---------------------------------\n")

    def __str__(self):
        """Provide a string representation for easier debugging and direct print statements."""
        return (f"KeplerOrbit(a={self._a}, b={self._b}, e={self._e}, p={self._p}, "
                f"e_angle={self._e_angle}, i={self._i}, Ω={self._OMEGA}, ω={self._w}, "
                f"M={self._M}, t={self._t}, t0={self._t0}, ν={self._v}, E={self._E})")


orbit = KeplerParameterHolder(a=7000000,b=234134,M=85,unit='degrees')

print(orbit)