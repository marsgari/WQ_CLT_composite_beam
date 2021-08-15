from math import sqrt, cosh, sinh, pi
from wq_clt_input import *


class WQBeam:
    def __init__(self):
        self.b_t = TOP_FLANGE[0]  # top flange width [mm]
        self.t_t = TOP_FLANGE[1]  # top flange thickness [mm]
        self.b_b = BOTTOM_FLANGE[0]  # bottom flange width [mm]
        self.t_b = BOTTOM_FLANGE[1]  # bottom flange thickness [mm]
        self.h_b = WALL[0]  # wall height [mm]
        self.t_w = WALL[1]  # wall thickness [mm]
        self.E_1 = STEEL_YOUNG  # steel Young's modulus [MPa]
        self.rho = STEEL_DENSITY  # steel density [kg/m^3]

        # Cross-sectional area
        self.A_1 = self.b_t * self.t_t + self.b_b * self.t_b + 2 * self.h_b * self.t_w

        # Static moment in relation to the bottom fiber [mm^3]
        self.S_bot = self.b_t * self.t_t * (self.t_b + self.h_b - self.t_t / 2) \
                     + 2 * self.h_b * self.t_w * (self.t_b + self.h_b / 2) + self.b_b * self.t_b * self.t_b / 2

        # Distances from the mass center to the top and bottom fibers [mm]
        self.y_bot = self.S_bot / self.A_1
        self.y_top = self.t_b + self.h_b - self.y_bot

        # Second moment of area [mm^4]
        self.I_1 = self.b_t * self.t_t ** 3 / 12 + self.b_t * self.t_t * (self.y_top - self.t_t / 2) ** 2 \
                   + self.b_b * self.t_b ** 3 / 12 + self.b_b * self.t_b * (self.y_bot - self.t_b / 2) ** 2 \
                   + 2 * (self.t_w * self.h_b ** 3 / 12 + self.t_w * self.h_b * (
                self.h_b / 2 + self.t_b - self.y_bot) ** 2)

        # Section moduli [mm^3]
        self.W_1top = - self.I_1 / self.y_top
        self.W_1bot = self.I_1 / self.y_bot


class CLT:
    def __init__(self):
        self.thicknesses = LAYERS  # thicknesses of layers from bottom to top [mm]
        self.orientation = ORIENTATION  # orientation of layers from bottom to top (0 if perpendicular to beam)
        self.h_l = sum(self.thicknesses)  # total thickness of CLT slab [mm]
        self.rho = CLT_DENSITY  # wood density [kg/m3]
        self.E_0 = CLT_YOUNG[0]  # Young's modulus along grains [MPa]
        self.E_90 = CLT_YOUNG[1]  # Young's modulus transverse to grains [MPa]
        self.G_9090 = ROLLING_MODULUS  # rolling shear modulus [MPa]
        self.gamma_used = GAMMA  # change to True if CLT bending stiffness is calculated with gamma-factors

        # Properties of the slab layer-by-layer
        self.E_long = []  # Young's moduli in longitudinal direction [MPa]
        self.E_tran = []  # Young's moduli in transverse direction [MPa]
        self.A = []  # areas [mm2]
        self.J = []  # second moments of area in relation to layer centroid [mm4]
        self.zeta = []  # distances from layer centroid to CLT centroid [mm]
        self.gamma_long = [1] * len(self.thicknesses)  # gamma-factors in longitudinal direction
        self.gamma_tran = [1] * len(self.thicknesses)  # gamma-factors in transverse direction

        # Properties of homogenized slab
        self.A_2 = 0  # area [mm2]
        self.I_2 = 0  # second moment of area [mm4]
        self.E_2_long = 0  # Young's modulus in longitudinal direction [MPa]
        self.E_2_tran = 0  # Young's modulus in transverse direction [MPa]
        self.EI_eff_2_long = 0  # bending stiffness in longitudinal direction
        self.EI_eff_2_tran = 0  # bending stiffness in transverse direction

    def calculate_properties(self, B_eff, span):
        """
        Calculate the properties of the CLT slab.
        Some properties depend of the effective width and the span of the WQ-beam.
        The properties are calculated REGARDLESS the orientation of the CLT slab in relation to the WQ-beam.
        The longitudinal (index _long) corresponds to the direction of the span of the slab.
        Transverse (index _tran) direction is perpendicular to longitudinal.
        """
        # Get lists with A, J zeta and E layer-by-layer
        self._make_lists(B_eff)
        self._clt_young_moduli()

        # Gamma factors
        # If the gamma method is not used, the default value gamma = 1 is used (see class constructor)
        if self.gamma_used:
            self._calculate_gamma(span)

        # EI_eff [N*mm^2]
        for r in range(0, len(self.thicknesses)):
            self.EI_eff_2_long += self.E_long[r] * self.J[r] + self.gamma_long[r] * self.zeta[r] ** 2 \
                                  * self.E_long[r] * self.A[r]
            self.EI_eff_2_tran += self.E_tran[r] * self.J[r] + self.gamma_tran[r] * self.zeta[r] ** 2 \
                                  * self.E_tran[r] * self.A[r]

        # CLT properties
        self.A_2 = B_eff * self.h_l  # [mm^2]
        self.I_2 = B_eff * self.h_l ** 3 / 12  # [mm^4]
        self.E_2_long = self.EI_eff_2_long / self.I_2  # [N/mm^2] = [MPa]
        self.E_2_tran = self.EI_eff_2_tran / self.I_2  # [N/mm^2] = [MPa]

    def _make_lists(self, B_eff):
        """Create the list with corresponding properties layer-by-layer starting from bottom"""
        layer_bottom = -self.h_l / 2  # distance to bottom of CLT slab [mm]
        for t in self.thicknesses:
            self.A.append(B_eff * t)
            self.J.append(B_eff * t ** 3 / 12)
            self.zeta.append(layer_bottom + t / 2)
            layer_bottom += t

    def _clt_young_moduli(self):
        """Create the list with the Young's moduli layer-by-layer starting from bottom"""
        for angle in self.orientation:
            if angle == 0:
                self.E_long.append(self.E_0)
                self.E_tran.append(self.E_90)
            else:
                self.E_long.append(self.E_90)
                self.E_tran.append(self.E_0)

    def _calculate_gamma(self, span):
        """Calculate gamma factors"""
        self.gamma_long.clear()
        self.gamma_tran.clear()
        L = span * 1000  # [mm]
        middle_layer = len(self.thicknesses) // 2

        for i in range(0, len(self.thicknesses)):
            # print("layer", i+1, ":", self.E_long[i], ",", self.E_tran[i])
            # Thickness of the adjacent layer closer to the CLT centroid [mm]
            if i < middle_layer:
                t_j = self.thicknesses[i + 1]
            else:
                t_j = self.thicknesses[i - 1]

            # Gamma factor
            if i == middle_layer:  # for the middle layer, gamma is always 1
                self.gamma_long.append(1)
                self.gamma_tran.append(1)
            elif self.orientation[i] == 0:  # longitudinal layer
                self.gamma_long.append(
                    1 / (1 + (pi ** 2 * self.E_long[i] * self.thicknesses[i]) / L ** 2 * t_j / self.G_9090))
                self.gamma_tran.append(1)
            else:  # transverse layer
                self.gamma_long.append(1)
                self.gamma_tran.append(
                    1 / (1 + (pi ** 2 * self.E_tran[i] * self.thicknesses[i]) / L ** 2 * t_j / self.G_9090))


class CompositeBeam:
    def __init__(self, beam, slab):
        # General properties
        self.L_1 = SPANS[0]  # span of the WQ-beam [m]
        self.L_2 = SPANS[1]  # span of the CLT slab [m]
        self.w = GAP  # gap between the CLT slab and the WQ-beam [mm]
        self.n_clt = SLAB_AMOUNT  # amount of CLT slabs
        self.rotated = ROTATED  # change to True if slabs are rotated by 90 degrees
        self.s = S  # adjusted parameter [-]

        # Loads
        self.q_h = HUMAN  # human induced load vor vibrations [kg/m^2]
        self.g = GRAVITY  # gravity acceleration [m/s^2]
        self.q = Q * self.n_clt / 2  # linear uniform load to the composite beam [kN/m]

        # WQ-beam
        self.beam = beam

        # CLT slabs
        self.clt = slab

        # Connector properties
        self.k_s = KS  # stiffness of one screw [parallel, perpendicular] to beam direction [N/mm]
        self.l_s = PERIOD  # period of connectors [mm]
        self.n_sc = SCREW_AMOUNT  # total amount of screws per connector (at both sides)

        # Properties of homogenized slab
        self.A_2 = 0  # area [mm2]
        self.I_2 = 0  # second moment of area [mm4]
        self.e_2 = []  # Young's moduli (layer-by-layer) in the beam direction [MPa]
        self.e_22 = []  # Young's moduli (layer-by-layer) perpendicular to the beam direction [MPa]
        self.E_2 = 0  # Young's modulus in the beam direction
        self.E_22 = 0  # Young's modulus perpendicular to the beam direction
        self.EI_eff_2 = 0  # bending stiffness in the beam direction
        self.EI_eff_22 = 0  # bending stiffness perpendicular to the beam direction

        # Distance between "faces" [mm]
        self.a = self.clt.h_l / 2 + self.beam.t_b - self.beam.y_bot

        # Shear factor
        if not self.rotated:
            self.k = self.k_s[0] * self.n_sc * self.a ** 2 / self.l_s
        else:
            self.k = self.k_s[1] * self.n_sc * self.a ** 2 / self.l_s

        # Effective width
        self.B_eff = self.L_1 / self.s * 1000 - 2 * self.w - 2 * self.beam.t_w - self.beam.b_t  # [mm]
        if self.n_clt == 1:
            self.B_eff = self.B_eff / 2

        # Update CLT properties
        self.clt.calculate_properties(self.B_eff, self.L_1)

        # Determine slab properties
        self._slab_properties()

        # Centroids of the composite beam [mm]
        self.y_1 = - self.E_2 * self.A_2 / (self.beam.E_1 * self.beam.A_1 + self.E_2 * self.A_2) * self.a
        self.y_2 = self.beam.E_1 * self.beam.A_1 / (self.beam.E_1 * self.beam.A_1 + self.E_2 * self.A_2) * self.a

        # Total bending stiffness of the composite beam [N*mm^2]
        self.EI = (self.beam.E_1 * self.beam.I_1) + (self.E_2 * self.I_2) + (
                self.y_1 ** 2 * self.beam.E_1 * self.beam.A_1) + self.y_2 ** 2 * self.E_2 * self.A_2

        # Steiner term [N*mm^2]
        self.EI_s = self.EI - self.beam.E_1 * self.beam.I_1 - self.E_2 * self.I_2

        # Coefficients
        self.alpha_1 = self.beam.E_1 * self.beam.I_1 / self.EI_s
        self.alpha_2 = self.E_2 * self.I_2 / self.EI_s
        self.alpha = self.alpha_1 + self.alpha_2
        self.beta = self.EI_s / (self.k * (self.L_1 * 1000) ** 2)
        self.lam = sqrt((1 + self.alpha) / (self.alpha * self.beta))

        # Calculated properties
        self.displacement = None  # vertical displacement of the composite beam at the point x [mm]
        self.sigma_wq = None  # stresses in WQ-beam [MPa]
        self.sigma_clt = None  # stresses in CLT beam [MPa]
        self.force = None  # force in screw [kN]
        self.f = None  # fundamental structural frequency [Hz]

    def _slab_properties(self):
        """Determine the properties of homogenized slab depending on the orientation of the slab"""
        self.A_2 = self.clt.A_2  # area [mm2]
        self.I_2 = self.clt.I_2  # second moment of area [mm4]
        if not self.rotated:  # normal case
            self.E_2 = self.clt.E_2_tran
            self.E_22 = self.clt.E_2_long
            self.EI_eff_2 = self.clt.EI_eff_2_tran
            self.EI_eff_22 = self.clt.EI_eff_2_long
            self.e_2 = self.clt.E_tran
            self.e_22 = self.clt.E_long
        else:  # rotated case
            self.E_2 = self.clt.E_2_long
            self.E_22 = self.clt.E_2_tran
            self.EI_eff_2 = self.clt.EI_eff_2_long
            self.EI_eff_22 = self.clt.EI_eff_2_tran
            self.e_2 = self.clt.E_long
            self.e_22 = self.clt.E_tran

    def vertical_displacement(self, x):
        """
        Return the vertical displacement of the composite beam at the point x
        :param x: coordinate of the measured point along the axis of the beam [m]
        """
        e = x / self.L_1
        self.displacement = self.q * self.L_1 ** 4 / self.EI * (e * (1 - 2 * e ** 2 + e ** 3) / 24 +
                                                                (e * (1 - e)) / (2 * self.alpha * self.lam ** 2)
                                                                - (cosh(self.lam / 2) - cosh(
                    self.lam * (1 - 2 * e) / 2)) /
                                                                (self.alpha * self.lam ** 4 * cosh(
                                                                    self.lam / 2))) * 10 ** 12

    def stresses_wq(self, x):
        """
        Return axial stresses at the top and bottom surfaces of the WQ-beam at the point x
        :param x: coordinate of the measured point along the axis of the beam [m]
        :return: list with stresses at the top and bottom surfaces [MPa]
        """
        # Relative coordinate of the measured point
        e = x / self.L_1

        # Internal moments and forces [kN*mm]
        M_1 = self._moment_i(self.alpha_1, e)
        M_s = self._moment_s(e)
        N_1 = M_s / self.a

        # Stresses [MPa]
        sigma_top = (N_1 / self.beam.A_1 + M_1 / self.beam.W_1top) * 1000
        sigma_bot = (N_1 / self.beam.A_1 + M_1 / self.beam.W_1bot) * 1000

        self.sigma_wq = {"top": sigma_top, "bot": sigma_bot}

    def _moment_i(self, alpha_i, e):
        """Return moment at face i [kN*mm]"""
        M_i = 1000 * self.q * self.L_1 ** 2 * alpha_i / (1 + self.alpha) * (e * (1 - e) / 2
                                                                            + 1 / (self.alpha * self.lam ** 2)
                                                                            * (cosh(self.lam / 2) - cosh(
                    (self.lam * (1 - 2 * e)) / 2)) / cosh(self.lam / 2))

        return M_i

    def _moment_s(self, e):
        M_s = 1000 * self.q * self.L_1 ** 2 / (1 + self.alpha) * (e * (1 - e) / 2
                                                                  - 1 / (self.lam ** 2) * (cosh(self.lam / 2) - cosh(
                    (self.lam * (1 - 2 * e)) / 2)) / cosh(self.lam / 2))
        return M_s

    def stresses_clt(self, x):
        """
        Return the axial stresses along the beam directions in every layer of the CLT slab at the point x
        :param x: coordinate of the measured point along the axis of the beam [m]
        :return sigma_clt: list with axial stresses [kN]
        """
        # Relative coordinate of the measured point
        e = x / self.L_1

        # Internal moments and forces [kN*mm]
        M_2 = self._moment_i(self.alpha_2, e)
        M_s = self._moment_s(e)
        N_2 = - M_s / self.a

        # EA_sum [N]
        EA_sum = 0
        for r in range(0, len(self.clt.thicknesses)):
            EA_sum += self.e_2[r] * self.clt.A[r]

        # Stress [N/mm^2 = MPa]
        self.sigma_clt = {}
        for r in range(0, len(self.clt.thicknesses)):
            sigma_bot = (self.e_2[r] / EA_sum * N_2 + self.e_2[r] / self.EI_eff_2 * M_2 * (
                -(self.clt.zeta[r] - self.clt.thicknesses[r] / 2))) * 1000
            sigma_top = (self.e_2[r] / EA_sum * N_2 + self.e_2[r] / self.EI_eff_2 * M_2 * (
                -(self.clt.zeta[r] + self.clt.thicknesses[r] / 2))) * 1000
            self.sigma_clt[r + 1] = {"top": sigma_top, "bot": sigma_bot}

    def screw_force(self, x):
        """
        Return the connector force per single screw at the point x
        :param x: coordinate of the measured point along the axis of the beam [m]
        """
        # Relative coordinate of the measured point
        e = x / self.L_1

        # Q_s [kN]
        Q_s = self.q * self.L_1 / (1 + self.alpha) * ((1 - 2 * e) / 2 - (sinh(self.lam * (1 - 2 * e) / 2))
                                                      / (self.lam * cosh(self.lam / 2)))

        # Connector force [kN]
        F_d = Q_s * self.l_s / self.a

        # Connector force per single screw [kN]
        self.force = F_d / self.n_sc

    def vibration(self):
        # Mass of the WQ-beam per m length [kg/m]
        m_1 = self.beam.A_1 * self.beam.rho / 10 ** 6

        # Mass of the CLT slab per m length [kg/m]
        m_2 = self.clt.h_l / 1000 * self.B_eff / 1000 * self.clt.rho

        # Human induced load [kg/m]
        m_h = self.B_eff / 1000 * self.q_h

        # Mass per unit length of the composite beam [kg/m]
        mu = m_1 + m_2 + m_h

        # Mass m [kg/m^2]
        m = self.clt.rho * self.clt.h_l / 1000 + self.q_h

        # Fundamental structural frequency of the wq beam [Hz]
        i = 1
        f_0 = 1 / (2 * pi) * sqrt(self.EI_s / (mu * self.L_1 ** 4)
                                  * (1 + self.alpha + self.alpha * self.beta * i ** 2 * pi ** 2)
                                  / (1 + self.beta * i ** 2 * pi ** 2) * i ** 4 * pi ** 4) / 1000

        # Second moment of inertia per m [mm^4/m]
        I_22 = self.I_2 / self.B_eff * 1000

        # Fundamental structural frequency of one CLT slab [Hz]
        f_1 = pi / (2 * self.L_2 ** 2) * sqrt(self.E_22 * I_22 / m) / 1000

        # Frequency of the total structural system [Hz]
        self.f = sqrt(1 / (1 / f_0 ** 2 + 1 / f_1 ** 2 + 1 / f_1 ** 2))

    def print_input_info(self):
        print("Initial data:")
        if self.rotated:
            print("Slabs are rotated by 90 degrees")
        if self.n_clt == 1:
            print("Single CLS slab, asymmetrical loading")
        if not self.clt.gamma_used:
            print("CLT properties are calculated ignoring Gamma-method")

        print("{:<12} {:>5} ".format("s", format(self.s, "0.1f")))
        print("{:<12} {:>5} {:<}".format("a", format(self.a, "0.1f"), "mm"))
        print("{:<12} {:>5} {:<}".format("B_eff", format(self.B_eff, "0.0f"), "mm"))

    def print_clt_properties(self):
        print()
        print("CLT properties:")
        print("{:<12} {:>5} {:<}".format("E_2", format(self.E_2, "0.0f"), "MPa"))
        print("{:<12} {:>5} {:<}".format("E_22", format(self.E_22, "0.0f"), "MPa"))

    def print_results(self):
        print()
        print("Results:")
        if self.displacement:
            print("{:<12} {:>5} {:<}".format("δ_wq", format(self.displacement, "0.1f"), "mm"))

        if self.sigma_wq:
            print("{:<12} {:>5} {:<}".format("σ_wq,top", format(self.sigma_wq["top"], ".0f"), "MPa"))
            print("{:<12} {:>5} {:<}".format("σ_wq,bot", format(self.sigma_wq["bot"], ".0f"), "MPa"))

        if self.sigma_clt:
            print("{:<12} {:>5} {:<}".format("σ_clt,top",
                                             format(self.sigma_clt[len(self.clt.thicknesses)]["top"], ".2f"), "MPa"))
            print("{:<12} {:>5} {:<}".format("σ_clt,bot", format(self.sigma_clt[1]["bot"], ".2f"), "MPa"))

        if self.force:
            print("{:<12} {:>5} {:<}".format("F_scr", format(self.force, "0.2f"), "kN"))

        if self.f:
            print("{:<12} {:>5} {:<}".format("f", format(self.f, "0.2f"), "Hz"))

    def print_clt_stresses(self):
        print()
        if self.sigma_clt:
            print("Stresses in CLT [MPa]:")
            print("{:<6} {:>6} {:>7}".format("Layer", "Top", "Bottom"))
            for layer in sorted(self.sigma_clt, reverse=True):
                print("{:<6} {:>6} {:>7}".format(layer, format(self.sigma_clt[layer]["top"], "0.2f"),
                                                 format(self.sigma_clt[layer]["bot"], "0.2f")))
        else:
            print("CLT stresses were not calculated and cannot be printed")
