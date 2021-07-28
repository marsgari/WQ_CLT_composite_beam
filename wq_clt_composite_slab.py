from math import sqrt, cosh, sinh, pi

# WQ properties
TOP_FLANGE = [280, 14]
BOTTOM_FLANGE = [490, 10]
WALL = [272, 5]

# CLT properties
LAYERS = [80, 40, 40, 40, 80]  # thicknesses of layers from bottom to top [mm]
ORIENTATION = [0, 90, 0, 90, 0]  # orientation of layers from bottom to top (0 if perpendicular to beam)
YOUNG = [11600, 730]  # Young's moduli in [longitudinal, transverse] directions [MPa]
ROLLING = 50  # rolling shear modulus [MPa]
GAMMA = False  # change to True if CLT bending stiffness is calculated with gamma-factors

# Composite slab properties
SPANS = [7, 7]  # spans [L1, L2] [m]
GAP = 20  # gap between the CLT slab and the WQ-beam [mm]
SLAB_AMOUNT = 2  # amount of CLT slabs
ROTATED = False  # change to True if slabs are rotated by 90 degrees
S = 2.0  # adjusted parameter [-]

# Connector properties
KS = [1330, 2960]  # stiffness of one screw [parallel, perpendicular] to beam direction [N/mm]
PERIOD = 600  # period of connectors [mm]
SCREW_AMOUNT = 2  # total amount of screws per connector (at both sides)

# Loads
Q = [3.0, 5.0]  # characteristic uniform [dead, live] loads [kN/m^2]
K = [1.15, 1.50]  # load factors for [dead, live] loads
HUMAN = 25  # human induced load vor vibrations [kg/m^2]


class WQBeam:
    def __init__(self):
        self.b_t = TOP_FLANGE[0]  # top flange width [mm]
        self.t_t = TOP_FLANGE[1]  # top flange thickness [mm]
        self.b_b = BOTTOM_FLANGE[0]  # bottom flange width [mm]
        self.t_b = BOTTOM_FLANGE[1]  # bottom flange thickness [mm]
        self.h_b = WALL[0]  # wall height [mm]
        self.t_w = WALL[1]  # wall thickness [mm]
        self.E_1 = 210000  # steel Young's modulus [MPa]
        self.rho = 7850  # steel density [kg/m^3]

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
        self.rho = 500  # wood density [kg/m3]
        self.E_0 = YOUNG[0]  # Young's modulus along grains [MPa]
        self.E_90 = YOUNG[1]  # Young's modulus transverse to grains [MPa]
        self.G_9090 = ROLLING  # rolling shear modulus [MPa]
        self.use_gamma = GAMMA  # change to True if CLT bending stiffness is calculated with gamma-factors

        # Properties of the slab layer-by-layer
        self.E_lon = []  # Young's moduli in the beam direction
        self.E_tran = []  # Young's moduli perpendicular to the beam direction
        self.A = []  # areas [mm2]
        self.J = []  # second moments of area in relation to own centroid [mm4]
        self.zeta = []  # distances from own centroid to the centroid of CLT

        # Properties of homogenized slab
        self.A_2 = 0  # area [mm2]
        self.I_2 = 0  # second moment of area [mm4]
        self.E_2 = 0  # Young's modulus in the beam direction
        self.E_22 = 0  # Young's modulus perpendicular to the beam direction
        self.EI_eff_2 = 0  # bending stiffness in the beam direction
        self.EI_eff_22 = 0  # bending stiffness perpendicular to the beam direction

    def calculate_properties(self, B_eff, L_1, rotated):
        # Get lists with A, J and zeta
        self._make_lists(B_eff)
        self.E_long, self.E_tran = self._clt_young_moduli(rotated)

        # Gamma factors
        gamma_2, gamma_22 = self._clt_gamma(L_1)

        # EI_eff [N*mm^2]
        for r in range(0, len(self.thicknesses)):
            if not self.use_gamma:
                self.EI_eff_2 += self.E_long[r] * self.J[r] + self.zeta[r] ** 2 * self.E_long[r] * self.A[r]
                self.EI_eff_22 += self.E_tran[r] * self.J[r] + self.zeta[r] ** 2 * self.E_tran[r] * self.A[r]
            else:
                self.EI_eff_2 += self.E_long[r] * self.J[r] + gamma_2[r] * self.zeta[r] ** 2 * self.E_long[r] * self.A[r]
                self.EI_eff_22 += self.E_tran[r] * self.J[r] + gamma_22[r] * self.zeta[r] ** 2 * self.E_tran[r] * self.A[r]

        # Update CLT properties
        self.A_2 = B_eff * self.h_l  # mm^2
        self.I_2 = B_eff * self.h_l ** 3 / 12  # mm^4
        self.E_2 = self.EI_eff_2 / self.I_2  # N/mm^2 = MPa
        self.E_22 = self.EI_eff_22 / self.I_2  # N/mm^2 = MPa

    def _make_lists(self, B_eff):
        lower = -self.h_l / 2  # distance to bottom of CLT layer [mm]
        for t in self.thicknesses:
            self.A.append(B_eff * t)
            self.J.append(B_eff * t ** 3 / 12)
            self.zeta.append(lower + t / 2)
            lower += t

    def _clt_young_moduli(self, rotated):
        E_lon = []
        E_tran = []
        for angle in self.orientation:
            if angle == 0:
                E_lon.append(self.E_90)
                E_tran.append(self.E_0)
            else:
                E_lon.append(self.E_0)
                E_tran.append(self.E_90)

        if not rotated:
            return [E_lon, E_tran]
        else:
            return [E_tran, E_lon]

    def _clt_gamma(self, L_1):
        gamma_long = []
        gamma_tran = []
        L = L_1 * 1000  # [mm]
        middle_layer = len(self.thicknesses) // 2
        for i in range(0, len(self.thicknesses)):
            # print("layer", i+1, ":", E_lon[i], ",", E_tran[i])
            if i < middle_layer:
                t_j = self.thicknesses[i + 1]
            else:
                t_j = self.thicknesses[i - 1]
            if i == middle_layer:
                gamma_long.append(1)
                gamma_tran.append(1)
            elif self.orientation[i] == 90:  # longitudinal layer
                gamma_long.append(
                    1 / (1 + (pi ** 2 * self.E_long[i] * self.thicknesses[i]) / L ** 2 * t_j / self.G_9090))
                gamma_tran.append(1)
            else:  # transverse layer
                gamma_long.append(1)
                gamma_tran.append(
                    1 / (1 + (pi ** 2 * self.E_tran[i] * self.thicknesses[i]) / L ** 2 * t_j / self.G_9090))

        return gamma_long, gamma_tran


class CompositeBeam:
    def __init__(self, beam, slab):
        # General properties
        self.L_1 = SPANS[0]  # span of the WQ-beam [m]
        self.L_2 = SPANS[0]  # span of the CLT slab [m]
        self.w = GAP  # gap between the CLT slab and the WQ-beam [mm]
        self.n_clt = SLAB_AMOUNT  # amount of CLT slabs
        self.rotated = ROTATED  # change to True if slabs are rotated by 90 degrees
        self.s = S  # adjusted parameter [-]

        # Loads
        self.q_dead = Q[0]  # characteristic uniform dead load [kN/m^2]
        self.q_live = Q[1]  # characteristic uniform live load [kN/m^2]
        self.k_dead = K[0]  # load factor for dead loads
        self.k_live = K[1]  # load factor for live loads
        self.q_h = HUMAN  # human induced load vor vibrations [kg/m^2]
        self.g = 10  # gravity acceleration [m/s^2]

        # WQ-beam
        self.beam = beam

        # CLT slab
        self.clt = slab

        # Connector properties
        self.k_s = KS  # stiffness of one screw [parallel, perpendicular] to beam direction [N/mm]
        self.l_s = PERIOD  # period of connectors [mm]
        self.n_sc = SCREW_AMOUNT  # total amount of screws per connector (at both sides)

        # Distance between "faces" [mm]
        self.a = self.clt.h_l / 2 + self.beam.t_b - self.beam.y_bot

        # Shear factor
        self.k = self._shear_factor()

        # Effective width
        self.B_eff = self._effective_width()

        # Update CLT properties
        self.clt.calculate_properties(self.B_eff, self.L_1, self.rotated)

        # Linear uniform load applied to the WQ-beam beam [kN/m]
        self.q = self._calculate_load()

        # Centroids of the composite beam [mm]
        self.y_1 = self.clt.E_2 * self.clt.A_2 / (self.beam.E_1 * self.beam.A_1 + self.clt.E_2 * self.clt.A_2) * self.a
        self.y_2 = - self.beam.E_1 * self.beam.A_1 / (self.beam.E_1 * self.beam.A_1 + self.clt.E_2 * self.clt.A_2) * self.a

        # Total bending stiffness of the composite beam [N*mm^2]
        self.EI = (self.beam.E_1 * self.beam.I_1) + (self.clt.E_2 * self.clt.I_2) + (
                    self.y_1 ** 2 * self.beam.E_1 * self.beam.A_1) + self.y_2 ** 2 * self.clt.E_2 * self.clt.A_2

        # Steiner term [N*mm^2]
        self.EI_s = self.EI - self.beam.E_1 * self.beam.I_1 - self.clt.E_2 * self.clt.I_2

        # Coefficients
        self.alpha_1 = self.beam.E_1 * self.beam.I_1 / self.EI_s
        self.alpha_2 = self.clt.E_2 * self.clt.I_2 / self.EI_s
        self.alpha = self.alpha_1 + self.alpha_2
        self.beta = self.EI_s / (self.k * (self.L_1 * 1000) ** 2)
        self.lam = sqrt((1 + self.alpha) / (self.alpha * self.beta))

    def _shear_factor(self):
        if not self.rotated:
            k = self.k_s[0] * self.n_sc * self.a ** 2 / self.l_s
        else:
            k = self.k_s[1] * self.n_sc * self.a ** 2 / self.l_s
        return k

    def _effective_width(self):
        B_eff = self.L_1 / self.s * 1000 - 2 * self.w - 2 * self.beam.t_w - self.beam.b_t  # [mm]
        if self.n_clt == 1:
            B_eff = B_eff / 2
        return B_eff

    def _calculate_load(self):
        """
        Return the total linear load applied to the WQ-beam.
        The load includes the weight of the WQ-beam and loads from the half of the span of the CLT slab.
        The command also calculates but not returns the loads [kN/m^2] applied separately to the WQ-beam
        and the CLT slab (for FEM)
        :return q: the total linear load applied to the WQ-beam [kN/m]
        """
        # Total dead load to CLT slabs [kN/m^2]
        q_dead_total = self.q_dead + self.clt.h_l / 1000 * self.clt.rho * self.g / 1000

        # Uniformly distributed load to the CLT slab [kN/m^2]
        q_clt = q_dead_total * self.k_dead + self.q_live * self.k_live

        # Mass of the WQ-beam per m length [kg/m]
        m_wq = self.beam.A_1 * self.beam.rho / 10 ** 6

        # Width of the top flange of the WQ-beam
        b_wq = self.beam.b_t + 2 * self.beam.t_w

        # Uniformly distributed load to the WQ-beam [kN/m^2]
        q_wq = q_clt * (b_wq + 2 * self.w) / b_wq + self.k_dead * m_wq * self.g / b_wq

        # Width of the applied load
        b_clt = self.L_2
        if self.n_clt == 1:
            b_clt = self.L_2 / 2

        # Linear uniform load applied to the WQ-beam beam [kN/m]
        q = q_clt * b_clt + self.k_dead * m_wq * self.g / 1000

        return q

    def vertical_displacement(self, x):
        """
        Return the vertical displacement of the composite beam at the point x
        :param x: coordinate of the measured point along the axis of the beam [m]
        :return v: vertical displacement of the composite beam [mm]
        """
        # Relative coordinate of the measured point
        e = x / self.L_1

        # Displacement of the composite beam [mm]
        v = self.q * self.L_1 ** 4 / self.EI * (e * (1 - 2 * e ** 2 + e ** 3) / 24 +
                                 (e * (1 - e)) / (2 * self.alpha * self.lam ** 2)
                                 - (cosh(self.lam / 2) - cosh(self.lam * (1 - 2 * e) / 2)) /
                                 (self.alpha * self.lam ** 4 * cosh(self.lam / 2))) * 10 ** 12

        return v

    def stresses_wq(self, x):
        """
        Return axial stresses at the top and bottom surfaces of the WQ-beam at the point x
        :param x: coordinate of the measured point along the axis of the beam [m]
        :return: list with stresses at the top and bottom surfaces [MPa]
        """
        # Relative coordinate of the measured point
        e = x / self.L_1

        # Moments [kN*mm]
        M_1 = self._moment_i(self.alpha_1, e)
        M_s = self._moment_s(e)

        # Axial force [kN]
        N_1 = M_s / self.a

        # Stresses [N/mm^2 = MPa]
        sigma_top = (N_1 / self.beam.A_1 + M_1 / self.beam.W_1top) * 1000
        sigma_bot = (N_1 / self.beam.A_1 + M_1 / self.beam.W_1bot) * 1000

        return {"top": sigma_top, "bot": sigma_bot}

    def _moment_i(self, alpha_i, e):
        """Return moment at face i [kN*mm]"""
        M_i = 1000 * self.q * self.L_1 ** 2 * alpha_i / (1 + self.alpha) * (e * (1 - e) / 2
                                                                            + 1 / (self.alpha * self.lam ** 2)
             * (cosh(self.lam / 2) - cosh((self.lam * (1 - 2 * e)) / 2)) / cosh(self.lam / 2))

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

        # Moments [kN*mm]
        M_2 = self._moment_i(self.alpha_2, e)
        M_s = self._moment_s(e)

        # Axial force [kN]
        N_2 = - M_s / self.a

        # EA_sum [N]
        EA_sum = 0
        for r in range(0, len(self.clt.thicknesses)):
            EA_sum += self.clt.E_long[r] * self.clt.A[r]

        # Stress [N/mm^2 = MPa]
        sigma_clt = {}
        for r in range(0, len(self.clt.thicknesses)):
            sigma_bot = (self.clt.E_long[r] / EA_sum * N_2 +
                              self.clt.E_long[r] / self.clt.EI_eff_2 * M_2 * (
                                  -(self.clt.zeta[r] - self.clt.thicknesses[r] / 2))) * 1000
            sigma_top = (self.clt.E_long[r] / EA_sum * N_2 +
                              self.clt.E_long[r] / self.clt.EI_eff_2 * M_2 * (
                                  -(self.clt.zeta[r] + self.clt.thicknesses[r] / 2))) * 1000
            sigma_clt[r+1] = {"top": sigma_top, "bot": sigma_bot}
        return sigma_clt

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
        F_0 = F_d / self.n_sc

        return F_0

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

        # Interger i
        i = 1

        # Fundamental structural frequency of the wq beam [Hz]
        f_0 = 1 / (2 * pi) * sqrt(self.EI_s / (mu * self.L_1 ** 4)
                                  * (1 + self.alpha + self.alpha * self.beta * i ** 2 * pi ** 2)
                                  / (1 + self.beta * i ** 2 * pi ** 2) * i ** 4 * pi ** 4) / 1000

        # Second moment of inertia per m [mm^4/m]
        I_22 = self.clt.I_2 / self.B_eff * 1000

        # Fundamental structural frequency of one CLT slab [Hz]
        f_1 = pi / (2 * self.L_2 ** 2) * sqrt(self.clt.E_22 * I_22 / m) / 1000

        # Frequency of the total structural system [Hz]
        f = sqrt(1 / (1 / f_0 ** 2 + 1 / f_1 ** 2 + 1 / f_1 ** 2))

        return f

    def print_stresses_clt(self, sigma_clt):
        """
        Print axial stresses in the CLT slab from the given list.
        Stresses are printed from the top layer to the bottom layer.
        :param sigma_clt: list with axial stresses [MPa]
        """
        n = len(self.clt.thicknesses)
        print("Stresses in CLT [MPa]:")
        for layer in sorted(sigma_clt, reverse=True):
            print("Layer ", layer, ": ", "σ_top = ", format(sigma_clt[layer]["top"], "0.2f"), sep="")
            print("         σ_bot =", format(sigma_clt[layer]["bot"], "0.2f"))


def main():

    # Make the composite beam
    beam = WQBeam()
    clt = CLT()
    composite_beam = CompositeBeam(beam, clt)

    # Coordinates of desired points
    x1 = composite_beam.L_1 / 2  # coordinate of middle of the span
    x2 = 0  # coordinate of the support

    # Calculate properties
    delta_wq = composite_beam.vertical_displacement(x1)  # calculate vertical displacements
    sigma_wq = composite_beam.stresses_wq(x1)  # calculate stresses in WQ-beam
    sigma_clt = composite_beam.stresses_clt(x1)  # calculate stresses in CLT slab
    F_scr = composite_beam.screw_force(x2)  # calculate force in screw
    f = composite_beam.vibration()  # calculate frequency

    # Print results
    if composite_beam.rotated:
        print("       Rotated!")
    if composite_beam.n_clt == 1:
        print("     Single slab!")
    print("{:>11} {:<} ".format("s =", format(composite_beam.s, "0.2f")))
    if 'delta_wq' in locals():
        print("{:>11} {:<} {:<}".format("δ_wq =", format(delta_wq, "0.1f"), "mm"))
    if 'sigma_wq' in locals():
        print("{:>11} {:<} {:<}".format("σ_wq,top =", format(sigma_wq["top"], ".0f"), "MPa"))
        print("{:>11} {:<} {:<}".format("σ_wq,bot =", format(sigma_wq["bot"], ".0f"), "MPa"))
    if 'sigma_clt' in locals():
        print("{:>11} {:<} {:<}".format("σ_clt,top =",
                                        format(sigma_clt[len(composite_beam.clt.thicknesses)]["top"], ".2f"), "MPa"))
        print("{:>11} {:<} {:<}".format("σ_clt,bot =", format(sigma_clt[1]["bot"], ".2f"), "MPa"))
    if 'F_scr' in locals():
        print("{:>11} {:<} {:<}".format("F_scr =", format(F_scr, "0.2f"), "kN"))
    if 'f' in locals():
        print("{:>11} {:<} {:<}".format("f =", format(f, "0.2f"), "Hz"))
    #composite_beam.print_stresses_clt(sigma_clt)


main()
