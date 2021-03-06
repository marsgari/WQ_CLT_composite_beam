# WQ-beam dimensions [mm]
TOP_FLANGE = [280, 14]
BOTTOM_FLANGE = [490, 10]
WALL = [272, 5]
STEEL_YOUNG = 210  # Young's modulus [GPa]
STEEL_DENSITY = 7850  # density [kg/m^3]

# CLT properties
LAYERS = [40, 40, 40, 40, 40, 40, 40]   # thicknesses of layers from bottom to top [mm]
ORIENTATION = [0, 0, 90, 0, 90, 0, 0]   # orientation of layers from bottom to top (0 if perpendicular to beam)
CLT_YOUNG = [11600, 730]                # Young's moduli in [longitudinal, transverse] directions [MPa]
ROLLING_MODULUS = 50                    # rolling shear modulus [MPa]
CLT_DENSITY = 500                       # density [kg/m3]
GAMMA = True                            # change to True if CLT bending stiffness is calculated with gamma-factors

# Composite slab properties
SPANS = [7, 7]      # spans of [WQ-beam, CLT slab] [m]
GAP = 20            # gap between the CLT slab and the WQ-beam [mm]
SLAB_AMOUNT = 2     # amount of CLT slabs
ROTATED = False     # change to True if slabs are rotated by 90 degrees
S = 1.65             # adjusted parameter [-]

# Connector properties
KS = [1330, 2960]   # stiffness of one screw [parallel, perpendicular] to beam direction [N/mm]
PERIOD = 600        # period of connectors [mm]
SCREW_AMOUNT = 2    # total amount of screws per connector (at both sides)

# Loads
Q = 89.0            # linear load to the composite beam [kN/m]
HUMAN = 25          # human induced load for vibrations [kg/m^2]
GRAVITY = 10        # gravity acceleration [m/s^2]
