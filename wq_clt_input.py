# WQ-beam dimensions [mm]
TOP_FLANGE = [280, 14]
BOTTOM_FLANGE = [490, 10]
WALL = [272, 5]

# CLT properties
LAYERS = [80, 40, 40, 40, 80]       # thicknesses of layers from bottom to top [mm]
ORIENTATION = [0, 90, 0, 90, 0]     # orientation of layers from bottom to top (0 if perpendicular to beam)
YOUNGS_MODULI = [11600, 730]        # Young's moduli in [longitudinal, transverse] directions [MPa]
ROLLING_MODULUS = 50                # rolling shear modulus [MPa]
GAMMA = True                        # change to True if CLT bending stiffness is calculated with gamma-factors

# Composite slab properties
SPANS = [7, 7]      # spans of [WQ-beam, CLT slab] [m]
GAP = 20            # gap between the CLT slab and the WQ-beam [mm]
SLAB_AMOUNT = 2     # amount of CLT slabs
ROTATED = False     # change to True if slabs are rotated by 90 degrees
S = 2.0             # adjusted parameter [-]

# Connector properties
KS = [1330, 2960]   # stiffness of one screw [parallel, perpendicular] to beam direction [N/mm]
PERIOD = 600        # period of connectors [mm]
SCREW_AMOUNT = 2    # total amount of screws per connector (at both sides)

# Loads
Q = 89.0            # linear load to the composite beam [kN/m]
HUMAN = 25          # human induced load for vibrations [kg/m^2]
