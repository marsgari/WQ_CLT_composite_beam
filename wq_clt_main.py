from wq_clt_classes import *


def main():

    # Make the composite beam
    beam = WQBeam()
    clt = CLT()
    composite_beam = CompositeBeam(beam, clt)

    # Coordinates of desired points [m]
    x1 = 3.5
    x2 = 0

    # Design of the composite beam
    composite_beam.vertical_displacement(x1)  # calculate vertical displacements
    composite_beam.stresses_wq(x1)  # calculate stresses in WQ-beam
    composite_beam.stresses_clt(x1)  # calculate stresses in CLT slab
    composite_beam.screw_force(x2)  # calculate force in screw
    composite_beam.vibration()  # calculate fundamental structural frequency

    # Print results
    composite_beam.print_input_info()
    composite_beam.print_clt_properties()
    composite_beam.print_results()
    composite_beam.print_clt_stresses()

    # Make prd-report
    composite_beam.make_report()


main()
