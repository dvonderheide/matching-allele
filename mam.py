import argparse
import os
import simbiofilm as sb
import numpy as np
from random import random

# There may be some warnings from libraries we use. Feel free to enable this
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)


class InfectAndEndController:
    def __init__(self, height, count, duration=10, biomass=False):
        """
        If biomass is given, interpret height as a biomass
        """
        self.infection_height = height
        self.infection_count = count
        self.end_time = 10000
        self.infection_length = duration
        self._use_mass = biomass

    def infect_point(self, space, time, phage, biomass, infectable):
        self.end_time = time + self.infection_length
        if self.infection_count > 1:
            sb.infect_point(space, time, self.infection_count, phage, infectable)
        else:
            sb.infect_at(space, time, phage, biomass, space.shape[1])

    def condition(self, space, time, phage, biomass, infectable):
        mass = sb.to_grid(biomass, "mass")
        if self._use_mass:
            return mass.sum() >= self.infection_height
        return (np.max(np.where(mass > 0)[0]) * space.dl) >= self.infection_height

    def end_condition(self, space, time):
        return time >= self.end_time
class Bacteria(sb.Bacteria):
    def phage_interacted(self, bacterium, phage):
        # see simbiofilm.behaviors.phage_interaction
        # Return true if we should remove the particle
        if not bacterium == self.with_id(bacterium.id):
            msg = "Individual does not appear to be in this container."
            raise RuntimeError(msg)
        return [bacterium.genotype == phage.genotype, bacterium.adsorbable]

    def clone(self, parent, **params):
        """Clones parent exactly."""
        parameters = dict(zip(parent.dtype.names, parent))
        if params:
            parameters.update(params)
        
        #parameters['genotype'] = np.random.choice(4, p=parent.mut_rate[parent.genotype])
        mut = random()

        if parent.genotype == 1 and mut < parent.mut_1:
            parameters['genotype'] = 0
            parameters['mu'] =  parent.mu_0
        elif mut < parent.mut_0:
            parameters['genotype'] = 1
            parameters['mu'] =  parent.mu_1

        return self.add_individual(parent.location, parameters)


def config():
    return sb.cfg_from_dict(
        {
            "general": {
                "output_frequency": 1,
                "nconnections": 1,
                "seed": 1235,
                "runtime": 20,
                "init_count": 150,
                "init_f": 0.95,
                "line_shove": False,
                "max_sim_time": 100,
                "3D": False,
                "connectivity": 0.05,
                "fixed_dt": 1 / 24 + 0.0000001,
                "targetone": True,
            },
            "space": {"dl": 3e-6, "shape": (75, 200), "well_mixed": False},
            "infection": {'count': 120, 'height': 20e-6, 'duration': 1.5},
            "substrate": {"max": 6, "diffusivity": 2e-5, "K": 1.18, "h": 15e-6},
            "erosion": {"biomass_rate": 7.5e8, "phage_rate": 2.4e12},
            "species1": {
                "density": 200e3,
                "mass": 1e-12,
                "division_mass": 1.333e-12,
                "impedance": 6,
                "mu": 24.5,
                "adhesion": 1,
                "resistant": False,
                "adsorbable": False,
                "yield_s": 0.495,
                "genotype": 1,
                "mut_1": 1,
                "mut_0": 0,
                "mu_0": 24.5,
                "mu_1": 24.5,
            },
            "species2": {
                "density": 200e3,
                "mass": 1e-12,
                "division_mass": 1.333e-12,
                "impedance": 6,
                "mu": 24.5,
                "adhesion": 1,
                "resistant": False,
                "adsorbable": False,
                "yield_s": 0.495,
                "genotype": 1,
                "mut_1": 1,
                "mut_0": 0,
                "mu_0": 24.5,
                "mu_1": 24.5,
            },
            "infected": {
                "density": 200e3,
                "mass": 1e-12,
                "division_mass": 1.333e-12,
                "impedance": 6,
                "adsorbable": False, # not sure if used
                "adhesion": 1,
            },
            "phage": {
                "diffusivity": 3.30e-6,
                "adsorption_rate": 70,
                "burst": 120,
                "incubation_period": 0.02,
                "adhesion": 1,
                "genotype": 1,
                "a": 0,
                "k": 10e-2,
            },
        }
    )


def setup(cfg, outdir="tmp"):
    """Do the thing."""

    space = sb.Space(cfg.space)
    sim = sb.Simulation(space, cfg)

    substrate = sb.Solute("solute", space, cfg.substrate)

    activeSpecies = [
        Bacteria("species1", space, cfg.species1, extragroups=['Susceptible']),
        Bacteria("species2", space, cfg.species2, extragroups=['Producer']),
    ]
    phage = sb.Phage("phage", space, cfg.phage)
    infected = sb.InfectedBacteria("infected", space, cfg.infected, cfg.phage)
    pairs = {phage: infected}

    sim.add_container(substrate, *activeSpecies, infected, phage)

    sb.inoculate_at(space, 0, activeSpecies[0], int(cfg.general.init_count * cfg.general.init_f))
    sb.inoculate_at(space, 0, activeSpecies[1], int(cfg.general.init_count * (cfg.general.init_f)))

    sb.initialize_bulk_substrate(space, 0, substrate)

    # Set up infection & end of simulation
    ic = InfectAndEndController(
        cfg.infection.height,
        int(cfg.infection.count),
        cfg.infection.duration,
        cfg.space.well_mixed,
    )
    sim.add_event(ic.infect_point, ic.condition, [phage, activeSpecies, activeSpecies[0]])
    # change infection duration in InfectAndEndController
    sim.add_end_condition(ic.end_condition, f"FIN-IC: {10} days after infection start")

    reactions = [sb.MonodRateReaction(substrate, sp) for sp in activeSpecies]

    sim.add_behavior(
        sb.DiffusionReaction(substrate, reactions),
        sb.Biomass_growth(reactions),
        sb.Bacterial_division(),
        sb.Lysis(pairs),
        sb.Erode_individuals(cfg.erosion.biomass_rate),
        sb.Phage_randomwalk(cfg.erosion.phage_rate),
        sb.Detach_biomass(),
        sb.Relax_biofilm(),
        sb.Phage_interaction(pairs),
    )

    sim.initialize(f"{outdir}/run1", [], cfg.general.output_frequency)
    return sim


def main():
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", metavar=("names", "values"), nargs=2)
    parser.add_argument("-o", metavar="output_dir")

    args = parser.parse_args()

    cfg = config()
    if args.p:
        cfg = sb.parse_params(args.p[0], args.p[1], cfg)
    name = args.o if args.o else f"local_runs/{sys.argv[0][:-3]}"
    sim = setup(cfg, name)
    try:
        sim.iterate(cfg.general.runtime, dt_max=0.02)
    except KeyboardInterrupt as kerr:
        sim.finish()
        raise kerr
    sim.finish()


if __name__ == "__main__":
    main()
