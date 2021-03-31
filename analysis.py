from make_parameters import Paramgroup as PG
import numpy as np
def externalinfluence():
    pg = PG()
    pg.default("general:init_f", [.5])
    pg.default("phage:generalist", [False])
    pg.default("phage:genotype", [1])
    pg.default("phage:burst", [120])
    pg.default("general:output_frequency", [100000])

    pg.default("phage_s:generalist", [False])
    pg.default("phage_s:genotype", [1])
    pg.default("phage_s:burst", [120])
    pg.add("freq", "phage:incubation_period", [.01,.02,.03,.05,.07,.1,.25,.5,1])
    #pg.add("freq", "phage:incubation_period", [0.05,.06,.07,.07,.08,.09,.1,.12,.15,.17,.2,.22,.25])
    #pg.add("freq", "phage:incubation_period", [0.02])
    pg.add("freq","general:init_f", [.05,.1,.3,.5,.7,.9,.95])

   
   


    pg.set_repeats(20)
    return pg

if __name__ == '__main__':
    externalinfluence().print()