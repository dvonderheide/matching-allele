from make_parameters import Paramgroup as PG
import numpy as np
def externalinfluence():
    pg = PG()
    pg.default("general:init_f", [.5])
    pg.default("general:output_frequency", [1])
    pg.default("phage:generalist", [True])
    pg.default("phage:genotype", [1])
    pg.default("phage:burst", [60])

    #pg.default("phage_s:generalist", [False])
    #pg.default("phage_s:genotype", [1])
    #pg.default("phage_s:burst", [120])
    
    pg.add("freq","general:init_f", [.05,.1,.5,.9,.95])
    pg.add("freq","general:init_f2", [.05,.1,.5,.9,.95])
   
   


    pg.set_repeats(30)
    return pg

if __name__ == '__main__':
    externalinfluence().print()