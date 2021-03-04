from make_parameters import Paramgroup as PG
import numpy as np
def externalinfluence():
    pg = PG()
    pg.default("phage:genotype", [1])
    pg.default("species1:genotype", [1])
    pg.default("species2:genotype", [1])
    pg.default("general:output_frequency", [1])
    pg.default("species1:mut_1", [0])
    pg.default("species1:mut_0", [0])

    pg.default("species2:mut_1", [0])
    pg.default("species2:mut_0", [0])

    pg.default("species1:mu_0", [24.5])
    pg.default("species1:mu_1", [24.5])

    pg.default("species2:mu_0", [24.5])
    pg.default("species2:mu_1", [24.5])


    pg.add('mam','species2:mut_1', [0,10e-6,10e-3,.1,.5,1])  
    pg.add('mam','species1:mut_1', [0,10e-6,10e-3,.1,.5,1])  
    #pg.add('mam','species2:mut_0', [0,1,10e-6,10e-3,10e-1])  

  
    #pg.add('cost','species1:mu_0', [14, 23, 24.5, 26]) 
    #pg.add('cost','species2:mu_0', [14, 23, 24.5, 26])

    #pg.add('cost','species2:mut_1', [0,1,10e-6])  
    #pg.add('cost','species2:mut_0', [0,1,10e-6])  


    pg.set_repeats(10)
    return pg

if __name__ == '__main__':
    
    externalinfluence().print()