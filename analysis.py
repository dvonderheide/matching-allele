from make_parameters import Paramgroup as PG
import numpy as np
def externalinfluence():
    pg = PG()
    pg.default("general:runtime", 20)
    
    pg.default('species2:genotype', [(1,1),(1,0),(0,1),(0,0)])
    pg.default('species1:genotype', [(1,1),(1,0),(0,1),(0,0)])
    pg.default('species2:mut_rate', [0, 10e-6, 10e-2,1])  
    pg.default('species1:mut_rate', [0, 10e-6, 10e-2,1])


    pg.set_repeats(10)
    return pg

if __name__ == '__main__':
    print(len(externalinfluence()))
    externalinfluence().print()