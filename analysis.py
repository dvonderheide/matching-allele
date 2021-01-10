from make_parameters import Paramgroup as PG
import numpy as np
def externalinfluence():
    pg = PG()
    pg.default("species2:genotype", [(1,1), (0,0)])
    pg.default("species1:genotype", [(1,1), (0,0)])
    pg.default("species2:mut_rate", 10e-6)
    pg.default("species1:mut_rate", 10e-6)
    
    pg.add('test_geno','species2:genotype', [(1,1),(1,0),(0,1),(0,0)])
    pg.add('test_geno','species1:genotype', [(1,1),(1,0),(0,1),(0,0)])
    pg.add('test_mutrat','species2:mut_rate', [0, 10e-6, 10e-2,1])  
    pg.add('test_mutrat','species1:mut_rate', [0, 10e-6, 10e-2,1])


    pg.set_repeats(10)
    return pg

if __name__ == '__main__':
    print(len(externalinfluence()))
    externalinfluence().print()