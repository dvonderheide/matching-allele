from make_parameters import Paramgroup as PG
import numpy as np
def externalinfluence():
    pg = PG()
    pg.default("general:runtime", 20)
    
    pg.add('geno1','species2:genotype', [(1,1),(1,0),(0,1),(0,0)])
    pg.add('geno2','species1:genotype', [(1,1),(1,0),(0,1),(0,0)])
    pg.add('mut_rate1','species2:mut_rate', [0, 10e-6, 10e-2,1])  
    pg.add('mut_rate2','species1:mut_rate', [0, 10e-6, 10e-2,1])


    pg.set_repeats(10)
    return pg

if __name__ == '__main__':
    print(len(externalinfluence()))
    externalinfluence().print()