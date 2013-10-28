import numpy as np
import matplotlib.pyplot as plt
import pdb

from groupy.order import *
from groupy.visualization import *

all_S2 = list()
for j in range(100):

    #atoms = np.random.randn(1000,3)
    atoms = np.ones((1000,3))
    directors = np.empty(shape=(100,3))

    for i, molecule in enumerate(np.split(atoms, 100)):
        I = calc_inertia_tensor(molecule)
        director = calc_director(I)
        directors[i] = director

    Q = calc_Q_tensor(directors)
    S2, director = calc_S2(Q)
    all_S2.append(S2)

plt.hist(all_S2)
plt.show()
