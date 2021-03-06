from __future__ import print_function

from copy import deepcopy
from random import seed, randint, shuffle

import numpy as np

from groupy.gbb import Gbb
from groupy.box import Box
from groupy.lattice import Lattice
import pdb


class Bilayer():
    """Class for initializing a bilayer. 

    """

    def __init__(self, lipids, n_x=10, n_y=10, 
                 area_per_lipid=1.0, solvent=None, solvent_per_lipid=None,
                 random_seed=12345,
                 mirror=True, solvent_density=1.0/1.661,
                 lipid_positions=None):
        """Constructor. Builds a bilayer with several tunable parameters.

        Args
        ----
            lipids : list((groupy.gbb, fraction))
                The lipid prototypes to go into the bilayer
            n_x : int
                The number of lipids in x
            n_y : int
                The number of lipids in y
            area_per_lipid : float
                area per lipid, determines box cross-sectional area
            solvent : groupy.gbb
                Solvent molecule
            solvent_per_lipid : int
                Number of solvent molecules per lipid
            random_seed : int
                Initialize random number generator
            mirror : bool
                If True, make the leaflets mirror images
            solvent_density : float
                Mass density of solvent to use
            lipid_positions : list(list(int))
                Positions of each component, used to phase-separated systems
        """

        # set constants from constructor call
        self.lipids = lipids
        self.solvent = solvent
        self.apl = area_per_lipid
        self.n_x = n_x
        self.n_y = n_y
        self.mirror = mirror
        self.random_seed = random_seed
        self.n_solvent_per_lipid = solvent_per_lipid
        self.molecules = []
        self.solvent_density = solvent_density
        self.lipid_positions = lipid_positions

        # parse the spacing_z and shift lipids as necessary
        for i, lipid in enumerate(lipids):
            lipid[0].shift_com_to_origin()
            lipid[0].translate(xyz=[0.0, 0.0, lipid[2]])         

        # a few calculations to figure things out
        self.n_lipids_per_layer = self.n_x * self.n_y
        self.n_solvent_per_layer = self.n_lipids_per_layer * self.n_solvent_per_lipid
        self.area = n_x * n_y * self.apl
        self.box = Box(lengths=np.array(
            [np.sqrt(self.area), 
             np.sqrt(self.area),
             1.0]))
        self.mask = Lattice()
        self.mask.grid_mask_2d(n_x, n_y, box=self.box)

        self._n_each_lipid_per_layer = list()

        # safety checks
        self.check_fractions()
        if lipid_positions:
            assert(len(lipids) == len(lipid_positions))
            for pos in range(n_x * n_y):
                found = False
                for comp_pos in lipid_positions:
                    if pos in comp_pos:
                        found = True
                        break
                if found == False:
                    print('Not all lattice positions accounted for')
                    exit()

        # for finding location of water boxes
        self.top_bound = -np.inf
        self.bottom_bound = np.inf

        # assemble the lipid layers
        seed(self.random_seed)   
        if lipid_positions:
            top_layer, new_lipid_labels = self.create_layer(lipid_labels=lipid_positions,
                random_mixing=False)
            bottom_layer, bottom_lipid_labels = self.create_layer(
                lipid_labels=lipid_positions, random_mixing=False, 
                flip_orientation=True)
        else: 
            top_layer, top_lipid_labels = self.create_layer()
            if self.mirror == True:   # ugly
                bottom_layer, bottom_lipid_labels = self.create_layer(
                        lipid_labels=new_lipid_labels,
                        flip_orientation=True)
            else:
                bottom_layer, bottom_lipid_labels = self.create_layer(
                        flip_orientation=True)
        for lipid in top_layer:
            self.molecules.append(lipid)
        for lipid in bottom_layer:
            self.molecules.append(lipid)
        if self.n_solvent_per_lipid > 0:
            solvents = self.add_solvent()
            for solvent in solvents:
                self.molecules.append(solvent)

    def check_fractions(self):
        frac_sum = 0
        for lipid in self.lipids:
            frac_sum += lipid[1]
        assert frac_sum == 1.0, 'Bilayer builder error: Lipid fractions do not add up to 1.'

    def n_each_lipid_per_layer(self):
        if self._n_each_lipid_per_layer:
            return self._n_each_lipid_per_layer

        self._n_each_lipid_per_layer = []
        for lipid in self.lipids[:-1]:
            self._n_each_lipid_per_layer.append(
                    int(round(lipid[1] * self.n_lipids_per_layer)))
        # TODO: give warning if frac*n different than actual
        # rounding errors may make this off by 1, so just do total - whats_been_added
        self._n_each_lipid_per_layer.append(
                self.n_lipids_per_layer - sum(self._n_each_lipid_per_layer))
        assert len(self._n_each_lipid_per_layer) == len(self.lipids)
        return self._n_each_lipid_per_layer

    @property
    def n_each_lipid_per_layer(self):
        import pdb
        if self._n_each_lipid_per_layer:
            return self._n_each_lipid_per_layer

        self._n_each_lipid_per_layer = []
        for lipid in self.lipids[:-1]:
            self._n_each_lipid_per_layer.append(
                    int(round(lipid[1] * self.n_lipids_per_layer)))
        # TODO: give warning if frac*n different than actual
        # rounding errors may make this off by 1, so just do total - whats_been_added
        self._n_each_lipid_per_layer.append(
                self.n_lipids_per_layer - sum(self._n_each_lipid_per_layer))
        assert len(self._n_each_lipid_per_layer) == len(self.lipids)
        return self._n_each_lipid_per_layer

    def create_layer(self, lipid_labels=None, flip_orientation=False, 
                     random_mixing=True):
        """
        Parameters
        ----------
        lipid_labels : list(int)
            An integer label for each lipid
        flip_orientation : bool
            If True, rotate lipids aboux x-axis (to create an opposite layer)
        random_mixing : bool
            If True, place the lipids randomly on a lattice
        """
        layer = []
        if not lipid_labels:
            lipid_labels = range(self.n_lipids_per_layer)
            shuffle(lipid_labels)
        lipids_placed = 0
        for i, n_of_lipid_type in enumerate(self.n_each_lipid_per_layer):
            if self.lipid_positions:
                assert(len(self.lipid_positions[i]) == n_of_lipid_type)
            current_type = self.lipids[i][0]
            for j, n_this_lipid_type in enumerate(range(n_of_lipid_type)):
                new_lipid = deepcopy(current_type)
                if flip_orientation == True:
                    new_lipid.rotate(angles=np.array([np.pi, 0.0, 0.0]))

                # Move to point on mask
                if random_mixing == True:
                    random_index = lipid_labels[lipids_placed]
                    position = self.mask.points[random_index]
                    new_lipid.translate(position)
                else:
                    position = self.mask.points[self.lipid_positions[i][j]]
                    new_lipid.translate(position)

                # see if lipid bounds change
                if np.amax(new_lipid.xyz[:, 2]) > self.top_bound:
                    self.top_bound = np.amax(new_lipid.xyz[:, 2])
                if np.amin(new_lipid.xyz[:, 2]) < self.bottom_bound:
                    self.bottom_bound = np.amin(new_lipid.xyz[:, 2])
                    
                layer.append(new_lipid)
                lipids_placed += 1
        return layer, lipid_labels

    def add_solvent(self):
        # just do this in a very basic manner for now
        # TODO: add in smart solvate like mBuild
        # first find the xy area of the box
        solvent_list = []
        box_area = self.box.lengths[0] * self.box.lengths[1]
        MW_solvent = np.sum(self.solvent.masses)
        water_box_z = self.n_solvent_per_layer * MW_solvent / (
                box_area * self.solvent_density)
        top_box_min = self.top_bound + (0.5 * float(self.n_x) / self.box.lengths[0])
        bottom_box_max = self.bottom_bound - (0.5 * float(self.n_x) / self.box.lengths[0])
        top_water_box = Box(mins=np.array(
                           [self.box.mins[0],
                            self.box.mins[1],
                            top_box_min]), 
                           maxs=np.array(
                               [self.box.maxs[0],
                                self.box.maxs[1],
                                top_box_min + water_box_z]))
        bottom_water_box = Box(
                mins=np.array(
                    [self.box.mins[0],
                     self.box.mins[1],
                     bottom_box_max - water_box_z]),
                maxs=np.array(
                    [self.box.maxs[0],
                     self.box.maxs[1],
                     bottom_box_max]))
        self.box = Box(mins=np.array(
            [self.box.mins[0],
             self.box.mins[1],
             bottom_water_box.mins[2]]),
            maxs = np.array(
            [self.box.maxs[0],
             self.box.maxs[1],
             top_water_box.maxs[2]]))
        n_water_x = int((top_water_box.lengths[0]**2.0 * self.n_solvent_per_layer / (
                top_water_box.lengths[2] * top_water_box.lengths[1]))**(1.0/3.0))
        n_water_y = int((top_water_box.lengths[1] * self.n_solvent_per_layer / (
            top_water_box.lengths[2] * n_water_x))**0.5)
        # add 1 for rounding purposes, stop adding when n_added = n_solvent_per_layer
        # TODO: will this always work? test to make sure
        n_water_z = 1 + int(self.n_solvent_per_layer / (n_water_x * n_water_y))
        top_water_points = Lattice()
        top_water_points.grid_mask_3d(n_water_x, n_water_y, n_water_z,
                box=top_water_box)
        bottom_water_points = Lattice()
        bottom_water_points.grid_mask_3d(n_water_x, n_water_y, n_water_z,
                box=bottom_water_box)
        self.solvent.shift_com_to_origin()
        for point in top_water_points.points[:self.n_solvent_per_layer]:
            t = deepcopy(self.solvent)
            t.translate(point)
            solvent_list.append(t)
        for point in bottom_water_points.points[:self.n_solvent_per_layer]:
            t = deepcopy(self.solvent)
            t.translate(point)
            solvent_list.append(t)
        return solvent_list
