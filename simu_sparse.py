import numpy as np
import pylab as plt
import pdb

BIG_PATCH = [[200, 200], [599, 599]]
SMALL_PATCH = [[3, 3], [6, 6]]

STANDARD_PARAM = {'w': 600,
                  'l_tot': 2500,
                  'r_tot': 100,
                  'p_on': 0.1,
                  'p_off': 10**(-3),
                  'patch': BIG_PATCH,
                  'D0': 0.1,
                  'D1': 0.1}


def brown_equilibrium_evolution(simu_param, nbr_steps):
    delta_t = 10
    simulation = brownian_motion(**simu_param)
    C_t = {0: len(simulation.C)}
    L_t = {0: len(simulation.L)}
    R_t = {0: len(simulation.R)}
    for t in range(nbr_steps):
        simulation.brownian_motion(delta_t)
        C_t[t * delta_t] = len(simulation.C)
        L_t[t * delta_t] = len(simulation.L)
        R_t[t * delta_t] = len(simulation.R)
        if t % (nbr_steps / 100.) == 0:
            print (float(t) / nbr_steps) * 100, " % have been done. C = ", len(simulation.C)
    return C_t, L_t, R_t


class brownian_motion(object):

    def __init__(self, w, l_tot, r_tot, p_on, p_off, patch, D0, D1, dt=1, dx=2):
        self.w = int(w)
        self.dt = float(dt)
        self.dx = float(dx)
        self.l_tot = int(l_tot)
        self.r_tot = int(r_tot)
        self.p_on = float(p_on)
        self.p_off = float(p_off)
        self.patch = list(patch)
        self.D0 = float(D0)
        self.D1 = float(D1)
        self.C = {}
        self.L = {i: coor for i, coor in
                  enumerate(self.generate_init_postions(l_tot))}
        self.R = {i: coor for i, coor in
                  enumerate(self.generate_init_postions(r_tot))}

    def is_in_patch(self, coordinates):
        cond1 = coordinates[0] >= self.patch[0][0]
        cond2 = coordinates[0] <= self.patch[1][0]
        cond3 = coordinates[1] >= self.patch[0][1]
        cond4 = coordinates[1] <= self.patch[1][1]
        return cond1 and cond2 and cond3 and cond4

    def generate_new_coordinates(self, coordinates):
        neighbours = [[0, 1], [1, 0], [0, -1], [-1, 0]]
        move = neighbours[np.random.choice(range(4))]
        new_coordinates = list(coordinates)
        new_coordinates[0] += move[0]
        new_coordinates[0] = max(
            min(new_coordinates[0], self.w - 1), 0)
        new_coordinates[1] += move[1]
        new_coordinates[1] = max(
            min(new_coordinates[1], self.w - 1), 0)
        return new_coordinates

    def diffusion_step(self, molecules_coordinates):
        proba_factor = 4 * self.dt / (self.dx)**2
        for molecule in molecules_coordinates.keys():
            coordinates = molecules_coordinates[molecule]
            D_ij = self.D0
            if self.is_in_patch(coordinates):
                D_ij = self.D1
            if np.random.rand() < D_ij * proba_factor:
                new_coor = self.generate_new_coordinates(coordinates)
                if new_coor not in molecules_coordinates.values():
                    molecules_coordinates[molecule] = tuple(new_coor)
        return molecules_coordinates

    def brownian_motion_one_step(self):
        # print "L \n", self.L
        # print "R \n", self.R
        # print "C \n", self.C
        self.L = self.diffusion_step(self.L)
        self.R = self.diffusion_step(self.R)
        self.C = self.diffusion_step(self.C)
        delta_positive_C = []
        delta_negative_C = []
        for reactive in set(self.L.values()).intersection(self.R.values()):
            if np.random.rand() < self.p_on:
                delta_positive_C.append(reactive)
        for product in self.C.values():
            if np.random.rand() < self.p_off:
                delta_negative_C.append(product)
        for new_prod in delta_positive_C:
            try:
                self.C[max(self.C.keys()) + 1] = new_prod
            except ValueError:
                self.C[0] = new_prod
            self.L = {key: value for key, value in self.L.items()
                      if value != new_prod}
            self.R = {key: value for key, value in self.R.items()
                      if value != new_prod}
        for new_reac in delta_negative_C:
            self.L[max(self.L.keys()) + 1] = new_reac
            self.R[max(self.R.keys()) + 1] = new_reac
            self.C = {key: value for key, value in self.C.items()
                      if value != new_reac}

    def quick_display(self):
        line_coor_l = [np.array(self.L.values())[:, 0],
                       np.array(self.L.values())[:, 1]]
        plt.scatter(line_coor_l[0], line_coor_l[1], color='b')
        line_coor_r = [np.array(self.R.values())[:, 0],
                       np.array(self.R.values())[:, 1]]
        plt.scatter(line_coor_r[0], line_coor_r[1], color='r')
        try:
            line_coor_c = [np.array(self.C.values())[:, 0],
                           np.array(self.C.values())[:, 1]]
            plt.scatter(line_coor_c[0], line_coor_c[1], color='g')
        except IndexError:
            None
        plt.show()

    def generate_init_postions(self, nbr_points):
        coor_list = []
        while len(coor_list) < nbr_points:
            x_new = np.random.choice(self.w)
            y_new = np.random.choice(self.w)
            new_mol = [x_new, y_new]
            if new_mol not in coor_list:
                coor_list.append(tuple(new_mol))
        return coor_list

    def brownian_motion(self, nbr_steps):
        for t in range(nbr_steps):
            # if t % (nbr_steps / 10.) == 0:
            #     print (float(t) / nbr_steps) * 100, " % have been done."
            # self.quick_display()
            self.brownian_motion_one_step()


class obstacles(object):

    def __init__(self, w, l_tot, r_tot, p_on, p_off, patch, D0, dt=1, dx=2):
        self.w = int(w)
        self.dt = float(dt)
        self.dx = float(dx)
        self.l_tot = int(l_tot)
        self.r_tot = int(r_tot)
        self.p_on = float(p_on)
        self.p_off = float(p_off)
        self.patch = list(patch)
        self.D0 = float(D0)
        self.C = {}
        self.L = {i: coor for i, coor in
                  enumerate(self.generate_init_postions(l_tot))}
        self.R = {i: coor for i, coor in
                  enumerate(self.generate_init_postions(r_tot))}

    def generate_new_coordinates(self, coordinates):
        neighbours = [[0, 1], [1, 0], [0, -1], [-1, 0]]
        move = neighbours[np.random.choice(range(4))]
        new_coordinates = list(coordinates)
        new_coordinates[0] += move[0]
        new_coordinates[0] = max(
            min(new_coordinates[0], self.w - 1), 0)
        new_coordinates[1] += move[1]
        new_coordinates[1] = max(
            min(new_coordinates[1], self.w - 1), 0)
        return new_coordinates

    def diffusion_step(self, molecules_coordinates, obst_list):
        proba = (4 * self.dt / (self.dx)**2) * self.D0
        for molecule in molecules_coordinates.keys():
            coordinates = molecules_coordinates[molecule]
            if np.random.rand() < proba:
                new_coor = self.generate_new_coordinates(self, coordinates)
                if new_coor not in obst_list:
                    molecules_coordinates[molecule] = tuple(new_coor)
        return molecules_coordinates

    def brownian_motion_one_step(self, patch, model):
        self.L = self.diffusion_step(self.L, patch, model)
        self.R = self.diffusion_step(self.R, patch, model)
        self.C = self.diffusion_step(self.C, patch, model)
        delta_positive_C = []
        delta_negative_C = []
        for reactive in set(self.L).intersection(self.R):
            if np.random.rand() < self.p_on:
                delta_positive_C.append(reactive)
        for product in self.C.values():
            if np.random.rand() < self.p_off:
                delta_negative_C.append(product)
        for new_prod in delta_positive_C:
            self.C[len(self.C + 1)] = new_prod
            self.L = {key: value for key, value in self.L.items()
                      if value != new_prod}
            self.R = {key: value for key, value in self.R.items()
                      if value != new_prod}
        for new_reac in delta_negative_C:
            self.L[len(self.C + 1)] = new_reac
            self.R[len(self.C + 1)] = new_reac
            self.C = {key: value for key, value in self.C.items()
                      if value != new_reac}

    def brownian_motion(self, nbr_steps, patch, model):
        if model == 'obstacles':
            None
        for t in range(nbr_steps):
            # if t % (nbr_steps / 10.) == 0:
            #     print (float(t) / nbr_steps) * 100, " % have been done."
            self.brownian_motion_one_step(patch, model)

    def quick_display(self, molecules_coordinates):
        line_argument = [np.array(molecules_coordinates)[:, 0],
                         np.array(molecules_coordinates)[:, 1]]
        plt.scatter(line_argument[0], line_argument[1])

    def generate_init_postions(self, nbr_points):
        coor_list = []
        while len(coor_list) < nbr_points:
            x_new = np.random.choice(self.w)
            y_new = np.random.choice(self.w)
            new_mol = [x_new, y_new]
            if new_mol not in coor_list:
                coor_list.append(tuple(new_mol))
        return coor_list

    def generate_obst_list(self, patch):
        nbr_obst = self.models['obstacles']['obst_density']
        nbr_obst *= (patch[1][0] - patch[0][0] + 1) * \
            (patch[1][1] - patch[0][1] + 1)
        nbr_obst = int(nbr_obst)
        obst_list = []
        while len(obst_list) < nbr_obst:
            x_new = np.random.choice(range(patch[0][0], patch[1][0]))
            y_new = np.random.choice(range(patch[0][1], patch[1][1]))
            new_obstacle = [x_new, y_new]
            if new_obstacle not in obst_list:
                obst_list.append(tuple(new_obstacle))
        return obst_list

    #
    # def display_state(self):
    #     fig = plt.figure()
    #     ax1 = fig.add_subplot(211)
    #     ax1.matshow(self.L)
    #     ax2 = fig.add_subplot(212)
    #     ax2.matshow(self.R)
    #     ax3 = fig.add_subplot(222)
    #     ax3.matshow(self.C)
    #     plt.show()


def time_evolution_MSD(w, dt, dx, l_tot, r_tot,
                       nbr_steps, patch, models, model):
    cell = simulation(w, dt, dx, l_tot, r_tot, models, p_on=0, p_off=0)
    L_init = dict(cell.L)
    evolution = {}
    for t in range(nbr_steps):
        cell.brownian_motion(1, patch, model)
        # pdb.set_trace()
        R2 = np.mean([(L_init[mol][0] - cell.L[mol][0])**2 +
                      (L_init[mol][1] - cell.L[mol][1])**2
                      for mol in cell.L.keys()])
        evolution[t] = R2
        if t % (nbr_steps / 100.) == 0:
            print (float(t) / nbr_steps) * 100, " % have been done."
    return evolution
