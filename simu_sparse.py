import numpy as np
import pylab as plt
import pdb
from bisect import bisect
from scipy.stats import pareto
import os
import pandas as pd


def small_movie(simulation, nbr_steps):
    for i in range(nbr_steps):
        simulation.brownian_motion(1)
        name = "figure_" + str(i) + ".png"
        simulation.quick_display(fig_name=name)
    os.system("avconv -r 1 -i /tmp/figure_%d.png test_ctrw.mp4")


def binary_search(a, x):
    pos = bisect(a, x)          # find insertion position
    if pos == 0:
        return False
    else:
        return (True if a[pos - 1] == x else False)


def time_evolution_MSD(simulation, simu_param, nbr_steps):
    delta_t = 100
    cell = simulation(**simu_param)
    L_init = dict(cell.L)
    evolution = pd.DataFrame({"R2": []})
    for t in range(nbr_steps):
        cell.brownian_motion(delta_t)
        R2 = np.mean([(L_init[mol][0] - cell.L[mol][0])**2 +
                      (L_init[mol][1] - cell.L[mol][1])**2
                      for mol in cell.L.keys()])
        evolution.loc[t * delta_t] = R2
        if t % (nbr_steps / 100.) == 0:
            print (float(t) / nbr_steps) * 100, " % have been done."
    evolution.to_csv("MSD_obstacles/data_2.csv")


def time_evolution_MSD_CTRW(simulation, simu_param, nbr_steps):
    delta_t = 100
    cell = simulation(**simu_param)
    L_init = {i: cell.L[i]['coordinates'] for i in cell.L.keys()}
    evolution = pd.DataFrame({"R2": []})
    for t in range(nbr_steps):
        cell.brownian_motion(delta_t)
        R2 = np.mean([(L_init[mol][0] - cell.L[mol]['coordinates'][0])**2 +
                      (L_init[mol][1] - cell.L[mol]['coordinates'][1])**2
                      for mol in cell.L.keys()])
        evolution.loc[t * delta_t] = R2
        if t % (nbr_steps / 100.) == 0:
            print (float(t) / nbr_steps) * 100, " % have been done."
    evolution.to_csv("MSD_CTRW/data_3.csv")


def brown_equilibrium_evolution(simulation, nbr_steps):
    delta_t = 1
    C_t = {0: len(simulation.C)}
    L_t = {0: len(simulation.L)}
    R_t = {0: len(simulation.R)}
    for t in range(nbr_steps):
        simulation.brownian_motion(delta_t)
        C_t[t * delta_t] = len(simulation.C)
        L_t[t * delta_t] = len(simulation.L)
        R_t[t * delta_t] = len(simulation.R)
        if t % (nbr_steps / 100.) == 0:
            print (float(t) / nbr_steps) * \
                100, " % have been done. C = ", len(simulation.C)
    return C_t, L_t, R_t


def gen_obst_list(obst_density, patch):
    size = (patch[1][0] - patch[0][0] + 1) * (patch[1][1] - patch[0][1] + 1)
    length = int(np.sqrt(size))
    nbr_obst = int(size * obst_density)
    obst_list = [(a / length, a % length)
                 for a in np.random.choice(range(size), nbr_obst, replace=False)]
    obst_list.sort()
    obst_list = np.array(obst_list) + np.array(patch[0])
    return list([tuple(a) for a in obst_list])


class obstacles(object):

    def __init__(self, w, l_tot, r_tot, p_on, p_off, patch, D0, D1, obst_list, dt=1, dx=2):
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
        self.obst_list = list(obst_list)
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
                new_coor = tuple(self.generate_new_coordinates(coordinates))
                if new_coor not in molecules_coordinates.values():
                    if not binary_search(self.obst_list, new_coor):
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

    def quick_display(self, fig_name="figure.png"):
        plt.clf()
        ax = plt.axes(xlim=(-1, self.w + 1), ylim=(-1, self.w + 1))
        ax.set_xticks(range(-1, self.w + 1))
        ax.set_yticks(range(-1, self.w + 1))
        if len(self.L) != 0:
            line_coor_l = [np.array(self.L.values())[:, 0],
                           np.array(self.L.values())[:, 1]]
            ax.scatter(line_coor_l[0], line_coor_l[1], color='b',
                       marker='s', s=80)
        if len(self.R) != 0:
            line_coor_r = [np.array(self.R.values())[:, 0],
                           np.array(self.R.values())[:, 1]]
            ax.scatter(line_coor_r[0], line_coor_r[1], color='r',
                       marker='s', s=80)
        if len(self.C) != 0:
            line_coor_c = [np.array(self.C.values())[:, 0],
                           np.array(self.C.values())[:, 1]]
            ax.scatter(line_coor_c[0], line_coor_c[1], color='g',
                       marker='s', s=120)
        if len(self.obst_list) != 0:
            line_coor_o = [np.array(self.obst_list)[:, 0],
                           np.array(self.obst_list)[:, 1]]
            ax.scatter(line_coor_o[0], line_coor_o[1], color='k',
                       marker='s', s=80)
        plt.grid(which='both')
        plt.savefig("/tmp/" + fig_name)

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


class CTRW(object):

    def __init__(self, w, l_tot, r_tot, p_on, p_off, patch, D0, alpha, tau_c, obst_list=[], dt=1, dx=2):
        self.w = int(w)
        self.dt = float(dt)
        self.dx = float(dx)
        self.l_tot = int(l_tot)
        self.r_tot = int(r_tot)
        self.p_on = float(p_on)
        self.p_off = float(p_off)
        self.patch = list(patch)
        self.D0 = float(D0)
        self.alpha = float(alpha)
        self.tau_c = float(tau_c)
        self.obst_list = list(obst_list)
        self.C = {}
        self.L = {i: coor
                  for i, coor in enumerate(self.generate_init_pos_time(l_tot))}
        self.R = {i: coor
                  for i, coor in enumerate(self.generate_init_pos_time(r_tot))}

    def draw_rdn_powerlaw(self):
        drawed = pareto.rvs(self.alpha)
        drawed /= (self.dt**(-self.alpha) - self.tau_c **
                   (-self.alpha))**(1 / (self.alpha + 1))
        if drawed <= self.tau_c:
            return drawed
        else:
            return self.draw_rdn_powerlaw()

    def draw_rdn_marton(self):
        return self.dt + (1 - np.random.rand())**(-1 / (self.alpha - 1))

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

    def diffusion_step(self, molecules_coordinates, obst_list=[]):
        proba_factor = 4 * self.dt / (self.dx)**2
        for molecule in molecules_coordinates.keys():
            molecules_coordinates[molecule]['res_time'] -= self.dt
            if (molecules_coordinates[molecule]['res_time'] + self.dt <= 0
                    or not self.is_in_patch(molecules_coordinates[molecule]['coordinates'])):
                coordinates = molecules_coordinates[molecule]['coordinates']
                D_ij = self.D0
                if np.random.rand() < D_ij * proba_factor:
                    new_coor = self.generate_new_coordinates(coordinates)
                    if new_coor not in molecules_coordinates.values():
                        if not binary_search(self.obst_list, new_coor):
                            molecules_coordinates[molecule] = {'coordinates': tuple(new_coor),
                                                               'res_time': self.draw_rdn_powerlaw()}
        return molecules_coordinates

    def brownian_motion_one_step(self, obst_list=[]):
        # print "L \n", self.L
        # print "R \n", self.R
        # print "C \n", self.C
        self.L = self.diffusion_step(self.L, self.obst_list)
        self.R = self.diffusion_step(self.R, self.obst_list)
        self.C = self.diffusion_step(self.C, self.obst_list)
        delta_positive_C = []
        delta_negative_C = []
        L_coor = [val['coordinates'] for val in self.L.values()]
        R_coor = [val['coordinates'] for val in self.R.values()]
        for reactive in set(L_coor).intersection(R_coor):
            if np.random.rand() < self.p_on:
                delta_positive_C.append({'coordinates': reactive,
                                         'res_time': self.draw_rdn_powerlaw()})
        for product in self.C.values():
            if np.random.rand() < self.p_off:
                delta_negative_C.append(product)
        for new_prod in delta_positive_C:
            try:
                self.C[max(self.C.keys()) + 1] = new_prod
            except ValueError:
                self.C[0] = new_prod
            self.L = {key: value for key, value in self.L.items()
                      if value['coordinates'] != new_prod['coordinates']}
            self.R = {key: value for key, value in self.R.items()
                      if value['coordinates'] != new_prod['coordinates']}
        for new_reac in delta_negative_C:
            self.L[max(self.L.keys()) + 1] = new_reac
            self.R[max(self.R.keys()) + 1] = new_reac
            self.C = {key: value for key, value in self.C.items()
                      if value['coordinates'] != new_reac['coordinates']}

    def quick_display(self, fig_name="figure.png"):
        plt.clf()
        ax = plt.axes(xlim=(-1, self.w + 1), ylim=(-1, self.w + 1))
        ax.set_xticks(range(-1, self.w + 1))
        ax.set_yticks(range(-1, self.w + 1))
        if len(self.L) != 0:
            L_coor = [val['coordinates'] for val in self.L.values()]
            line_coor_l = [np.array(L_coor)[:, 0],
                           np.array(L_coor)[:, 1]]
            ax.scatter(line_coor_l[0], line_coor_l[1], color='b',
                       marker='s', s=80)
        if len(self.R) != 0:
            R_coor = [val['coordinates'] for val in self.R.values()]
            line_coor_r = [np.array(R_coor)[:, 0],
                           np.array(R_coor)[:, 1]]
            ax.scatter(line_coor_r[0], line_coor_r[1], color='r',
                       marker='s', s=80)
        if len(self.C) != 0:
            C_coor = [val['coordinates'] for val in self.C.values()]
            line_coor_c = [np.array(C_coor)[:, 0],
                           np.array(C_coor)[:, 1]]
            ax.scatter(line_coor_c[0], line_coor_c[1], color='g',
                       marker='s', s=120)
        if len(self.obst_list) != 0:
            None
        plt.grid(which='both')
        plt.savefig("/tmp/" + fig_name)

    def generate_init_pos_time(self, nbr_points):
        coor_list = []
        while len(coor_list) < nbr_points:
            x_new = np.random.choice(self.w)
            y_new = np.random.choice(self.w)
            new_mol = [x_new, y_new]
            if new_mol not in coor_list:
                coor_list.append({'coordinates': tuple(new_mol),
                                  'res_time': self.draw_rdn_powerlaw()
                                  })
        return coor_list

    def brownian_motion(self, nbr_steps):
        for t in range(nbr_steps):
            self.brownian_motion_one_step()


BIG_PATCH = [[0, 0], [799, 799]]
SMALL_PATCH = [[0, 10], [30, 30]]
VBIG_PATCH = [[0, 0], [4999, 4999]]

OBSTACLE_LIST = gen_obst_list(0.35, BIG_PATCH)

DISPLAY_PARAM = {'w': 800,
                 'l_tot': 300,
                 'r_tot': 70,
                 'p_on': 0.1,
                 'p_off': 10**(-3),
                 'patch': SMALL_PATCH,
                 'D0': 1,
                 'D1': 1,
                 'obst_list': OBSTACLE_LIST}

STANDARD_PARAM_CTRW = {'w': 800,
                       'l_tot': 1500,
                       'r_tot': 100,
                       'p_on': 0.1,
                       'p_off': 10**(-3),
                       'patch': BIG_PATCH,
                       'D0': 1,
                       'alpha': 0.8,
                       'tau_c': 5 * 10**3}

STANDARD_PARAM = {'w': 600,
                  'l_tot': 2500,
                  'r_tot': 100,
                  'p_on': 0.1,
                  'p_off': 10**(-3),
                  'patch': BIG_PATCH,
                  'D0': 0.1,
                  'D1': 0.1}

MSD_PARAM_OBSTACLES = {'w': 800,
                       'l_tot': 300,
                       'r_tot': 0,
                       'p_on': 0,
                       'p_off': 0,
                       'patch': BIG_PATCH,
                       'D0': 1,
                       'D1': 1,
                       'obst_list': OBSTACLE_LIST}

MSD_PARAM_CTRW = {'w': 800,
                  'l_tot': 300,
                  'r_tot': 0,
                  'p_on': 0,
                  'p_off': 0,
                  'patch': BIG_PATCH,
                  'D0': 1,
                  'alpha': 0.8,
                  'tau_c': 5 * 10**4}


def big_simu(simulation, model_param):
    simu = simulation(**model_param)
    C_t, L_t, R_t = brown_equilibrium_evolution(simu, 10**6)
    data = pd.DataFrame({'C': C_t, 'L': L_t, 'R': R_t})
    data.to_csv("long_time_study/data.csv")


def big_simu_CTRW(simulation, model_param):
    for tau in [50, 200, 800, 1600, 3200, 6400, 128000, 25600, 51200, 102400]:
        model_param["tau_c"] = tau
        simu = simulation(**model_param)
        C_t, L_t, R_t = brown_equilibrium_evolution(simu, 10**5)
        data = pd.DataFrame({'C': C_t, 'L': L_t, 'R': R_t})
        data.to_csv("equilibrium_vs_tau/data_" + str(tau) + "_2.csv")
        print "Simulation with tau_c = " + str(tau) + " done"


def res_time_effects(simulation, model_param):
    simu = simulation(**model_param)
    data = pd.DataFrame({'tau_moy': [], "card_tau_<0": []})
    for t in range(100000):
        # pdb.set_trace()
        R_mean = np.mean([val['res_time'] for val in simu.R.values()])
        L_mean = np.mean([val['res_time'] for val in simu.L.values()])
        C_mean = np.mean([val['res_time'] for val in simu.C.values()])
        if np.isnan(R_mean):
            R_mean = 0
        if np.isnan(L_mean):
            L_mean = 0
        if np.isnan(C_mean):
            C_mean = 0
        tau_moy = 1 / 3. * (R_mean + L_mean + C_mean)
        move_R = sum([val['res_time'] <= 0 for val in simu.R.values()])
        move_L = sum([val['res_time'] <= 0 for val in simu.L.values()])
        move_C = sum([val['res_time'] <= 0 for val in simu.C.values()])
        move_tot = move_R + move_L + move_C
        data.loc[t] = {'tau_moy': tau_moy, "card_tau_<0": move_tot}
        simu.brownian_motion(1)
        if t % (100000 / 100.) == 0:
            print (float(t) / 100000) * 100, " % have been done. "
    data.to_csv("res_time_study/res_time_simulation_tau=100000.csv")
    return data

# t1 = time.clock()
# OBSTACLES_LIST = gen_obst_list(0.35, VBIG_PATCH)
# t2 = time.clock()
# print "obstacles list generated in: ", t2 - t1, "s contains ",
# len(OBSTACLES_LIST), " elements"
