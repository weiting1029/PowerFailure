import time
from functools import partial

import numpy as np
import scipy.linalg as la
import scipy.sparse.linalg as sla
from numpy.linalg import eig
from scipy.integrate import odeint
from scipy.sparse import diags
from os import fork, getpid
from disturbancesGnr import normaldisturbances
from violationChecking import globalcheck

import multiprocessing as mp

import pandas as pd


# import scipy.sparse.csr_matrix.multiply as sp_matmul


class PowerNetworkSolver(object):
    def __init__(self, int_theta, int_omega, A, L, ngnr, D, M, K, OMEGA):
        super(PowerNetworkSolver, self).__init__()
        self.int_omega = int_omega
        self.int_theta = int_theta
        self.A = A
        self.L = L
        self.ngnr = ngnr
        self.D = D
        self.M = M
        self.K = K
        self.OMEGA = OMEGA

    # @nb.njit()
    def kuramoto2nd(self, X, t):
        n = self.ngnr
        syn_omega = np.sum(self.OMEGA) / np.sum(self.D)
        OMEGA0 = np.reshape(self.OMEGA - self.D * syn_omega, (1, n))
        theta = X[0:n] - syn_omega * t
        dtheta = np.reshape(X[n:2 * n] - syn_omega, (1, n))
        matrix1 = np.repeat(np.reshape(theta, (1, n)), n, axis=0)
        matrix2 = matrix1 - np.transpose(matrix1)
        sinmatrix = np.sin(matrix2)
        domega = np.multiply(1 / self.M, -np.reshape(self.D, (1, n)) * dtheta + OMEGA0 - self.K * np.sum(
            np.multiply(sinmatrix, self.A), axis=0))
        return np.append(dtheta, domega)

    def solkuramoto(self, sol0, dt):
        return odeint(self.kuramoto2nd, sol0, dt)

    def analytical_solkuramoto(self, sol0, dt):
        n = self.ngnr
        sub_matrix1 = sla.inv(diags(self.M, format='csc').multiply(-sla.inv(diags(self.D, format='csc'))))
        sub_matrix2 = -self.K * sla.inv(diags(self.M, format='csc')) @ self.L
        sub_matrix3 = np.concatenate((sub_matrix1.toarray(), sub_matrix2), axis=1)
        sub_matrix4 = np.concatenate((np.eye(n), np.zeros((n, n))), axis=1)
        big_matrix = np.concatenate((sub_matrix3, sub_matrix4), axis=0)
        eigVals, eigVecs = eig(big_matrix)
        # Lambda = np.diag(eigVals)
        solutions = np.zeros([len(dt), 2 * n])
        for i in range(len(dt)):
            solutions[i, :] = eigVecs @ np.diag(np.exp(dt[i] * eigVals)) @ la.inv(eigVecs) @ sol0
            # print(np.exp(dt[i]*eigVals))
        #     solutions[i, :] = np.diag(np.exp(dt[i]*eigVals)) @ sol0
        return solutions

    def getDotOmega(self, theta, omega, nn):
        n = self.ngnr
        syn_omega = np.sum(self.OMEGA) / np.sum(self.D)
        OMEGA0 = np.reshape(self.OMEGA - self.D * syn_omega, (1, n))
        sol_domega = np.zeros((nn + 1, n))
        for i in range(nn + 1):
            dtheta = omega[i, :]
            matrix1 = np.repeat(np.reshape(theta[i, :], (1, n)), n, axis=0)
            matrix2 = matrix1 - np.transpose(matrix1)
            sinmatrix = np.sin(matrix2)
            sol_domega[i, :] = np.multiply(1 / self.M, -np.reshape(self.D, (1, n)) * dtheta + OMEGA0 - self.K * np.sum(
                np.multiply(sinmatrix, self.A), axis=0))
        return sol_domega

    # @nb.njit(nopython=True, parallel=True)
    def Simulation(self, check_times, thres, t, nn, disturbances):
        KK = len(disturbances)
        starttime = time.time()
        n = self.ngnr
        dt = np.linspace(0, t, nn + 1)
        vcheck_omega = np.zeros(n)
        vcheck_domega = np.zeros(n)
        vcheck_any = np.zeros(n)
        mcheck_omega = np.zeros((n, n))
        mcheck_domega = np.zeros((n, n))
        mcheck_any = np.zeros((n, n))

        # disturbances = dtb_gnr(n, KK, sigma)

        for k in range(KK):
            sol0 = np.pad(disturbances[k], (n, 0), 'constant', constant_values=(0, 0))
            vec_sol = self.solkuramoto(sol0, dt)  # get numerical tracks for theta and omega
            sol_domega = self.getDotOmega(vec_sol[:, :n], vec_sol[:, n:], nn)  # get numerical tracks for omega dot
            omega_domega = np.concatenate((vec_sol[:, n:], sol_domega), axis=1)  # create a vector of [omega,domega]
            all_check = globalcheck(omega_domega, check_times, thres, n, nn)
            vcheck_omega += all_check[0]
            vcheck_domega += all_check[1]
            vcheck_any += all_check[2]

            mcheck_omega += np.outer(all_check[0], all_check[0])
            mcheck_domega += np.outer(all_check[1], all_check[1])
            mcheck_any += np.outer(all_check[2], all_check[2])

        vcheck_omega = vcheck_omega / KK
        vcheck_domega = vcheck_domega / KK
        vcheck_any = vcheck_any / KK
        mcheck_omega = mcheck_omega / KK
        mcheck_domega = mcheck_domega / KK
        mcheck_any = mcheck_any / KK

        stoptime = time.time()
        totaltime = stoptime - starttime
        return {'vcheck_omega': vcheck_omega, 'vcheck_domega': vcheck_domega, 'vcheck_any': vcheck_any,
                'mcheck_omega': mcheck_omega, 'mcheck_domega': mcheck_domega, 'mcheck_any': mcheck_any,
                'totaltime': "{:2.2}sec".format(totaltime)}

    def analytical_Simulation(self, check_times, thres, t, nn, disturbances):
        KK = len(disturbances)
        starttime = time.time()
        n = self.ngnr
        dt = np.linspace(0, t, nn + 1)
        vcheck_omega = np.zeros(n)
        vcheck_domega = np.zeros(n)
        vcheck_any = np.zeros(n)
        mcheck_omega = np.zeros((n, n))
        mcheck_domega = np.zeros((n, n))
        mcheck_any = np.zeros((n, n))

        # disturbances = dtb_gnr(n, KK, sigma)

        for k in range(KK):
            sol0 = np.pad(disturbances[k], (0, n), 'constant', constant_values=(0, 0))
            vec_sol = self.analytical_solkuramoto(sol0, dt)  # get analytical tracks for theta and omega
            sol_domega = self.getDotOmega(vec_sol[:, n:], vec_sol[:, :n], nn)  # get numerical tracks for omega dot
            omega_domega = np.concatenate((vec_sol[:, :n], sol_domega), axis=1)  # create a vector of [omega,domega]
            all_check = globalcheck(omega_domega, check_times, thres, n, nn)
            vcheck_omega += all_check[0]
            vcheck_domega += all_check[1]
            vcheck_any += all_check[2]

            mcheck_omega += np.outer(all_check[0], all_check[0])
            mcheck_domega += np.outer(all_check[1], all_check[1])
            mcheck_any += np.outer(all_check[2], all_check[2])

        vcheck_omega = vcheck_omega / KK
        vcheck_domega = vcheck_domega / KK
        vcheck_any = vcheck_any / KK
        mcheck_omega = mcheck_omega / KK
        mcheck_domega = mcheck_domega / KK
        mcheck_any = mcheck_any / KK

        stoptime = time.time()
        totaltime = stoptime - starttime
        return {'vcheck_omega': vcheck_omega, 'vcheck_domega': vcheck_domega, 'vcheck_any': vcheck_any,
                'mcheck_omega': mcheck_omega, 'mcheck_domega': mcheck_domega, 'mcheck_any': mcheck_any,
                'totaltime': "{:2.2}sec".format(totaltime)}

    #
    def parallelized_Simulation(self, check_times, thres, t, nn, sigma, KK):
        n = self.ngnr
        disturbances = normaldisturbances(n, KK, sigma)
        # KK = len(disturbances)
        dt = np.linspace(0, t, nn + 1)
        vcheck_omega = np.zeros(n)
        vcheck_domega = np.zeros(n)
        vcheck_any = np.zeros(n)
        mcheck_omega = np.zeros((n, n))
        mcheck_domega = np.zeros((n, n))
        mcheck_any = np.zeros((n, n))

        for k in range(KK):
            sol0 = np.pad(disturbances[k], (0, n), 'constant', constant_values=(0, 0))
            vec_sol = self.analytical_solkuramoto(sol0, dt)  # get analytical tracks for theta and omega
            sol_domega = self.getDotOmega(vec_sol[:, n:], vec_sol[:, :n], nn)  # get numerical tracks for omega dot
            omega_domega = np.concatenate((vec_sol[:, :n], sol_domega), axis=1)  # create a vector of [omega,domega]
            all_check = globalcheck(omega_domega, check_times, thres, n, nn)
            vcheck_omega += all_check[0]
            vcheck_domega += all_check[1]
            vcheck_any += all_check[2]

            mcheck_omega += np.outer(all_check[0], all_check[0])
            mcheck_domega += np.outer(all_check[1], all_check[1])
            mcheck_any += np.outer(all_check[2], all_check[2])

        vcheck_omega = vcheck_omega / KK
        vcheck_domega = vcheck_domega / KK
        vcheck_any = vcheck_any / KK
        mcheck_omega = mcheck_omega / KK
        mcheck_domega = mcheck_domega / KK
        mcheck_any = mcheck_any / KK

        return vcheck_any

    def k_simulation(self, check_times, thres, t, nn, disturbances, k):
        n = self.ngnr
        dt = np.linspace(0, t, nn + 1)

        sol0 = np.pad(disturbances[k], (0, n), 'constant', constant_values=(0, 0))
        vec_sol = self.analytical_solkuramoto(sol0, dt)  # get analytical tracks for theta and omega
        sol_domega = self.getDotOmega(vec_sol[:, n:], vec_sol[:, :n], nn)  # get numerical tracks for omega dot
        omega_domega = np.concatenate((vec_sol[:, :n], sol_domega), axis=1)  # create a vector of [omega,domega]
        all_check = globalcheck(omega_domega, check_times, thres, n, nn)
        vcheck_omega = all_check[0]
        vcheck_domega = all_check[1]
        vcheck_any = all_check[2]

        mcheck_omega = np.outer(all_check[0], all_check[0])
        mcheck_domega = np.outer(all_check[1], all_check[1])
        mcheck_any = np.outer(all_check[2], all_check[2])

        df = pd.DataFrame({'RoCoF': vcheck_domega, 'AFV': vcheck_omega,
                           'AV': vcheck_any})
        # print("I'm process", getpid())
        return df

    #
    def parallelized_analytical_sml(self, check_times, thres, t, nn, disturbances, pool):
        # pool = mp.Pool(8)
        # print(mp.cpu_count() - 4)
        KK = len(disturbances)
        parallel_function = partial(self.k_simulation, check_times, thres, t, nn, disturbances)
        results = pool.map(parallel_function, [k for k in range(KK)])
        # pool.close()
        # pool.join()
        return results
