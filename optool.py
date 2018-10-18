#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = 'zhihan'

from scipy.optimize import linprog
import matplotlib.pyplot as plt


# convert S to b

def cb(S):
    B = len(S)
    b = []
    num = 0
    while num < B:
        b.insert(num, 0)
        num += 1
    return b


# convert kcat,Et to Vmax

def kc2Vm(kc, E):
    Vm = []
    num = 0
    while num < len(kc):
        Vm.insert(num, kc[num] * E[num])
        num += 1
    return Vm


# convert bounds to list

def convert(bounds):
    B = []
    num = 0
    while num < len(bounds):
        B.insert(num, (0, bounds[num]))
        num += 1
    B = tuple(B)
    return B


# search for the enzyme

def op(c, A, kc, E, Enzyme):
    bounds = kc2Vm(kc, E)
    B0 = convert(bounds)
    EN = []
    Output = []
    b = cb(A)

    # calculate the initial output
    res = linprog(c, A_eq=A, b_eq=b, bounds=B0)
    f0 = -1 * res.fun
    fi = -1 * res.fun

    # set times of calculation circle
    r = len(bounds) - 1
    while r > 0:

        # calculate gradient
        i = len(bounds) - 1
        while i >= 0:
            bounds[i] += B0[i][1]
            B = convert(bounds)
            res = linprog(c, A_eq=A, b_eq=b, bounds=B)
            f = res.fun * -1

            if f > f0:
                EN.append(Enzyme[i])
                print('Enzyme:', Enzyme[i])
                Output.append(f / fi)
                print(' ')
                break
            else:
                bounds[i] -= B0[i][1]
            i -= 1

        # generate the output
        f0 = f

        # descent
        n = 100
        while n > 0:
            bounds[i] += B0[i][1]
            B = convert(bounds)
            res = linprog(c, A_eq=A, b_eq=b, bounds=B)
            f = res.fun * -1
            if f == f0:
                break
            n -= 1
        r -= 1

        # generate the output again
        f0 = f

    plt.plot(Output)
    plt.ylabel('Output Magnification')
    plt.xlabel(EN)
    plt.show()

    return bounds


# predict output

def prd(Km, c, A, E, bounds):
    # M-M equation
    v = bounds[0] * (E[0] / (Km[0] + E[0]))
    bounds[0] = v
    B = convert(bounds)
    b = cb(A)
    res = linprog(c, A_eq=A, b_eq=b, bounds=B)
    output = res.fun * -1
    print('Output=', output)
