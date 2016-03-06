#!/usr/bin/env python
# -*- coding: utf-8 -*-
from sympy import sin, cos, Matrix, symbols, lambdify
from optparse import OptionParser
from math import radians as rad
import numpy as np


np.set_printoptions(suppress=True)  # Disable scientific notation for numpy


def getInit(xa1, ya1, xa2, ya2, EO1, EO2, f):
    """Compute initial values of unknown parameters"""
    X1, Y1, Z1 = EO1
    X2, Y2, Z2 = EO2

    B = np.sqrt((X2 - X1)**2 + (Y2 - Y1)**2)    # The baseline
    pa = ya1 - ya2                              # The parallax

    H = (Z1 + Z2) / 2

    # Compute arbitrary horizontal coordinates with formula 8-5~8-7
    XA = B * (xa1 / pa)
    YA = B * (ya1 / pa)
    ZA = H - (B * f) / pa

    # Compute the transformation parameters between
    # arbitrary and true object coordinate system
    a = np.cos(np.arctan2((X2 - X1), (Y2 - Y1)))
    b = np.sin(np.arctan2((X2 - X1), (Y2 - Y1)))
    Tx = X1
    Ty = Y1

    # Transform the horizontal coordinates of arbitrary object point
    # and use the result as initial values
    XA2 = a * XA - b * YA + Tx
    YA2 = a * YA + b * XA + Ty

    return XA2, YA2, ZA


def getM(Omega, Phi, Kappa):
    """Compute rotation matrix M"""
    M = np.matrix([
        [
            cos(Phi)*cos(Kappa),
            sin(Omega)*sin(Phi)*cos(Kappa) + cos(Omega)*sin(Kappa),
            -cos(Omega)*sin(Phi)*cos(Kappa) + sin(Omega)*sin(Kappa)],
        [
            -cos(Phi)*sin(Kappa),
            -sin(Omega)*sin(Phi)*sin(Kappa) + cos(Omega)*cos(Kappa),
            cos(Omega)*sin(Phi)*sin(Kappa) + sin(Omega)*cos(Kappa)],
        [
            sin(Phi),
            -sin(Omega)*cos(Phi),
            cos(Omega)*cos(Phi)]
        ])

    return M


def getEqn(IO, EO, PT, pt):
    """List observation equations"""
    f, xo, yo = IO
    XL, YL, ZL, Omega, Phi, Kappa = EO
    XA, YA, ZA = PT
    xa, ya = pt

    M = getM(Omega, Phi, Kappa)

    r = M[0, 0] * (XA - XL) + M[0, 1] * (YA - YL) + M[0, 2] * (ZA - ZL)
    s = M[1, 0] * (XA - XL) + M[1, 1] * (YA - YL) + M[1, 2] * (ZA - ZL)
    q = M[2, 0] * (XA - XL) + M[2, 1] * (YA - YL) + M[2, 2] * (ZA - ZL)

    F = Matrix([xa - xo + f * (r / q), ya - yo + f * (s / q)])
    return F


def spaceIntersection(inputFile="input.txt", s=rad(5)):
    """Perform a space intersection"""
    # Define symbols
    EO = symbols("XL YL ZL Omega Phi Kappa")  # Exterior orienration parameters
    PT = symbols("XA YA ZA")    # Object point coordinates
    pt = symbols("xa ya")       # Image coordinates

    # Read observables from txt file
    fin = open(inputFile)
    lines = fin.readlines()
    fin.close()

    f = float(lines[0])     # The focal length in mm
    EO1, EO2 = np.vsplit(np.array(map(
        lambda x: x.split()[1:], lines[1:3])).astype(np.double), 2)

    EO1, SigEO1 = map(lambda x: x.flatten(), np.hsplit(EO1, 2))
    EO2, SigEO2 = map(lambda x: x.flatten(), np.hsplit(EO2, 2))

    # Convert the angles from degree to radian
    EO1[3:] = map(lambda x: rad(x), EO1[3:])
    EO2[3:] = map(lambda x: rad(x), EO2[3:])

    xa1, ya1, xa2, ya2, = map(lambda x: float(x), lines[3].split()[1:])

    # Compute initial values
    X0 = np.matrix(getInit(xa1, ya1, xa2, ya2, EO1[:3], EO2[:3], f)).T

    print "Initial Values:\n Param\tValue"
    print "   XA\t%.6f" % X0[0, 0]
    print "   YA\t%.6f" % X0[1, 0]
    print "   ZA\t%.6f" % X0[2, 0]
    print

    # Define variable for inerior orienration parameters
    IO = f, 0, 0

    # List and linearize observation equations
    F = getEqn(IO, EO, PT, pt)
    JFx = F.jacobian(PT)
    JFl = F.jacobian(EO)    # Jacobian matrix for observables

    # Create lambda function objects
    FuncJFx = lambdify((EO+PT), JFx)
    FuncJFl = lambdify((EO+PT), JFl)
    FuncF = lambdify((EO+PT+pt), F)

    # Define weight matrix
    err = np.append(SigEO1, SigEO2)     # Error vector
    W = np.matrix(np.diag(s**2 / err**2))

    dX = np.ones(1)                              # Initial value for iteration

    # Iteration process
    while abs(dX.sum()) > 10**-12:
        # Compute coefficient matrix and constants matrix
        A = np.matrix(np.zeros((4, 12))).astype(np.double)
        B = np.matrix(np.zeros((4, 3))).astype(np.double)
        f = np.matrix(np.zeros((4, 1))).astype(np.double)

        # Deifne array of constants
        val1 = np.concatenate(
            (EO1, np.array(X0).flatten(), np.array([xa1, ya1])))
        val2 = np.concatenate(
            (EO2, np.array(X0).flatten(), np.array([xa2, ya2])))

        A[:2, :6] = FuncJFl(*val1[:-2])
        A[2:, 6:] = FuncJFl(*val2[:-2])
        B[:2] = FuncJFx(*val1[:-2])
        B[2:] = FuncJFx(*val2[:-2])
        f[:2] = -FuncF(*val1)
        f[2:] = -FuncF(*val2)

        # Solve the unknown parameters
        Qe = (A * W * A.T).I
        We = Qe.I
        N = (B.T * We * B)                  # Compute normal matrix
        t = (B.T * We * f)                  # Compute t matrix
        dX = N.I * t                        # Compute unknown parameters
        V = W.I * A.T * We * (f - B * dX)   # Compute residual vector

        X0 += dX            # Update initial values

        # Compute error of unit weight
        res = (V.T * W * V)[0, 0]
        s0 = (res / (B.shape[0] - B.shape[1]))**0.5

    # Compute other informations
    SigmaXX = s0**2 * N.I
    paramStd = np.sqrt(np.diag(SigmaXX))
    XA, YA, ZA = np.array(X0).flatten()

    # Output results
    print "Object point coordinates:"
    print (" %9s %11s %11s") % ("Parameter", "Value", "Std.")
    print " %-10s %11.6f %11.6f" % ("XL", XA, paramStd[0])
    print " %-10s %11.6f %11.6f" % ("YL", YA, paramStd[1])
    print " %-10s %11.6f %11.6f" % ("ZL", ZA, paramStd[2])
    print "\nError of unit weight : %.6f" % s0


def main():
    parser = OptionParser(usage="%prog [-i] [-s]", version="%prog 0.1")

    # Define options
    parser.add_option(
        "-i", "--input",
        help="read input data from FILE, the default value is \"input.txt\"",
        metavar="FILE")
    parser.add_option(
        "-s", "--sigma",
        type="float",
        dest="s",
        help="define a priori error, the default value is 0.001 m",
        metavar="N")

    # Instruct optparse object
    (options, args) = parser.parse_args()

    # Define default values if there are nothing given by the user
    if not options.input:
        options.input = "input.txt"

    if not options.s:
        options.s = 0.001

    spaceIntersection(inputFile=options.input, s=options.s)

    return 0


if __name__ == '__main__':
    main()
