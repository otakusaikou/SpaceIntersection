#!/usr/bin/env python
# -*- coding: utf-8 -*-
from sympy import sin, cos, Matrix, symbols, lambdify
from optparse import OptionParser
import numpy as np
import pandas as pd


np.set_printoptions(suppress=True)  # Disable scientific notation for numpy


def getInit(xa, ya, EO, f):
    """Compute initial values of unknown parameters"""
    xa1, xa2 = xa.ravel()
    ya1, ya2 = ya.ravel()

    X1, Y1, Z1 = EO[0, :]
    X2, Y2, Z2 = EO[1, :]

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


def spaceIntersection(inputFile, s):
    """Perform a space intersection"""
    # For I.O.
    with open(inputFile) as fin:
        f = float(fin.readline())   # The focal length in mm

    # For E.O.
    # xp yp XL YL ZL O P K SigXL SigYL SigZL SigO SigP SigK
    data = pd.read_csv(
        inputFile,
        delimiter=' ',
        usecols=range(1, 15),
        names=[str(i) for i in range(14)],
        skiprows=1)

    EO, SigEO = np.hsplit(data.values[:, 2:], 2)

    # Convert from degrees to radians
    EO[:, 3:] = np.radians(EO[:, 3:])
    SigEO[:, 3:] = np.radians(SigEO[:, 3:])

    # For image points
    xa, ya = np.hsplit(data.values[:, :2], 2)

    # Compute initial values
    X0 = np.matrix(getInit(xa[:2], ya[:2], EO[:2, :3], f)).T

    # print "Initial Values:\n Param\tValue"
    # print "   XA\t%.6f" % X0[0, 0]
    # print "   YA\t%.6f" % X0[1, 0]
    # print "   ZA\t%.6f" % X0[2, 0]
    # print

    # Define variable for inerior orienration parameters
    IO = f, 0, 0

    # Define symbols
    EOs = symbols("XL YL ZL Omega Phi Kappa")   # E.O. parameters
    PTs = symbols("XA YA ZA")       # Object point coordinates
    pts = symbols("xa ya")          # Image coordinates

    # Define weight matrix
    err = SigEO.ravel()     # Error vector
    W = np.matrix(np.diag(s**2 / err**2))
    Q = W.I

    # List and linearize observation equations
    F = getEqn(IO, EOs, PTs, pts)
    JFx = F.jacobian(PTs)
    JFl = F.jacobian(EOs)       # Jacobian matrix for observables

    # Create lambda function objects
    FuncJFx = lambdify((EOs+PTs), JFx, 'numpy')
    FuncJFl = lambdify((EOs+PTs), JFl, 'numpy')
    FuncF = lambdify((EOs+PTs+pts), F, 'numpy')

    numPt = len(data)

    # Create observable array as argument of function objects
    l = np.zeros((numPt, 11))
    l[:, :6] = EO
    l[:, 6:9] = X0[:, :].T
    l[:, 9] += xa.ravel()
    l[:, 10] += ya.ravel()

    dX = np.ones(1)                              # Initial value for iteration

    # Iteration process
    lc = 0          # Loop count
    dRes = 1.       # Termination criteria
    res = 1.        # Initial value of residual
    while dRes > 10**-12 and lc < 20:
        # Compute coefficient matrix and constants matrix
        A = np.zeros((2 * numPt, len(err)))
        B = np.zeros((2 * numPt, 3))

        Ai = FuncJFl(*np.hsplit(l, 11)[:-2])
        Bi = FuncJFx(*np.hsplit(l, 11)[:-2])
        F0 = np.matrix(-FuncF(*np.hsplit(l, 11)).T.reshape(-1, 1))

        for i in range(numPt):
            A[2*i:2*(i+1), 6*i:6*(i+1)] = Ai[:, :, i].reshape(2, 6)
            B[2*i:2*(i+1), :] = Bi[:, :, i].reshape(2, 3)

        A = np.matrix(A)
        B = np.matrix(B)

        # Solve the unknown parameters
        AT = A.T.copy()
        Qe = (A * Q * AT)
        We = Qe.I
        N = (B.T * We * B)                  # Compute normal matrix
        t = (B.T * We * F0)                 # Compute t matrix
        dX = N.I * t                        # Compute unknown parameters
        V = Q * AT * We * (F0 - B * dX)     # Compute residual vector

        X0 += dX            # Update initial values
        l[:, 6:9] += dX[:, :].T

        # Update termination criteria
        if lc > 1:
            dRes = abs(((V.T * W * V)[0, 0]/res) - 1)
        res = (V.T * W * V)[0, 0]

        # Compute sigma0
        s0 = (res / (B.shape[0] - B.shape[1]))**0.5

        lc += 1

    # Compute other informations
    SigmaXX = s0**2 * N.I
    paramStd = np.sqrt(np.diag(SigmaXX))
    XA, YA, ZA = np.array(X0).flatten()
    # Print outputs/Ground Coordinate 
    print ("Object point coordinates:")
    print(XA)
    print(YA)
    print(ZA)
    
     return XA, YA, ZA

s = 0.005
input = "./input4.txt"
spaceIntersection(input, s)

    # Output results
    # print "Object point coordinates:"
    # print (" %9s %11s %11s") % ("Parameter", "Value", "Std.")
    # print " %-10s %11.6f %11.6f" % ("XL", XA, paramStd[0])
    # print " %-10s %11.6f %11.6f" % ("YL", YA, paramStd[1])
    # print " %-10s %11.6f %11.6f" % ("ZL", ZA, paramStd[2])
    # print "\nSigma0 : %.6f" % s0


""" def main():
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
        options.input = "./input.txt"

    if not options.s:
        options.s = 0.005

    spaceIntersection(inputFile=options.input, s=options.s)

    return 0


if __name__ == '__main__':
    main()
"""
