import numpy as np
from .insidetri import insidetri

def shadowAnaly(x, y, z, barC, delta, L_gw):
    """
    Checks for mesh elements shadowed from a direction by others.

    Parameters:
        x (np.ndarray): X coordinates of the vertices of the N faces (3xN).
        y (np.ndarray): Y coordinates of the vertices of the N faces (3xN).
        z (np.ndarray): Z coordinates of the vertices of the N faces (3xN).
        barC (np.ndarray): Barycentre of each mesh element (3xN).
        delta (np.ndarray): Angles between flow and surface normal for each mesh element (1xN).
        L_gw (np.ndarray): Rotation matrix from wind to geometric coordinate frame.

    Returns:
        np.ndarray: Indices of the panels that are shadowed.
    """
    pAw = np.dot(L_gw.T, np.array([x[0, :], y[0, :], z[0, :]]))
    pBw = np.dot(L_gw.T, np.array([x[1, :], y[1, :], z[1, :]]))
    pCw = np.dot(L_gw.T, np.array([x[2, :], y[2, :], z[2, :]]))
    barCw = np.dot(L_gw.T, barC)

    xW = np.vstack([pAw[0, :], pBw[0, :], pCw[0, :]])
    xWmax = np.max(xW, axis=0)
    xWmin = np.min(xW, axis=0)

    indB = np.where(delta * 180 / np.pi > 90.0001)[0]
    indF = np.where(delta * 180 / np.pi <= 90.0001)[0]

    minXwF = np.min(xWmin[indF])
    maxXwB = np.max(xWmin[indB])

    indFPot = np.where(xWmin[indF] - maxXwB < 1e-5)[0]
    indBPot = np.where(xWmax[indB] - minXwF > 1e-5)[0]

    yWC = np.vstack([pAw[1, indB[indBPot]], pBw[1, indB[indBPot]], pCw[1, indB[indBPot]]])
    zWC = np.vstack([pAw[2, indB[indBPot]], pBw[2, indB[indBPot]], pCw[2, indB[indBPot]]])

    shadPan = np.zeros(len(barCw[0]), dtype=int)
    tolB = 1e-5

    for i in indFPot:
        transYcord = yWC - barCw[1, indF[i]]
        yCoord = np.abs(np.sum(np.sign(transYcord), axis=0))
        yC_chang = np.where(yCoord < 3)[0]

        if len(yC_chang) == 0:
            continue

        transZcord = zWC[:, yC_chang] - barCw[2, indF[i]]
        zCoord = np.abs(np.sum(np.sign(transZcord), axis=0))
        zC_chang = np.where(zCoord < 3)[0]

        if len(zC_chang) == 0:
            continue

        for idx in zC_chang:
            p1 = np.array([transYcord[0, idx], transZcord[0, idx]])
            p2 = np.array([transYcord[1, idx], transZcord[1, idx]])
            p3 = np.array([transYcord[2, idx], transZcord[2, idx]])
            points = np.zeros((2, len(zC_chang)))

            inside = insidetri(p1.reshape(2, -1), p2.reshape(2, -1), p3.reshape(2, -1), points)
            if inside.size > 0:
                if np.any((barCw[0, indB[indBPot[yC_chang[idx]]]] - barCw[0, indF[i]]) > tolB):
                    shadPan[indF[i]] = 1

    return np.where(shadPan > 0)[0]
