import numpy as np
from .insidetri import insidetri

def shadowAnaly(x, y, z, barC, delta, L_gw):
    """
    Identifies triangles shadowed by others (occlusion detection).

    Returns:
        np.ndarray: Indices of shadowed panels.
    """
    # Koordinatentransformation
    pAw = L_gw.T @ np.array([x[0, :], y[0, :], z[0, :]])
    pBw = L_gw.T @ np.array([x[1, :], y[1, :], z[1, :]])
    pCw = L_gw.T @ np.array([x[2, :], y[2, :], z[2, :]])
    barCw = L_gw.T @ barC

    xW = np.vstack([pAw[0], pBw[0], pCw[0]])
    xWmin = np.min(xW, axis=0)
    xWmax = np.max(xW, axis=0)

    indB = np.where(delta * 180 / np.pi > 90.0001)[0]
    indF = np.where(delta * 180 / np.pi <= 90.0001)[0]

    indFPot = np.where(xWmin[indF] - np.max(xWmin[indB]) < 1e-5)[0]
    indBPot = np.where(xWmax[indB] - np.min(xWmin[indF]) > 1e-5)[0]

    # Y/Z-Überschneidung vorbereiten
    yWC = np.vstack([pAw[1, indB[indBPot]], pBw[1, indB[indBPot]], pCw[1, indB[indBPot]]])
    zWC = np.vstack([pAw[2, indB[indBPot]], pBw[2, indB[indBPot]], pCw[2, indB[indBPot]]])

    shadPan = np.zeros(len(barCw[0]), dtype=int)
    tolB = 1e-5

    for i in indFPot:
        transY = yWC - barCw[1, indF[i]]
        yOverlap = np.abs(np.sum(np.sign(transY), axis=0)) < 3
        zOverlap = np.abs(np.sum(np.sign(zWC[:, yOverlap] - barCw[2, indF[i]]), axis=0)) < 3

        for idx in np.where(zOverlap)[0]:
            # Teste Punkt im projizierten Dreieck
            p1 = transY[0:2, idx]
            p2 = transY[1:3, idx]
            p3 = transY[2:4, idx]
            pt = np.zeros((2, 1))

            if insidetri(p1.reshape(2, 1), p2.reshape(2, 1), p3.reshape(2, 1), pt).any():
                if barCw[0, indB[indBPot[yOverlap][idx]]] - barCw[0, indF[i]] > tolB:
                    shadPan[indF[i]] = 1

    return np.where(shadPan > 0)[0]

