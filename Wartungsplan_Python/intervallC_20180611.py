import itertools
import numpy as np

## General used variables
combintervall = None
combCP = None
combCELmax = None
combmplan = None


def calcCost(numcomp, q10, n, mplanred):
    global combmplan

    #% Wartungsschritte und -kosten pro Zeitschritt
    combmplan = np.zeros((np.power(int(np.max(q10)), numcomp), n, numcomp))
    combCGES = np.zeros((np.power(int(np.max(q10)), numcomp), n))
    CGES = np.zeros((np.power(int(np.max(q10)), numcomp), 1))


    #combmplan = np.zeros((numcomp * int(np.max(q10)), n, numcomp))
    #combCGES = np.zeros((numcomp * int(np.max(q10)), n))
    #CGES = np.zeros((numcomp * int(np.max(q10)), 1))


    i4 = 0
    icost = np.zeros((numcomp, 1))

    if numcomp == 3:
        for bb in range(len(mplanred)):
            for cc in range(len(mplanred)):
                for dd in range(len(mplanred)):
                    for icost4 in range(len(mplanred)):
                        for ee in range(n):
                            vector = [mplanred[0][bb][ee], mplanred[1][cc][ee], mplanred[2][dd][ee]]
                            combmplan[i4][ee] = vector
                            cell_cisize = len(combintervall)
                            for i6 in range(cell_cisize):
                                check = [np.sum(a == combintervall[1][i6]) for a in combmplan[i4][ee]]  #  ismember(combmplan{i4,ee},cell_ci{1,i6})
                                if sum(check) == numcomp:
                                    combCGES[i4][ee] = sum(combCP[:][i6]) + combCELmax[1][i6]

                                check = 0

                        #% Kosten der Wartungsstrategie = Summe der Kosten f�r
                        #% Zeitschritte
                        CGES[i4][1] = sum(combCGES[i4][:])
                        i4 = i4+1

    elif numcomp == 4:
        for bb in range(len(mplanred[0])):
            for cc in range(len(mplanred[0])):
                for dd in range(len(mplanred[0])):
                    for icost4 in range(len(mplanred[0])):
                        for ee in range(n):
                            vector = [mplanred[0][bb][ee], mplanred[1][cc][ee], mplanred[2][dd][ee], mplanred[3][icost4][ee]]
                            for ncomp in range(numcomp):
                                combmplan[i4][ee][ncomp] = vector[ncomp]
                            cell_cisize = len(combintervall)
                            check = 0#[0, 0]
                            for i6 in range(cell_cisize):
                                if all(combmplan[i4][ee] == combintervall[i6]) == True:
                                    combCGES[i4][ee] = np.sum(combCP[:][i6]) + combCELmax[0][i6]
                                    check = 0

                    # Kosten der Wartungsstrategie = Summe der Kosten für
                    # Zeitschritte
                        CGES[i4][0] = np.sum(combCGES[i4][:])
                        i4 = i4+1
    return CGES



def calcInterval(CP, CEL):

    # Wartungsfälle durch Kombination für Zeitschritt definieren
    # Alle möglichen Wartungsaktivitäten pro Zeitschritt
    a = [[1, 0], [1, 0]]

    global combintervall
    global combCP
    global combCELmax


    combintervall = list(itertools.product(*a, repeat=2))
    asize = [len(combintervall[0]), len(combintervall)]
    #cell_ci = [1][asize[1][2]]
    #cell_ci = np.zeros((asize[0], asize[1]))
    for ci in range(len(combintervall)):
        combintervall[ci] = list(reversed(combintervall[ci]))

    # Kosten für Komponentenwartung
    combCP = np.zeros((asize[1], asize[0]))
    combCEL = np.zeros((asize[1], asize[0]))
    for a1 in range(asize[1]):
        for a2 in range(asize[0]):
            if combintervall[a1][a2] == 1:
                combCP[a1][a2] = CP[a2]
                combCEL[a1][a2] = CEL[a2]



    # Arrays mit präventiven Wartungskosten sowie Kosten für die Ersatzleistung
    combCELsize = len(combCEL)
    combCELmax = np.zeros((1, combCELsize))
    for c4 in range(len(combCEL)):
        if np.sum(combCEL[:][c4]) == 0:
            combCELmax[0][c4] = 0
        elif sum(combCEL[:][c4]) == 1:
            for c5 in range(combCELsize):
                if combCEL[1][c5] > 0:
                    combCELmax[0][c4] = combCEL[1][c5]
                    break
                elif combCEL[1][c5] > 0:
                    combCELmax[0][c4] = combCEL[2][c5]
                    break
                elif combCEL[2][c5] > 0:
                    combCELmax[0][c4] = combCEL[3][c5]
                    break
        elif np.sum(combCEL[:][c4]) > 1:
            combCELmax[0][c4] = max(combCEL[:][c4])
