import numpy as np
import matplotlib.pyplot as plt
import itertools
import Wartungsplan_Python.intervallC_20180611 as intervallC
import Wartungsplan_Python.input_nonlinear as input

n = 8               # Anzahl Zeitschritte
skalierung = 1000   # Betriebsstunden pro Zeitschritt

X0 = np.linspace(0, (n-1) * skalierung, n)


############################################
# # # # #Konstanten

# Entscheidungsvariable, wann gewartetet werden soll
# Ist Maß für die Risikobereitschaft des Entscheiders
drisk = 1 * input.CP


###########################################
# # # # # Eigene Funktionen

def mat2cell_porting(array):
    vectorAsArray = np.zeros((array[0].size, len(array)))

    for index in range(0, array[0].size):

        vectorAsArray[index][0] = array[0][index]
        vectorAsArray[index][1] = array[1][index]

    return vectorAsArray


def calcG(x, i):

    y = 1 - np.exp(-x/input.T[i])

    return y


##########################################
# # # # # Überführungsversuch

# mplan = np.zeros((numcomp, n * n, 2))
mplan = np.zeros((input.numcomp, np.power(2, n), n))
m = None
for j in range(input.numcomp):
    L = np.ones((2, n))
    L[1] = np.zeros(n)  # example input
    L2 = mat2cell_porting(L)
    if j == 0:
        mplan[j] = list(itertools.product(*L2))
        # TODO: Was wurde hier gemacht?
        # mplan=shiftdim(combvec(L2{:})',-1);
        # : mplan = permute(mplan,[2 3 1]);
        m = len(mplan[0])
    elif j >= 1:
        mplan[j] = mplan[0]

for j in range(input.numcomp):
    for i in range(len(mplan[0])):
        mplan[j][i] = list(reversed(mplan[j][i]))


# Ausfallwahrscheinlichkeiten G und Erwartungswert der Kosten EWC
G = np.zeros((input.numcomp, m, n))
EWC = np.zeros((input.numcomp, m, n))
EWCGES = np.zeros((m, n))
EWCEL = np.zeros((input.numcomp, m, n))
q10 = np.zeros((1, input.numcomp))
A = np.zeros((input.numcomp, m, n))


for jj in range(input.numcomp):
    for ll in range(m):
        for mm in range(n):
            # Wenn erster Zeitschritt Wartung, Ausfallwahrscheinlichkeit = 0
            if mm == 0:
                if mplan[jj][ll][mm] == 1:
                    G[jj][ll][mm] = 0
                    EWC[jj][ll][mm] = 0
                    EWCEL[jj][ll][mm] = G[jj][ll][mm] * input.CEL[jj]

                elif mplan[jj][ll][mm] == 0:
                    for b1 in range(n):
                        mplan[jj][ll][b1] = np.Inf

            # Ab zweitem Zeitschritt
            elif mm >= 1:
                if mplan[jj][ll][mm] == 1:
                    G[jj][ll][mm] = 0
                    EWC[jj][ll][mm] = 0
                    EWCEL[jj][ll][mm] = 0
                elif mplan[jj][ll][mm] == 0:
                    # Differenz der Wahrscheinlichkeit im Intervall addieren
                    G[jj][ll][mm] = G[jj][ll][mm-1] + calcG((mm-(mm-1))*skalierung, jj)

                    EWC[jj][ll][mm] = G[jj][ll][mm] * input.CR[jj]
                    EWCEL[jj][ll][mm] = G[jj][ll][mm] * input.CEL[jj]

            # Wenn Ausfallwahrscheinlichkeit innerhalb der
            # Simulation >= 1 dann komplette Reihe = 1 (0 <= G <= 1)
            if G[jj][ll][mm] > 1:
                for b2 in range(n):
                    mplan[jj][ll][b2] = np.Inf

            # Wenn Erwartungswert der reaktiven Kosten den Schwellwert
            # Überschreitet, dann Zeile unbrauchbar machen
            if EWC[jj][ll][mm] + EWCEL[jj][ll][mm] >= drisk[jj]:
                for b3 in range(n):
                    mplan[jj][ll][b3] = np.Inf

    # Kosten durch Stillstand
    A[jj] = mplan[jj] * input.CEL[jj]
    # reduziertes Array mplan mplanred ohne Zeilen mit inf erstellen
    for i3 in range(m):
        for i4 in range(n):
            if mplan[jj][i3][i4] == np.Inf:

                q10[0][jj] = q10[0][jj] + 1
                break

    q10[0][jj] = m - q10[0][jj]

mplanred = np.ones((input.numcomp, int(np.max(q10)), n))
tmpRow = np.ones((int(np.max(q10)), n))
infFound = False
for jj2 in range(input.numcomp):
    c3 = 0
    for c1 in range(len(mplan[jj2])):
            if np.Inf not in (mplan[jj2][c1][:]):
                for c2 in range(n):
                    mplanred[jj2][c3][c2] = mplan[jj2][c1][c2]
                c3 = c3 + 1


            # Erwartungswert für Gesamtkosten durch spontanen Ausfall für
# alle Komponenten bestimmen
for ii in range(m):
    for nn in range(n):
        EWCGES[ii][nn] = np.sum(EWC[:, ii, nn])  # [:][ii][nn])


# Kosten der Wartungsfolgen bestimmen
# Kosten für jeweilige Wartungskombination bestimmen
intervallC.calcInterval(input.CP, input.CEL)

# Vorübergehende Berechnung zum Anzeigen
CGES = intervallC.calcCost(input.numcomp, q10, n, mplanred)
ergebnis = np.min(CGES)

rowIndex = []
#for entry in CGES:
row = (np.where(CGES == ergebnis))[0]
    #row.extend(entry)
        #row.extend(entry)

Y0 = intervallC.combmplan[row[3]]
Y0 = np.swapaxes(Y0, 0, 1)


# # # # # Anzeige-Abschnitt

currentBarHeight = np.zeros(n)

for p in range(input.numcomp):
    plt.bar(X0, Y0[p], 1500, bottom=currentBarHeight)
    currentBarHeight = currentBarHeight + Y0[p]

plt.xticks(X0, X0)
plt.yticks(range(input.numcomp + 1))
plt.legend(('Comp1', 'Comp2', 'Comp3', 'Comp4'))
plt.show()
