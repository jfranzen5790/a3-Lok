import numpy as np

# Anzahl Komponenten
numcomp = 4

# Gesamte Betriebszeit
# In Jahren
TGESA = 5
# In Tagen(220 Arbeitstage pro Jahr)
TGESD = TGESA * 220
# In Stunden bei 8 Arbeitsstunden pro Tag
arbeitszeit = 24
TGES = TGESD * arbeitszeit

# Schrittweite
# sw = 100;

# Mittlere Reparaturdauer der jeweiligen Komponente, Stunden
MTTR = [60, 10, 30, 40]

# Stillstand
# präventiv in Tagen
MTTRD = [56, 10, 28, 30]

# Verögerung im reaktiven Fall in Tagen
TV = 0.3 * np.asarray(MTTRD)

# Gesamte Ausfallzeit(reaktiv) in Tagen

DTRd = [MTTRD[0] + TV[0], MTTRD[1] + TV[1], MTTRD[2] + TV[2], MTTRD[3] + TV[3]]
# Gesamte Ausfallzeit( = Stillstandszeit) in Stunden
DTRh = DTRd * 24

# Kosten für Instandsetzung / Ersetzung Komponente
CIS = [100, 4000, 100, 2000]
# Werkersatz pro Stunde
M = 100

# Kosten für die Wartung = Reparaturdauer * Werkersatz
CC = [MTTR[0] * M, MTTR[1] * M, MTTR[2] * M, MTTR[3] * M]
# Pauschalkosten für Werkstattaufenthalt
CWS = 500

# Kostensatz Ausfall des Gesamtsystems pro Tag
# Miete pro Monat für Ersatzlok(G2000)
CELM = 25000
# Pro Tag
CELd = CELM / 30
# Pro Stunde
CELh = CELd / 24

# Kosten für Ersatzleistung
# CA = [DTRh(1) * CELh DTRh(2) * CELh];

CP = np.zeros(numcomp)
# Fixe Kosten bei präventiver Wartung
for i in range(numcomp):
    CP[i] = CIS[i] + CC[i] + (CELd * MTTRD[i])
CR = CP * 1.3
CEL = [CELd * DTRd[0], CELd * DTRd[1], CELd * DTRd[2], CELd * DTRd[3]]


# Geforderte Verfügbarkeit
AVAILIBILITY = 1

# Weibull - Parameter
b = [1, 1, 1, 1]
T = [10000, 5000, 10000, 2500]
