global MTTR MTTRD TV DTRd DTRh CIS CR CP b CEL T n skalierung decisionfactor;

%Anzahl Komponenten
numcomp=4;

%Zeitschritte
n=50;

%Größe Population
PopSize=50;

arbeitszeit=8;
% In Betriebsstunden
TGESh=6000;
%Später rausnehmen!
TGESd=TGESh/arbeitszeit;

%
decisionfactor=1.0;

skalierung=80; %Betriebsstunden pro Zeitschritt


% Mittlere Reparaturdauer der jeweiligen Komponente, Stunden
MTTR=[40 10 30 40 4 15 25]; 


%Stillstand
%präventiv in Tagen
MTTRD=[56 10 28 30 30 10 15];

% Verögerung im reaktiven Fall in Tagen
TV=0.3*MTTRD;

% Gesamte Ausfallzeit (reaktiv) in Tagen 

DTRd=[MTTRD(1)+TV(1) MTTRD(2)+TV(2) MTTRD(3)+TV(3) MTTRD(4)+TV(4) MTTRD(5)+TV(5) MTTRD(6)+TV(6) MTTRD(7)+TV(7)];
%Gesamte Ausfallzeit (=Stillstandszeit) in Stunden
DTRh=DTRd*24;


% Kosten für Instandsetzung/Ersetzung Komponente
CIS=[15000 4000 1000 20000 3000 25000 15000];
%Werkersatz pro Stunde
M=100;

%Kosten für die Wartung = Reparaturdauer * Werkersatz
CC=[MTTR(1) MTTR(2) MTTR(3) MTTR(4) MTTR(5) MTTR(6) MTTR(7)]*M;
%Pauschalkosten für Werkstattaufenthalt
CWS=500;

%Kostensatz Ausfall des Gesamtsystems pro Tag
% Miete pro Monat für Ersatzlok (G2000)
CELM=25000;
% Pro Tag
CELd=CELM/30;
% Pro Stunde
CELh=CELd/24;

% Fixe Kosten bei präventiver Wartung
for i=1:numcomp
    CP(i)= CIS(i) + CC(i) + (CELd * MTTRD(i));
end
CR=CP*1.3;
CEL=[DTRd(1) DTRd(2) DTRd(3) DTRd(4) DTRd(5) DTRd(6) DTRd(7)]*CELd;


%Geforderte Verfügbarkeit
AVAILIBILITY=1;

%Weibull-Parameter
b=[1 1 1 1 1 1 1];
T=[10000 5000 7500 2500 6000 8000 3000];
%T=[10000 5000 7500];