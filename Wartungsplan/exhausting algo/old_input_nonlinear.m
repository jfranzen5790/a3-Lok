global MTTR MTTRD TV DTRd DTRh CIS CR b T;

%Anzahl Komponenten
numcomp=4;

%Gesamte Betriebszeit
% In Stunden bei 8 Arbeitsstunden pro Tag
arbeitszeit=8;
% In Betriebsstunden
TGESh=6000;
% % In Tagen (220 Arbeitstage pro Jahr)
TGESd=TGESh/arbeitszeit;

% TGES=TGESD*arbeitszeit;

% % Schrittweite
% sw=100;

% Mittlere Reparaturdauer der jeweiligen Komponente, Stunden
MTTR=[60 10 30 40];


%Stillstand
%pr�ventiv in Tagen
MTTRD=[56 10 28 30];

% Ver�gerung im reaktiven Fall in Tagen
TV=0.3*MTTRD;

% Gesamte Ausfallzeit (reaktiv) in Tagen 

DTRd=[MTTRD(1)+TV(1) MTTRD(2)+TV(2) MTTRD(3)+TV(3) MTTRD(4)+TV(4)];
%Gesamte Ausfallzeit (=Stillstandszeit) in Stunden
DTRh=DTRd*24;


% Kosten f�r Instandsetzung/Ersetzung Komponente
CIS=[100 4000 100 2000];
%Werkersatz pro Stunde
M=100;

%Kosten f�r die Wartung = Reparaturdauer * Werkersatz
CC=[MTTR(1)*M MTTR(2)*M MTTR(3)*M MTTR(4)*M];
%Pauschalkosten f�r Werkstattaufenthalt
CWS=500;

%Kostensatz Ausfall des Gesamtsystems pro Tag
% Miete pro Monat f�r Ersatzlok (G2000)
CELM=25000;
% Pro Tag
CELd=CELM/30;
% Pro Stunde
CELh=CELd/24;
 
% %Kosten f�r Ersatzleistung
% CA=[DTRh(1)*CELh DTRh(2)*CELh];

% Fixe Kosten bei pr�ventiver Wartung
for i=1:numcomp
    CP(i)= CIS(i) + CC(i) + (CELd * MTTRD(i));
end
CR=CP*1.3;
CEL=[CELd*DTRd(1) CELd*DTRd(2) CELd*DTRd(3) CELd*DTRd(4)];


%Geforderte Verf�gbarkeit
AVAILIBILITY=1;

%Weibull-Parameter
b=[1 1 1 1];
T=[10000 5000 10000 2500];