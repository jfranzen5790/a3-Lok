global MTTR MTTRD TV DTRd DTRh CIS CR b T;



%Gesamte Betriebszeit
% In Jahren
TGESA=5;
% In Tagen (220 Arbeitstage pro Jahr)
TGESD=TGESA*220;
% In Stunden bei 8 Arbeitsstunden pro Tag
arbeitszeit=24;
TGES=TGESD*arbeitszeit;

% Schrittweite
sw=100;

% Mittlere Reparaturdauer der jeweiligen Komponente, Stunden
MTTR=[60 10];


%Stillstand
%präventiv in Tagen
MTTRD=[56 10];

% Verögerung im reaktiven Fall in Tagen
TV=0.3*MTTRD;

% Gesamte Ausfallzeit (reaktiv) in Tagen 

DTRd=[MTTRD(1)+TV(1) MTTRD(2)+TV(2)];
%Gesamte Ausfallzeit (=Stillstandszeit) in Stunden
DTRh=DTRd*24;


% Kosten für Instandsetzung/Ersetzung Komponente
CIS=[3000 4000];
%Werkersatz pro Stunde
M=100;

%Kosten für die Wartung = Reparaturdauer * Werkersatz
CC=[MTTR(1)*M MTTR(2)*M];
%Pauschalkosten für Werkstattaufenthalt
CWS=500;

%Kostensatz Ausfall des Gesamtsystems pro Tag
% Miete pro Monat für Ersatzlok (G2000)
CELM=25000;
% Pro Tag
CELd=CELM/30;
% Pro Stunde
CELh=CELd/24;
 
% %Kosten für Ersatzleistung
% CA=[DTRh(1)*CELh DTRh(2)*CELh];

% Fixe Kosten bei präventiver Wartung
for i=1:length(CIS)
    CP(i)= CIS(i) + CC(i) + (CELd * MTTRD(i));
end
CR=CP*1.3;
CEL=[CELd*DTRd(1) CELd*DTRd(2)];


%Geforderte Verfügbarkeit
AVAILIBILITY=1;

%Weibull-Parameter
b=[1 1];
T=[5000 7000];