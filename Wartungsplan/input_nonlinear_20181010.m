global avhelp MTTR MTTRD TV DTRd DTRh CIS CR CP b CEL T n skalierung decisionfactor;


avhelp=1.0;
%Anzahl Komponenten
numcomp=7;

%Zeitschritte
n=150;

%Größe Population
PopSize=50;

% Faktor für Entscheidung, ob gewartet werden soll
decisionfactor=1.0;

skalierung=250; %Betriebsstunden pro Zeitschritt

%Komponentenzuordnung
%Komponente 1: Schwingungsdämpfer
%Komponente 2: Achslagerführung
%Komponente 3: Bremsbeläge
%Komponente 4: Gelenkwelle
%Komponente 5: Radreifen
%Komponente 6: Voith-Getriebe
%Komponente 7: Motor



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

%% Daten zum statischen Wartungsplan - Wartungsstufen


%Komponentenzuordnung
%Komponente 1: Schwingungsdämpfer
%Komponente 2: Achslagerführung
%Komponente 3: Bremsbeläge
%Komponente 4: Gelenkwelle
%Komponente 5: Radreifen
%Komponente 6: Voith-Getriebe
%Komponente 7: Dieselmotor


%W1
tW1h=20; %Häufigkeit
%W2
tW2h=500; %Häufigkeit
compW2=[1 0 1 0 1 0 1]; %Prüfung
dtW2d=3;
%Zeitstrahl Wartungsstufe W2
inspplanW2=zeros(1,n);
tpW2=tW2h/skalierung;
cW2=1; %counter
for i=1:length(inspplanW2)
    if cW2 == tpW2 
        inspplanW2(1,i)=1;
        cW2=1;
    else
        cW2=cW2+1;
    end
end

%W3
tW3h=1000; %Häufigkeit
compW3=[1 0 1 1 1 1 1]; %Prüfung
dtW3d=6;
%Zeitstrahl Wartungsstufe W3
inspplanW3=zeros(1,n);
tpW3=tW3h/skalierung;
cW3=1; %counter
for i=1:length(inspplanW3)
    if cW3 == tpW3 
        inspplanW3(1,i)=1;
        cW3=1;
    else
        cW3=cW3+1;
    end
end

%W4
tW4h=2000; %Häufigkeit
compW4=[1 1 1 1 1 1 1]; %Prüfung
dtW4d=12;
%Zeitstrahl Wartungsstufe W4
inspplanW4=zeros(1,n);
tpW4=tW4h/skalierung;
cW4=1; %counter
for i=1:length(inspplanW4)
    if cW4 == tpW4 
        inspplanW4(1,i)=1;
        cW4=1;
    else
        cW4=cW4+1;
    end
end


%W5
tW5h=12000; %Häufigkeit 
compW5=[1 1 1 1 1 1 1]; %Prüfung
dtW5d=25;
%Zeitstrahl Wartungsstufe W5
inspplanW5=zeros(1,n);
tpW5=tW5h/skalierung;
cW5=1; %counter
for i=1:length(inspplanW5)
    if cW5 == tpW5 
        inspplanW5(1,i)=1;
        cW5=1;
    else
        cW5=cW5+1;
    end
end


%W6
tW6h=32000; %Häufigkeit
compW6=[1 1 1 1 1 1 1]; %Prüfung
dtW6d=50;
%Zeitstrahl Wartungsstufe W6
inspplanW6=zeros(1,n);
tpW6=tW6h/skalierung;
cW6=1; %counter
for i=1:length(inspplanW6)
    if cW6 == tpW6 
        inspplanW6(1,i)=1;
        cW6=1;
    else
        cW6=cW6+1;
    end
end
%Anzahl Wartungsstufen
nws=5; %eigentlich 6, aber WS 1 nicht berücksichtigt

%Inspektionsplan
inspplan=zeros(6-1, n);
inspplan(1,:)=inspplanW2(1,:);
inspplan(2,:)=inspplanW3(1,:);
inspplan(3,:)=inspplanW4(1,:);
inspplan(4,:)=inspplanW5(1,:);
inspplan(5,:)=inspplanW6(1,:);
%bereinigen um doppelte Einträge
sizeinspplan=size(inspplan);
for i=1:sizeinspplan(1,2)
    for j=1:sizeinspplan(1,1)
        if inspplan(nws-j+1,i) == 1
            inspplan(:,i) = 0;
            inspplan(nws-j+1,i) = 1;
        end
    end
end
%Faktorisieren für Anzeige
inspplan_dspl=inspplan;
for i = 1:nws
    inspplan_dspl(i,:)=inspplan_dspl(i,:)*(i+1);
end
            
            
    
    

%% Print des Ergebnisses        
X0=linspace(0,(n-1)*skalierung,n);
figure
hold on
% for i=1:length(drawvec)
%     patch('XData', [drawvec(i)-1 drawvec(i)+ergAred(2,drawvec(i))-1 drawvec(i)+ergAred(2,drawvec(i))-1 drawvec(i)-1]*skalierung, 'YData', [0 0 numcomp numcomp], 'EdgeColor', 'none', 'FaceColor', 'red', 'FaceAlpha', 0.5);
% end
bar(X0,inspplan_dspl', 'stacked');
legend('WS2','WS3','WS4','WS5','WS6');

