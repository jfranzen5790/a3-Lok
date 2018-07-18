global skalierung numcomp asize cell_ci CGES n C transmat G EWC A TGESh;

%Anzahl Zeitschritte
n=21;
skalierung=80; %Betriebsstunden pro Zeitschritt
%Umrechnung in Jahre für Geldwert
arbeitszeit=8;
TGESh=skalierung*n;
TGESd=TGESh/arbeitszeit;
TGESa=TGESd/220;
C=zeros(1,n);
transmat=zeros(numcomp,n);
G=zeros(numcomp,n);
EWC=zeros(numcomp,n);
A=zeros(numcomp,n);
ergC=zeros(numcomp,n);
ergEWC=zeros(numcomp,n);
ergG=zeros(numcomp,n);
ergA=zeros(numcomp,n);
ergwartung=ones(1,numcomp);
ergAred=zeros(2,n);


input_nonlinear;
generatecombinations;

% Kosten für obig bestimmte Wartungskombinationen pro Zeitschritt
combCP=zeros(asize(1,1),asize(1,2));
combCEL=zeros(asize(1,1),asize(1,2));
combMTTRD=zeros(asize(1,1),asize(1,2));
for a1=1:asize(1,1)
    for a2=1:asize(1,2)
        if combintervall(a1,a2) == 1
            combCP(a1,a2)=CP(a1);
            combCEL(a1,a2)=CEL(a1);
            combMTTRD(a1,a2)=MTTRD(a1);
        end
    end
end
% Arrays mit präventiven Wartungskosten sowie Kosten für die Ersatzleistung
combCELsize=size(combCEL);
combMTTRDsize=size(combMTTRD);
combCELmax=zeros(1,combCELsize(1,2));
combMTTRDmax=zeros(1,combMTTRDsize(1,2));
for c4=1:length(combCEL)
    if sum(combCEL(:,c4))==0
        combCELmax(1,c4)=0;
        combMTTRDmax(1,c4)=0;
    elseif sum(combCEL(:,c4))==1
        for c5=1:combCELsize(1,1)
            if combCEL(1,c5)>0
                combCELmax(1,c4)=combCEL(1,c5);
                break
            elseif combCEL(2,c5)>0
                combCELmax(1,c4)=combCEL(2,c5);
                break
            elseif combCEL(3,c5)>0
                combCELmax(1,c4)=combCEL(3,c5);
                break
            end
        end
    elseif sum(combCEL(:,c4))>1
        combCELmax(1,c4)=max(combCEL(:,c4));
    end
    combMTTRDmax(1,c4)=max(combMTTRD(:,c4));
end
CGES=zeros(1,length(cell_ci));
for i2=1:length(cell_ci)
    CGES(1,i2)=combCELmax(1,i2)+sum(combCP(:,i2));
end

%% Optimierung
lb=ones(1,n);
ub=ones(1,n);
ub(1,1:n)=asize(1,2);
ub(1,1)=1;
lb(1,n)=1;
ub(1,n)=1;

nonlcon=[];
% IntCon=1:n;

% options = gaoptimset('PlotFcns',@cost);

IP=zeros(7,n);
IP(1,:)=[1 16 16 16 16 8 16 16 16 16 5 16 16 16 16 16 16 7 16 16 1];
IP(2,:)=[1 16 4 1 3 16 16 12 16 16 8 16 16 16 16 16 16 16 16 16 1];
IP(3,:)=[1 16 16 7 16 16 16 16 11 16 8 16 5 16 16 16 16 16 16 16 1];
IP(4,:)=[1 16 16 16 16 8 16 16 10 16 16 16 16 16 6 16 14 16 16 16 1];
IP(5,:)=[1 16 16 16 4 16 16 16 16 16 8 16 16 4 16 8 16 16 16 16 1];
IP(6,:)=[1 16 5 12 16 16 14 16 12 16 16 16 16 13 16 1 16 16 16 16 1];
IP(7,:)=[1 16 16 16 8 16 15 16 16 16 8 16 12 16 16 8 16 16 16 16 1];

options = optimoptions('gamultiobj','FunctionTolerance', 10e-8, 'MaxStallGenerations', 200, 'InitialPopulation',IP);

[X,fval] = gamultiobj(@cost,n,[],[],[],[],lb, ub, nonlcon, options);

%% Auswertung des Optimierungsergebnisses
% Ergebnis in Vektor und Wartungsfolge wandeln
Xint=uint8(X(2,:));
Xintsize=size(Xint);
erg=zeros(numcomp,Xintsize(1,2));
for j1=1:Xintsize(1,2)
    ergvec=cell_ci{1,Xint(j1)};
    for j2=1:numcomp
        erg(j2,j1)=ergvec(1,j2);
    end
end
% Kosten und Ausfallwahrscheinlichkeiten für Zeitschritte bestimmen

erglifetime=zeros(numcomp,n);

for j3=1:n
    for j4=1:numcomp
        if j3 == 1
            if erg(j4,j3)== 1
                erglifetime(j4,j3)=0;
                ergG(j4,j3)=0;
                ergA(j4,j3) = 0; %erste Wartung darf nicht mitzählen
            elseif erg(j4,j3) == 0
                ergG(j4,j3)=calcG(j3*skalierung,j4);
                ergA(j4,j3) = 0;
                erglifetime(j4,j3)=j3*skalierung;
            end
        elseif (1 < j3) && (j3 < n)
            if erg(j4,j3)== 1
                ergG(j4,j3)=0;
                ergwartung(1,j4)=j3;
                ergA(j4,j3) = MTTRD(1,j4);
                erglifetime(j4,j3)=0;
            elseif erg(j4,j3) == 0
                ergG(j4,j3)=calcG((j3-ergwartung(1,j4))*skalierung,j4);
                erglifetime(j4,j3)=(j3-ergwartung(1,j4))*skalierung;
            end
        elseif j3 == n
            if erg(j4,j3)== 1
                ergG(j4,j3)=0;
                ergwartung(1,j4)=j3;
                ergA(j4,j3) = 0;
                erglifetime(j4,j3)=0;
            elseif erg(j4,j3) == 0
                ergG(j4,j3)=calcG((j3-ergwartung(1,j4))*skalierung,j4);
                erglifetime(j4,j3)=(j3-ergwartung(1,j4))*skalierung;
            end
        end
            ergEWC(j4,j3)=(ergG(j4,j3)*(CEL(1,j4)+CR(1,j4)))/CP(j4);
    end
end

for j5=1:n
    ergAred(1,j5)=max(ergA(:,j5));
    ergAred(2,j5)=ergAred(1,j5)*8;
end
avai=((TGESh-sum(ergAred(2,:)))/TGESh)*100;
disp(avai);

        
%% Print des Ergebnisses        
X0=linspace(0,(n-1)*skalierung,n);
figure
hold on
patch('XData', [0 0 500 500], 'YData', [0 4 4 0], 'EdgeColor', 'none', 'FaceColor', 'red', 'FaceAlpha', 0.5);
bar(X0,erg', 'stacked');
avaistr=num2str(avai);
fprintf('Verfügbarkeit: %s Prozent', avaistr);

% function f=objfcn(x)
%     f = (x(1) + x(2))^2 - x(3);
% end

% %% Outputfkt.
% 
% function stop = outfun()
%     stop = false;
% end


%% Nichtlineare Nebenbedingungen
% 
% function [c, ceq]= nlcon(x)
%     if x(1) == 1
%         c(1) = x(16)-EWC(1,1);
%     end
%     ceq = [];
% end

%% Zielfkt. Optimierung
function y=cost(x)
    global A G EWC CEL CR asize cell_ci CP CGES n C numcomp transmat skalierung MTTRD TGESh;    
    %Umwandeln in Vektor
    transmat=zeros(numcomp, n);
    G=zeros(numcomp, n); % Ausfallwahrscheinlichkeiten
    A=zeros(numcomp, n); % Ausfallzeiten
    EWC=zeros(numcomp, n); %Erwartungswert der Audfallkosten
    avhelp=1; % Hilfsvariable, um Ergebnis ungültig zu machen
    lifetime = zeros(numcomp, n); % Lebensdauer beim jeweiligen Zeitschritt
    Ared=zeros(2,n); % Red. Stillstandszeiten
    trans=cell(1,n);
    x=round(x);
    for i1=1:n
        trans{1,i1}=cell_ci{1,x(i1)};
        transvec=cell_ci{1,x(i1)};
        for i7=1:numcomp
            transmat(i7,i1)=transvec(1,i7);        
        end
        
    end
    %Kosten für jeden Zeitpunkt bestimmen
    for i2=1:n
        for i3=1:asize(1,2)
            if trans{1,i2} == cell_ci{1,i3}
                C(1,i2)= CGES(1,i3)*(1.02*floor(i2/8));
                break
            elseif trans{1,i2} ~= cell_ci{1,i3}
                C(1,i2) = 0;
            end
        end
    end
    
    %Ausfallwahrscheinlichkeiten bestimmten
    for i5=1:numcomp
        wartung=ones(1,numcomp);
        for i4=1:n
            if i4==1
                if transmat(i5,i4) == 0
                    lifetime(i5,i4)=i4*skalierung;
                    G(i5,i4)=calcG(lifetime(i5,i4),i5);
                elseif transmat(i5,i4) == 1
                    G(i5,i4) = 0;
                    lifetime(i5,i4)=0;
                end
            elseif i4 >= 2
                if transmat(i5,i4) == 0
                    lifetime(i5,i4)=(i4-wartung(1,i5))*skalierung;
                    G(i5,i4)=calcG(lifetime(i5,i4),i5);
%                     G(i5,i4)=calcG((i4-wartung(1,i5))*skalierung,i5);
                elseif transmat(i5,i4) == 1
                    lifetime(i5,i4)=0;
                    G(i5,i4) = 0;
                    wartung(1,i5)=i4;
                end
            end
            if G(i5,i4) == 0
                EWC(i5,i4) = 0;
            elseif G(i5,i4) > 0
                EWC(i5,i4)=G(i5,i4)*(CEL(1,i5)+CR(1,i5));
            end
%             %Prüfen, ob EW der Kosten > als Präventive Kosten sind
%             if EWC(i5,i4)>=0
%                 C(1,i4)=10000000;
%                 disp('Hallo');
%             end
        end
    end
    for i6=1:numcomp
        for i8=1:n
            EWC(i6,i8)=G(i6,i8)*(CEL(1,i6)+CR(1,i6));
                % erster Durchlauf "System neu" darf nicht berücksichtigt
                % werden
                if i8 == 1
                    A(i6,i8)=0;
                elseif (1 < i8) && (i8 < n)
                    if transmat(i6,i8) == 0
                        A(i6,i8)=0;
                    elseif transmat(i6,i8) == 1
                        A(i6,i8)=MTTRD(1,i6);
                    end
                elseif i8 == n
                    A(i6,i8)=0;
                end
        end
    end
    for j1=1:numcomp
        for j2=1:n
            if EWC(j1,j2)>=CP(1,j1)
                C(1,j2)=inf;
                avhelp=0;
%                 x(1,n+1)=-1;
            end
        end
    end
    
    for j5=1:n
        Ared(1,j5)=max(A(:,j5));
        %Zeitschritte umrechnen
        Ared(2,j5)=ceil(Ared(1,j5)*8/skalierung);
    end
    if avhelp == 0
        y(2)=inf;
    elseif avhelp >= 0
        y(2)=1-((TGESh-sum(Ared(1,:))*8)/TGESh);
    end
% %     nahe beeinanderliegende Lösung rausfiltern
%     for j6=1:n-max(Ared(2,:))
%         if sum(transmat(:,j6+1:j6+max(Ared(2,:)))) > 0
%             C(1,j6)=inf;
%                 
%         end
%     end
    ma=0;
    for j5=1:numcomp    
        ma = ma + sum(transmat(j5,:));
    end
    y(3) = ma;
    y(1) = sum(C(1,:));
end


%% Fkt. für Best. der Ausfallwahrscheinlichkeit zum Zeitpunkt x
function aw = calcG(x,i)
    global T;
    aw = 1-exp(-x/T(i));
end


%% Notizen

% Labels mit Informationen zur optimalen Lösung

% Anzahl der Komponenten erhöhen

% Güte der gefundenen Lösung beurteilen

% Ausfallzeit des Fahrzeugs wird bisher nicht auf alle Komponenten
% übertragen

% Stillstände als Zeitfenster markieren/einblenden

% Unterscheidung technische Verfügbarkeit/betriebliche Verfügbarkeit

