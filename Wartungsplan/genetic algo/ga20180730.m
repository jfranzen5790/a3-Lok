global skalierung numcomp asize cell_ci CGES n Cmax arbeitszeit transmat G EWC A TGESh;

%Anzahl Zeitschritte
n=100;
skalierung=80; %Betriebsstunden pro Zeitschritt
%Umrechnung in Jahre für Geldwert
arbeitszeit=8;
TGESh=skalierung*n;
TGESd=TGESh/arbeitszeit;
TGESa=TGESd/220;
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
% maximal mögliche Kosten (zur Normierung)
Cmax=n*max(CGES);


%% Optimierung
lb=ones(1,n);
ub=ones(1,n);
ub(1,1:n)=asize(1,2);
ub(1,1)=1;
% ub(1,n)=1;

nonlcon=[];
% IntCon=1:n;

% options = gaoptimset('PlotFcns',@cost);

initPop;
options = optimoptions('ga', 'PlotFcn', 'gaplotbestf', 'FunctionTolerance', 10e-8, 'PopulationSize', 200,  'MaxStallGenerations', 1000, 'InitialPopulation',IP);
% options = optimoptions('ga', 'PlotFcn', 'gaplotbestf', 'FunctionTolerance', 10e-8, 'MaxStallGenerations', 200);
%  options = optimoptions('gamultiobj','FunctionTolerance', 10e-4, 'MaxStallGenerations', 200);

    IntCon=(1:n);
    [X,fval] = ga(@cost,n,[],[],[],[],lb, ub, [], IntCon, options);
    
%% Auswertung des Optimierungsergebnisses
% Ergebnis in Vektor und Wartungsfolge wandeln
% Xint=uint8(X(1,:));

[ergtrans, ergtransmat] = calctransmat(X);
ergC = calcCost(ergtrans);
% Stillstandzeiten pro Zeitschritt
[ergA, ergAred]=calcdt(ergtransmat);

% Lebensdauern bestimmen    
erglifetime=calclifetime(ergtransmat, ergAred);    

% Ausfallwahrscheinlichkeiten bestimmten
[ergG,ergEWC] = calcGmat(ergtransmat,erglifetime,1);
for i=1:n
    for j=1:numcomp
        ergEWC(j,i)=ergEWC(j,i)/CP(1,j);
    end
end
sumAred=sum(ergAred(2,:));
techav=calctechav(sumAred);
ergCGES=sum(ergC(1,:));
 
Xintsize=size(X);
erg=zeros(numcomp,Xintsize(1,2));
for j1=1:Xintsize(1,2)
    ergvec=cell_ci{1,X(j1)};
    for j2=1:numcomp
        erg(j2,j1)=ergvec(1,j2);
    end
end
% % Kosten und Ausfallwahrscheinlichkeiten für Zeitschritte bestimmen
% 
% erglifetime=zeros(numcomp,n);
% 
% for j3=1:n
%     for j4=1:numcomp
%         if j3 == 1
%             if erg(j4,j3)== 1
%                 erglifetime(j4,j3)=0;
%                 ergG(j4,j3)=0;
%                 ergA(j4,j3) = 0; %erste Wartung darf nicht mitzählen
%             elseif erg(j4,j3) == 0
%                 ergG(j4,j3)=calcG(j3*skalierung,j4);
%                 ergA(j4,j3) = 0;
%                 erglifetime(j4,j3)=j3*skalierung;
%             end
%         elseif (1 < j3) && (j3 < n)
%             if erg(j4,j3)== 1
%                 ergG(j4,j3)=0;
%                 ergwartung(1,j4)=j3;
%                 ergA(j4,j3) = MTTRD(1,j4);
%                 erglifetime(j4,j3)=0;
%             elseif erg(j4,j3) == 0
%                 ergG(j4,j3)=calcG((j3-ergwartung(1,j4))*skalierung,j4);
%                 erglifetime(j4,j3)=(j3-ergwartung(1,j4))*skalierung;
%             end
%         elseif j3 == n
%             if erg(j4,j3)== 1
%                 ergG(j4,j3)=0;
%                 ergwartung(1,j4)=j3;
%                 ergA(j4,j3) = 0;
%                 erglifetime(j4,j3)=0;
%             elseif erg(j4,j3) == 0
%                 ergG(j4,j3)=calcG((j3-ergwartung(1,j4))*skalierung,j4);
%                 erglifetime(j4,j3)=(j3-ergwartung(1,j4))*skalierung;
%             end
%         end
%             ergEWC(j4,j3)=(ergG(j4,j3)*(CEL(1,j4)+CR(1,j4)))/CP(j4);
%     end
% end
% 
% for j5=1:n
%     ergAred(1,j5)=max(ergA(:,j5));
%     ergAred(2,j5)=ergAred(1,j5)*8;
% end
% % avai=((TGESh-sum(ergAred(2,:)))/TGESh)*100;
% % disp(avai);
idraw=1;
for i=1:n
    if sum(transmat(:,i)) > 0
        drawvec(idraw)=i;
        idraw=idraw+1;
    end
end

        
%% Print des Ergebnisses        
X0=linspace(0,(n-1)*skalierung,n);
figure
subplot(2,1,1)
hold on
for i=1:length(drawvec)
    patch('XData', [drawvec(i)-1 drawvec(i)+ergAred(2,drawvec(i))-1 drawvec(i)+ergAred(2,drawvec(i))-1 drawvec(i)-1]*skalierung, 'YData', [0 0 numcomp numcomp], 'EdgeColor', 'none', 'FaceColor', 'red', 'FaceAlpha', 0.5);
end
bar(X0,ergtransmat', 'stacked');
strtext1=sprintf('LCC: %0.0f €', round(ergCGES));
strtext2=sprintf('Verfügbarkeit: %0.2f %%', techav*100);
text(TGESh-2000,-3,strtext1);
text(TGESh-2000,-2,strtext2);
% avaistr=num2str(avai);

subplot(2,1,2)
labels={'in Betrieb', 'außer Betrieb'};
pie([techav, 1-techav],labels)

%Ergebnis in 


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
    global A G EWC CP n C numcomp skalierung trans transmat TGESh Cmax;    
    %Umwandeln in Vektor
%     transmat=zeros(numcomp, n);
%     G=zeros(numcomp, n); % Ausfallwahrscheinlichkeiten
%     A=zeros(numcomp, n); % Ausfallzeiten
%     EWC =zeros(numcomp, n); %Erwartungswert der Audfallkosten
    avhelp=1; % Hilfsvariable, um Ergebnis ungültig zu machen
%     trans=cell(1,n);
%     x=round(x);

    %Wartungsplan aus Eingangsvariablen erstellen und Kosten pro Schritt
    %bestimmen
    [trans, transmat] = calctransmat(x);
    C = calcCost(trans);
    % Stillstandzeiten pro Zeitschritt
    [A, Ared]=calcdt(transmat);
    
% Lebensdauern bestimmen    
    lifetime=calclifetime(transmat, Ared);    
    
    % Ausfallwahrscheinlichkeiten bestimmten
    [G,EWC] = calcGmat(transmat,lifetime,1);
    
    
    for i=1:numcomp
        for j=1:n
            if EWC(i,j)>CP(1,i)
%                 C(1,j)=inf;
                avhelp=0;
            end
        end
    end
    
    gueltig=0; %Hilfsvariable zur Überprüfung, ob die Stillstandszeiten eingehalten werden
    for j5=1:n
        gueltig=0;
        dauer=0;
        if Ared(2,j5)>=1
            dauer=Ared(2,j5);
        end
        if dauer >= 2
            for j6=1:dauer -1
                if j5+j6 <= n
                    gueltig=gueltig+Ared(1,j5+j6);
                end
            end
        end
        if gueltig > 0
            avhelp=0; %Schleife abbrechen, wenn ungültige Reihe
            break;
        end
    end

    
    
%     ma=0;
%     if avhelp == 0
%         y(1)=inf;
%         y(2)=1;
%         y(3)=inf;
%     elseif avhelp >= 1
%         y(1) = sum(C(1,:));
%         y(2) = 1-((TGESh-(sum(Ared(1,:))*8))/TGESh);
%         % Anzahl der Stillstände minieren
%         for i=1:n
%                 ma = ma + sum(transmat(:,i));
%         end
%         y(3) = ma;
%         if y(2) < 0
%             disp("Negative Zahl");
%         end
%     end
    if avhelp == 0
        y=1.0;
    elseif avhelp == 1
        sumC=sum(C(1,:));
        sumAred=sum(Ared(2,:));
        %Ausfallzeit bestimmen
        techav=calctechav(sumAred);
        
        y=0.7*(sumC/Cmax)+0.3*(1-techav);
    end
end
%% Fkt. für Best. der Ausfallwahrscheinlichkeit zum Zeitpunkt x
function aw = calcG(x,i)
    global T;
    aw = 1-exp(-x/T(i));
end

%techn. Verfügbarkeit bestimmen
function tav=calctechav(sAred)
    global skalierung TGESh;
    dtcycle=sAred*skalierung;
    tav=((TGESh-dtcycle)/TGESh);
end

% Zeitstrahl berechnen
function [transcell,transmat2] = calctransmat(x)
    global n numcomp cell_ci;
    transcell=cell(1,n);
    transmat2=zeros(numcomp,n);
    for i=1:n
        transcell{1,i}=cell_ci{1,x(i)};
        transvec=cell_ci{1,x(i)};
        for j=1:numcomp
            transmat2(j,i)=transvec(1,j);        
        end   
    end
end

function C = calcCost(transcell)
    global n cell_ci CGES asize skalierung arbeitszeit;
    C=zeros(1,n);
    for i=1:n
        for j=1:asize(1,2)
            if transcell{1,i} == cell_ci{1,j}
                C(1,i)= CGES(1,j)*1.02^(i*(skalierung/(200*arbeitszeit))); %Zinssatz berücksichtigen
                break
            elseif transcell{1,i} ~= cell_ci{1,j}
                C(1,i) = 0;
            end
        end
    end
end

function [Amat,Aredmat]= calcdt(transmat2)
    global n numcomp skalierung MTTRD
    Amat=zeros(numcomp,n);
    Aredmat=zeros(2,n);
    for i=1:numcomp
        for j=1:n
            % erster Durchlauf "System neu" darf nicht berücksichtigt
                % werden
                if j == 1
                    Amat(i,j)=0;
                elseif (1 < j) && (j < n)
                    if transmat2(i,j) == 0
                        Amat(i,j)=0;
                    elseif transmat2(i,j) == 1
                        Amat(i,j)=MTTRD(1,i);
                    end
                elseif j == n
                    Amat(i,j)=0;
                end
        end
    end
    
    % Anzahl Zeitschritte Wartung pro Zeitschritt
    for i=1:n
        Aredmat(1,i)=max(Amat(:,i));
        %Zeitschritte umrechnen
        Aredmat(2,i)=ceil(Aredmat(1,i)*8/skalierung);
    end
end

function lt = calclifetime(transmat2, Aredmat)
    global n numcomp skalierung
    lt=zeros(numcomp,n);
    for i = 1:numcomp
        wartung=ones(1,numcomp);
        jlt=1;
        skipi=0;
        while jlt <= n
          
            if jlt==1
                if transmat2(i,jlt) == 0
                    lt(i,jlt)=jlt*skalierung;
                elseif transmat2(i,jlt) == 1
                    lt(i,jlt)=0;
                    wartung(1,i)=jlt;
                end
            elseif jlt > 1
                if sum(transmat2(:,jlt)) == 0
                    if transmat2(i,jlt-skipi) == 0
                        lt(i,jlt)=lt(i,wartung(1,i))+(jlt-wartung(1,i))*skalierung;
                    elseif transmat2(i,jlt-skipi) == 1
                        lt(i,jlt)=(jlt-wartung(1,i))*skalierung;
                    end
                    skipi=0;
        
                elseif sum(transmat2(:,jlt)) > 0
                    skipi=Aredmat(2,jlt);
                    wartung(1,i)=jlt-1+skipi;
                    for k=1:Aredmat(2,jlt)
                        if jlt-1+k<n
                            if transmat2(i,jlt) == 1
                                lt(i,jlt-1+k)=0;
                            elseif transmat2(i,jlt) == 0
                                lt(i,jlt-1+k)=lt(i,jlt-1);
                            end
                        end
                    end
                end
            end
            if skipi == 0
                jlt=jlt+1;
            elseif skipi > 0
                jlt=jlt+skipi;
            end
        end
        jlt=1;
    end
end

function [G, EWC] = calcGmat(transmat2, lt,factC)
    global n numcomp CEL CR 
    G=zeros(numcomp,n);
    EWC=zeros(numcomp,n);
    for i=1:numcomp
%         wartung=ones(1,numcomp);
        for j=1:n
            if j==1
                if transmat2(i,j) == 0
                    G(i,j)=calcG(lt(i,j),i);
                elseif transmat2(i,j) == 1
                    G(i,j) = 0;
                end
            elseif j >= 2
                if transmat2(i,j) == 0
                    G(i,j)=calcG(lt(i,j),i);
                elseif transmat2(i,j) == 1
                    G(i,j) = 0;
                end
            end
            if G(i,j) == 0
                EWC(i,j) = 0;
            elseif G(i,j) > 0
                EWC(i,j)=G(i,j)*(CEL(1,i)+CR(1,i))/factC;
            end
        end
    end
end

%% Notizen

% Labels mit Informationen zur optimalen Lösung

% Anzahl der Komponenten erhöhen

% Güte der gefundenen Lösung beurteilen

% Stillstände als Zeitfenster markieren/einblenden

% Unterscheidung technische Verfügbarkeit/betriebliche Verfügbarkeit

% Ranken der Kombinationne