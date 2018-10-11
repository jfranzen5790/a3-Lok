global PopSize skalierung numcomp asize cell_ci CGES n Cmax arbeitszeit transmat G EWC A TGESh allCosts costList gueltigCounter ungueltigCounter iterationsTillNewValidValue iterationsTillNewValidValue_ListOfallValues longestRuntime allInputValues multipleInputValuesCounter;

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
allCosts = {};
costList = zeros(0);
gueltigCounter = 0;
ungueltigCounter = 0;


iterationsTillNewValidValue = 0;
iterationsTillNewValidValue_ListOfallValues = 0;
longestRuntime = 0;
allInputValues = {};
multipleInputValuesCounter = 0;

addpath("../");
input_nonlinear_20181010;
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

% 'CrossoverFraction', 0.8, 'MigrationFraction', 0.2, , 'PlotFcn', 'gaplotbestf'
%options = optimoptions('ga', 'FunctionTolerance', 10e-4, 'PopulationSize', 150,  'MaxStallGenerations', 100, 'CrossoverFraction', 0.8, 'MigrationFraction', 0.2, 'InitialPopulation',IP);
%options = optimoptions('ga', 'FunctionTolerance', 10e-4, 'PopulationSize', 150,  'MaxStallGenerations', 100, 'CrossoverFraction', 0.8, 'MigrationFraction', 0.2);
% options = optimoptions('ga', 'PlotFcn', 'gaplotbestf', 'FunctionTolerance', 10e-8, 'MaxStallGenerations', 200);
%  options = optimoptions('gamultiobj','FunctionTolerance', 10e-4, 'MaxStallGenerations', 200);

%Letzte funktionierende Version
%options = optimoptions('ga', 'FunctionTolerance', 10e-4, 'PopulationSize', 150,  'MaxStallGenerations', 100, 'CrossoverFraction', 0.8, 'MigrationFraction', 0.2);

%options = optimoptions('ga', 'FunctionTolerance', 10e-8, 'PopulationSize', 150,  'MaxStallGenerations', 200, 'CrossoverFraction', 0.8, 'MigrationFraction', 0.2);
options = optimoptions('ga', 'FunctionTolerance', 1e-15, 'PopulationSize', 50,  'MaxStallGenerations', 400, 'CrossoverFraction', 0.8, 'MigrationFraction', 0.2 , 'PlotFcn', 'gaplotbestf',  'InitialPopulation',IP);
%options = optimoptions('ga', 'FunctionTolerance', 1e-15, 'PopulationSize', 50,  'MaxStallGenerations', 400, 'CrossoverFraction', 0.8, 'MigrationFraction', 0.2 , 'PlotFcn', 'gaplotbestf');

    %nlcon = @mycon;
    nlcon = @mycon;
    
    IntCon=(1:n);
    %[X,fval] = ga(@cost,n,[],[],[],[], lb, ub, nlcon, IntCon, options);
    %[X,fval] = ga(@cost_angepasst,n,[],[],[],[], lb, ub, nlcon);%, IntCon, options);
    [X,fval] = ga(@cost_angepasst,n,[],[],[],[], lb, ub, nlcon, IntCon, options);
% 
%% Simulated Annealing
%saoptions=optimoptions(@simulannealbnd,'PlotFcns',{@saplotbestx,...
%          @saplotbestf,@saplotx,@saplotf}, 'FunctionTolerance', 10e-4);
%%Randbedingungen für Simulated Annealing
%ubsa=ones(1,n)*asize(1,2);
%ubsa(1,1)=1;
%lbsa=ones(1,n);
% Festsetzen der Variablen, welche durch GA als keine Wartung festgelegt
% wurden
%sizeIP=size(IP);
%for i=1:sizeIP(1,2)
%    if IP(1,i)==asize(1,2)
%        lbsa(1,i)=asize(1,2);
%    end
%end
%[X2, fval2] = simulannealbnd(@cost,X,lbsa,ubsa,saoptions);
%     
    
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


figure
hold on
plot(costList, ".")
gueltigRate = gueltigCounter / (ungueltigCounter + gueltigCounter);
disp("Rate der gültigen Werte: " + gueltigRate);
%Ergebnis in 


%% Nichtlineare Nebenbedingungen
% 
% function [c, ceq]= nlcon(x)
%     if x(1) == 1
%         c(1) = x(16)-EWC(1,1);
%     end
%     ceq = [];
% end





%% Notizen

% Anzahl der Komponenten erhöhen

% Güte der gefundenen Lösung beurteilen

% Unterscheidung technische Verfügbarkeit/betriebliche Verfügbarkeit

% Ranken der Kombinationne




