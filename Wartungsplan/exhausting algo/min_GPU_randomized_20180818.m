global PopSize skalierung numcomp asize cell_ci CGES n Cmax arbeitszeit transmat G EWC A TGESh allCosts costList gueltigCounter ungueltigCounter bestCost bestComb bestInput iterationsTillNewValidValue iterationsTillNewValidValue_ListOfallValues longestRuntime   x A G EWC CP n C numcomp decisionfactor skalierung trans transmat TGESh Cmax;
%%Testabschnitt zum Bestimmen der CSV-Werte

%for csvNCounter = 6:5
%    for skalierungCounter = 1:5
        addpath("../");
        input_nonlinear;
        %Anzahl Zeitschritte
%        n=csvNCounter;
        %n=2;
%        skalierung= skalierungCounter * 400; %Betriebsstunden pro Zeitschritt
        %skalierung= 2000; %Betriebsstunden pro Zeitschritt
        X0=linspace(0,(n-1)*skalierung,n);


        gueltigCounter = 0;
        ungueltigCounter = 0;
        iterationsTillNewValidValue = 0;
        iterationsTillNewValidValue_ListOfallValues = 0;
        longestRuntime = 0;
        gueltigRate = 0;
        
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


%% Neue reduzierte Berechnung von mplanred für größere n-Werte und ncomp Werte

randomCombination = zeros(1, n);
%amountMaxExpectedValues = 1000000;
amountMaxExpectedValues = 10000;
amountAllPossibleValues = (2^n)^numcomp;

if amountAllPossibleValues > amountMaxExpectedValues
    amountValues = amountMaxExpectedValues;
else
    amountValues = amountAllPossibleValues;
end

randomMaintenanceCombinations = zeros(amountValues, n);
maintenanceCombinationsCounter = 1;

while maintenanceCombinationsCounter <= amountValues
    for t = 1:n
       randomCombination(t) = randperm(2^numcomp,1);
    end 
   if   ~isempty(randomMaintenanceCombinations)
        if  ~ismember(randomCombination, randomMaintenanceCombinations, "rows")
            randomMaintenanceCombinations(maintenanceCombinationsCounter, :) = randomCombination;
            maintenanceCombinationsCounter = maintenanceCombinationsCounter + 1;
        else
           disp("Else"); 
        end
   end
end
        


%% calccost in abgespeckter Form

currentInput = zeros(0);
bestInput = zeros(0);
bestCost = 1;
currentCost = 1;

bestComb = zeros(1, n);
bestComb(:) = 2^numcomp;
currentComb = zeros(1, n);
currentComb(:) = 2^numcomp;

 TestArray_GPU = gpuArray(randomMaintenanceCombinations);
 testtest = bestCost;
 
 %TestArray_GPU = gpuArray(1:1000);
 
[bestCost] = arrayfun(@calculateCurrentSetsCost, TestArray_GPU, bestCost, x, A, G, EWC, CP, n, C, numcomp, decisionfactor, skalierung, trans, transmat, TGESh, Cmax);
z = gather(z);
disp("Test: ");
disp(test);
disp("z: ");
disp(z);

% for value = 1:amountValues
%     %combmplan{i4,xx}=vector;
%     currentInput = randomMaintenanceCombinations(value, :);
%     %currentInputIndexVersion;
%     currentCost = (cost(currentInput));
%     if currentCost < bestCost
%         bestCost = currentCost;
%         bestComb = currentInput;
%         bestInput = currentInput;
%     end
% end

% Minimale Kosten identifizieren

ergebnis=min(CGES);
row=find(CGES==ergebnis);
%%%%calcavailibility;
%Umrechnen
% k_1=floor(k/m);
% k_2(1)=k(1)-k_1(1)*m;
% Stillstand durch präventive Wartung
% A=zeros(length(k),1);
% for oo=1:length(k)
%     pp=k(oo);
%     A(oo)=sum(mplan1(pp,:))*MTTRD(1);
% end
%Y01perm=cell2mat(bestInput);
%Y0=zeros(length(Y01perm)/numcomp,numcomp);
%for i5=1:n
%    for i7=1:numcomp
%        Y0(i5,i7)=Y01perm(1,(i5-1)*numcomp+i7);
%    end
%end
% Y01=shiftdim(Y01perm,1);
% Y02perm=permute(mplan(row(1),:,2),[3 1 2]);
% Y02=shiftdim(Y02perm,1);
% Y0=zeros(length(Y01),2);
% for y0i=1:length(Y01)
%     Y0(y0i,1)= Y01(y0i);
%     Y0(y0i,2)= Y02(y0i);
% end
% G01perm=permute(G(column(1),:,1),[3 1 2]);
% G01=shiftdim(G01perm,1);
% G02perm=permute(G(row(1),:,2),[3 1 2]);
% G02=shiftdim(G02perm,1);
% EWC0=shiftdim(EWC(k(1),:,1),1);
%figure
% subplot(3,1,1);
% plot(X0,G01);
% if sum(G01) == 0
%     limx=1;
% elseif sum(G01) > 0
%     limx=max(G01);
% end
% if sum(G02) == 0
%     limx2=1;
% elseif sum(G02) > 0
%     limx2=max(G02);
% end
% 
% 
% axis([0 n*skalierung 0 limx]);
% subplot(3,1,2);
% % plot(X0,G02);
% axis([0 n*skalierung 0 limx2]);
%subplot(3,1,3);
%bar(X0,Y0, 'stacked');
% axis([0 n*skalierung 0 1.5]);


%% Auswertung des Optimierungsergebnisses
[ergtrans, ergtransmat] = calctransmat(bestComb);
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
 
Xintsize=size(bestComb);
erg=zeros(numcomp,Xintsize(1,2));
for j1=1:Xintsize(1,2)
    ergvec=cell_ci{1,bestComb(j1)};
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


%% CSV und gültige Raten

gueltigRate = gueltigCounter / (ungueltigCounter + gueltigCounter);
disp("Rate der gültigen Werte: " + gueltigRate);

figure
hold on
plot(costList, ".");
str = {strcat('N = ', num2str(n)),strcat('nComp = ', num2str(numcomp)), {'Rate gültiger Werte', strcat('= ', num2str(gueltigRate))}};
plot_axis5percent = length(costList)*0.05;
text([plot_axis5percent plot_axis5percent plot_axis5percent],[0.9 0.8 0.7],str)


csvMatrix = csvread('Rate_gültigeWerte_v2.dat');

csvSeptuple = [numcomp, n, skalierung, gueltigRate, longestRuntime, amountValues, amountAllPossibleValues];
csvMatrix(end + 1, :) = csvSeptuple;
csvwrite('Rate_gültigeWerte_v2.dat', csvMatrix);
type Rate_gültigeWerte_v2.dat;

csvAuswertungAktiviert = 0;

if csvAuswertungAktiviert == 1
    %% CSV Auswerungsbereich
    csvMatrix = csvread('Rate_gültigeWerte_v2.dat');
    %Beginn bei Index 2, da erste Zeile aus Nullen besteht!
    %barY = csvMatrix;
    for row = 1:size(csvMatrix, 1)
       
        % Alle Komponenten aus n=2
       if csvMatrix(row, 2) == 2
          %barY_n2 = [barY; csvMatrix(row, :)]; 
       end
    end
    
    %Vorerst händisch!
    %Nach N sortiert
    barY_n2(1) = csvMatrix(2, 4);
    barY_n2(2) = csvMatrix(11, 4);
    barY_n2(3) = csvMatrix(20, 4);
    barY_n2(4) = csvMatrix(29, 4);
    barY_n2(5) = csvMatrix(38, 4);
    
    barY_n3(1) = csvMatrix(3, 4);
    barY_n3(2) = csvMatrix(12, 4);
    barY_n3(3) = csvMatrix(21, 4);
    barY_n3(4) = csvMatrix(30,4);
    barY_n3(5) = csvMatrix(39, 4);
    
    barY_n4(1) = csvMatrix(4, 4);
    barY_n4(2) = csvMatrix(13, 4);
    barY_n4(3) = csvMatrix(22, 4);
    barY_n4(4) = csvMatrix(31, 4);
    barY_n4(5) = csvMatrix(40, 4);
    
    barY_n5(1) = csvMatrix(5, 4);
    barY_n5(2) = csvMatrix(14, 4);
    barY_n5(3) = csvMatrix(23, 4);
    barY_n5(4) = csvMatrix(32, 4);
    barY_n5(5) = csvMatrix(41, 4);
    
    barY_n6(1) = csvMatrix(6, 4);
    barY_n6(2) = csvMatrix(15, 4);
    barY_n6(3) = csvMatrix(24, 4);
    barY_n6(4) = csvMatrix(33, 4);
    barY_n6(5) = csvMatrix(42, 4);
    
    barY_n7(1) = csvMatrix(7, 4);
    barY_n7(2) = csvMatrix(16, 4);
    barY_n7(3) = csvMatrix(25, 4);
    barY_n7(4) = csvMatrix(34, 4);
    barY_n7(5) = csvMatrix(43, 4);
    
    barY_n8(1) = csvMatrix(8, 4);
    barY_n8(2) = csvMatrix(17, 4);
    barY_n8(3) = csvMatrix(26, 4);
    barY_n8(4) = csvMatrix(35, 4);
    barY_n8(5) = csvMatrix(44, 4);
    
    barY_n9(1) = csvMatrix(9, 4);
    barY_n9(2) = csvMatrix(18, 4);
    barY_n9(3) = csvMatrix(27, 4);
    barY_n9(4) = csvMatrix(36, 4);
    barY_n9(5) = csvMatrix(45, 4);
    
    barY_n10(1) = csvMatrix(10, 4);
    barY_n10(2) = csvMatrix(19, 4);
    barY_n10(3) = csvMatrix(28, 4);
    barY_n10(4) = csvMatrix(37, 4);
    barY_n10(5) = csvMatrix(46, 4);
    
    
    %Nach Komponenten sortiert
    barY_c3(1) = csvMatrix(2, 4);
    barY_c3(2) = csvMatrix(3, 4);
    barY_c3(3) = csvMatrix(4, 4);
    barY_c3(4) = csvMatrix(5, 4);
    barY_c3(5) = csvMatrix(6, 4);
    barY_c3(6) = csvMatrix(7, 4);
    barY_c3(7) = csvMatrix(8, 4);
    barY_c3(8) = csvMatrix(9, 4);
    barY_c3(9) = csvMatrix(10, 4);
    
    
    barY_c4(1) = csvMatrix(11, 4);
    barY_c4(2) = csvMatrix(12, 4);
    barY_c4(3) = csvMatrix(13, 4);
    barY_c4(4) = csvMatrix(14, 4);
    barY_c4(5) = csvMatrix(15, 4);
    barY_c4(6) = csvMatrix(16, 4);
    barY_c4(7) = csvMatrix(17, 4);
    barY_c4(8) = csvMatrix(18, 4);
    barY_c4(9) = csvMatrix(19, 4);
    
    barY_c5(1) = csvMatrix(20, 4);
    barY_c5(2) = csvMatrix(21, 4);
    barY_c5(3) = csvMatrix(22, 4);
    barY_c5(4) = csvMatrix(23, 4);
    barY_c5(5) = csvMatrix(24, 4);
    barY_c5(6) = csvMatrix(25, 4);
    barY_c5(7) = csvMatrix(26, 4);
    barY_c5(8) = csvMatrix(27, 4);
    barY_c5(9) = csvMatrix(28, 4);
    
    barY_c6(1) = csvMatrix(29, 4);
    barY_c6(2) = csvMatrix(30, 4);
    barY_c6(3) = csvMatrix(31, 4);
    barY_c6(4) = csvMatrix(32, 4);
    barY_c6(5) = csvMatrix(33, 4);
    barY_c6(6) = csvMatrix(34, 4);
    barY_c6(7) = csvMatrix(35, 4);
    barY_c6(8) = csvMatrix(36, 4);
    barY_c6(9) = csvMatrix(37, 4);
    
    
    barY_c7(1) = csvMatrix(38, 4);
    barY_c7(2) = csvMatrix(39, 4);
    barY_c7(3) = csvMatrix(40, 4);
    barY_c7(4) = csvMatrix(41, 4);
    barY_c7(5) = csvMatrix(42, 4);
    barY_c7(6) = csvMatrix(43, 4);
    barY_c7(7) = csvMatrix(44, 4);
    barY_c7(8) = csvMatrix(45, 4);
    barY_c7(9) = csvMatrix(46, 4);
    
    
    %Nach Komponenten sortiert und Differenz berechnet
    barY_c3_diff(1) = barY_c3(1) - barY_c3(2);
    barY_c3_diff(2) = barY_c3(2) - barY_c3(3);
    barY_c3_diff(3) = barY_c3(3) - barY_c3(4);
    barY_c3_diff(4) = barY_c3(4) - barY_c3(5);
    barY_c3_diff(5) = barY_c3(5) - barY_c3(6);
    barY_c3_diff(6) = barY_c3(6) - barY_c3(7);
    barY_c3_diff(7) = barY_c3(7) - barY_c3(8);
    barY_c3_diff(8) = barY_c3(8) - barY_c3(9);
    
    
    barY_c4_diff(1) = barY_c4(1) - barY_c4(2);
    barY_c4_diff(2) = barY_c4(2) - barY_c4(3);
    barY_c4_diff(3) = barY_c4(3) - barY_c4(4);
    barY_c4_diff(4) = barY_c4(4) - barY_c4(5);
    barY_c4_diff(5) = barY_c4(5) - barY_c4(6);
    barY_c4_diff(6) = barY_c4(6) - barY_c4(7);
    barY_c4_diff(7) = barY_c4(7) - barY_c4(8);
    barY_c4_diff(8) = barY_c4(8) - barY_c4(9);
    
    
    barY_c5_diff(1) = barY_c5(1) - barY_c5(2);
    barY_c5_diff(2) = barY_c5(2) - barY_c5(3);
    barY_c5_diff(3) = barY_c5(3) - barY_c5(4);
    barY_c5_diff(4) = barY_c5(4) - barY_c5(5);
    barY_c5_diff(5) = barY_c5(5) - barY_c5(6);
    barY_c5_diff(6) = barY_c5(6) - barY_c5(7);
    barY_c5_diff(7) = barY_c5(7) - barY_c5(8);
    barY_c5_diff(8) = barY_c5(8) - barY_c5(9);
    
    
    barY_c6_diff(1) = barY_c6(1) - barY_c6(2);
    barY_c6_diff(2) = barY_c6(2) - barY_c6(3);
    barY_c6_diff(3) = barY_c6(3) - barY_c6(4);
    barY_c6_diff(4) = barY_c6(4) - barY_c6(5);
    barY_c6_diff(5) = barY_c6(5) - barY_c6(6);
    barY_c6_diff(6) = barY_c6(6) - barY_c6(7);
    barY_c6_diff(7) = barY_c6(7) - barY_c6(8);
    barY_c6_diff(8) = barY_c6(8) - barY_c6(9);
    
    
    
    barY_c7_diff(1) = barY_c7(1) - barY_c7(2);
    barY_c7_diff(2) = barY_c7(2) - barY_c7(3);
    barY_c7_diff(3) = barY_c7(3) - barY_c7(4);
    barY_c7_diff(4) = barY_c7(4) - barY_c7(5);
    barY_c7_diff(5) = barY_c7(5) - barY_c7(6);
    barY_c7_diff(6) = barY_c7(6) - barY_c7(7);
    barY_c7_diff(7) = barY_c7(7) - barY_c7(8);
    barY_c7_diff(8) = barY_c7(8) - barY_c7(9);
    
    
    
    figure;
    barY_ncomp = [barY_c3; barY_c4; barY_c5; barY_c6; barY_c7];
    titles_ncomp = categorical({'NC=3', 'NC=4', 'NC=5', 'NC=6', 'NC=7'});
    bar(titles_ncomp, barY_ncomp);
    legend('N = 2','N = 3', 'N = 4', 'N = 5', 'N = 6', 'N = 7', 'N = 8', 'N = 9', 'N = 10');
    
    figure;
    barY_ncomp = [barY_c3_diff; barY_c4_diff; barY_c5_diff; barY_c6_diff; barY_c7_diff];
    titles_ncomp = categorical({'NC=3', 'NC=4', 'NC=5', 'NC=6', 'NC=7'});
    bar(titles_ncomp, barY_ncomp);
    legend('N = 2','N = 3', 'N = 4', 'N = 5', 'N = 6', 'N = 7', 'N = 8', 'N = 9', 'N = 10');
    
    figure;
    barY_n = [barY_n2; barY_n3; barY_n4; barY_n5; barY_n6; barY_n7; barY_n8; barY_n9];
    titles_n = categorical({'N = 2', 'N = 3', 'N = 4', 'N = 5', 'N = 6', 'N = 7', 'N = 8', 'N = 9'});
    bar(titles_n, barY_n);
    legend('NC = 3', 'NC = 4', 'NC = 5', 'NC = 6', 'NC = 7');
end
