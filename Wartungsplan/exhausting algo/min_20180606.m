global PopSize skalierung numcomp asize cell_ci CGES n Cmax arbeitszeit transmat G EWC A TGESh allCosts costList gueltigCounter ungueltigCounter bestCost bestComb bestInput iterationsTillNewValidValue iterationsTillNewValidValue_ListOfallValues longestRuntime;
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
        

        % Entscheidungsvariable, wann gewartetet werden soll
        % Ist Maß für die Risikobereitschaft des Entscheiders
        drisk=1*CP;
        %Array für Wartungsfolgen (erschöpfend), j:Anzahl Komponenten
        for j=1:numcomp    
            L = ones(1,n);
            L(2,:) = zeros(1,length(L)); %// example input
            L2 = mat2cell(L.', ones(1,size(L,2)), size(L,1)); %// step 1
            if j==1
            mplan=shiftdim(combvec(L2{:})',-1);
            mplan=permute(mplan,[2 3 1]);
            m=length(mplan);
            elseif j>=1
            mplan(:,:,j)=mplan(:,:,1);
            end    
        end
        mplanred = mplan;
        % %Zufällige Wartungsfolgen erzeugen
        % for ii=1:m
        %     for jj=1:n
        %       if mplan(ii,jj)>=0.5
        %           mplan(ii,jj)=1;
        %       elseif mplan(ii,jj)<0.5
        %           mplan(ii,jj)=0;
        %       end
        %     end
        % end

        %Ausfallwahrscheinlichkeiten G und Erwartungswert der Kosten EWC
        G=zeros(m,n,j);
        EWC=zeros(m,n,j);
        EWCGES=zeros(m,n);
        EWCEL=zeros(m,n,j);
        q10=zeros(1,numcomp);
        for jj=1:numcomp
            for ll=1:m
                for mm=1:n
                    %Wenn erster Zeitschritt Wartung, Ausfallwahrscheinlichkeit = 0;
                    if mm == 1
                        if mplan(ll,mm,jj) == 1
        %                     if mplan(ll,n,jj) == 0
        %                         mplan(ll,:,jj)=1;
        %                     elseif mplan(ll,n,jj) == 1
                            G(ll,mm,jj)=0;
                            EWC(ll,mm,jj)=0;
                            EWCEL(ll,mm,jj)=G(ll,mm,jj)*(CEL(jj)+CR(jj));
        %                    end
                        elseif mplan(ll,mm,jj) == 0
                            for b1=1:n%TODO: 
                                mplan(ll,b1,jj)=1;
                            end
            %                 G(ll,mm)=calcG(mm*skalierung);
            %                 EWC(ll,mm)=G(ll,mm)*CR(1);
            %                 EWCEL(ll,mm)=G(ll,mm)*CEL(1);
                        end
                    %Ab zweitem Zeitschritt    
                    elseif mm >= 2
                        if mplan(ll,mm,jj)==1
                            G(ll,mm,jj)=0;
                            EWC(ll,mm,jj)=0;
                            EWCEL(ll,mm,jj)=0;
                        elseif mplan(ll,mm,jj)==0
                            % Differenz der Wahrscheinlichkeit im Intervall addieren
                            G(ll,mm,jj)=G(ll,mm-1,jj)+calcG((mm-(mm-1))*skalierung,jj);
                            % G(ll,mm)=G(ll,mm-1)+0.01*skalierung;
                            EWC(ll,mm,jj)=G(ll,mm,jj)*CR(jj);
                            EWCEL(ll,mm,jj)=G(ll,mm,jj)*CEL(jj);
                        end           
                    end
                    %Wenn Ausfallwahrscheinlichkeit innerhalb der 
                    % Simulation >= 1 dann komplette Reihe = 1 (0 <= G <= 1)    
                    if G(ll,mm,jj)>1
                        for b2=1:n
                           %TODO: mplan(ll,b2,jj)=Inf;
                           mplan(ll,b2,jj)=1;
                        end
                    end
                    % Wenn Erwartungswert der reaktiven Kosten den Schwellwert
                    % überschreitet, dann Zeile unbrauchbar machen
                    if EWC(ll,mm,jj)+EWCEL(ll,mm,jj)>=drisk(jj)
                        for b3=1:n%TODO: 
                            mplan(ll,b3,jj)=1;
                        end
                    end
                end
            end
        %Kosten durch Stillstand
        %     A=zeros(m,n,jj);
        %     A(:,:,jj)=mplan(:,:,jj)*CEL(jj);
            % reduziertes Array mplan mplanred ohne Zeilen mit inf erstellen
            for i3=1:length(mplan(:,:,jj))%TODO: 
                if mplan(i3,:,jj) == 1
                    q10(1,jj)=q10(1,jj)+1;
                end
            end
            q10(1,jj)=length(mplan(:,:,jj))-q10(1,jj);
        %     for c1=1:m
        %         if mplan(c1,:,jj) ~= Inf
        %             for c2=1:n
        %                     mplanred(c3,c2,jj)=mplan(c1,c2,jj);
        %             end
        %             c3=c3+1;
        %         end
        %     end


        end

        %% Ursprüngliche Berechnung von mplanred
        %mplanred=ones(max(q10),n,numcomp);
        %for jj2=1:numcomp
        %    c3=1;
        %    for c1=1:length(mplan(:,:,jj2))
        %        if mplan(c1,:,jj2) ~= Inf
        %            for c2=1:n
        %                    mplanred(c3,c2,jj2)=mplan(c1,c2,jj2);
        %            end
        %            c3=c3+1;  
        %            gueltigCounter = gueltigCounter + 1;
        %        else
        %            ungueltigCounter = ungueltigCounter + 1;    
        %        end
        %    end
        %end
        %mplanred = mplan;
        
        %% Neue reduzierte Berechnung von mplanred für größere n-Werte und ncomp Werte
        
       %randomCombination = zeros(1, n);
       %mountValues = 1000000;
        
        %mplanred = zeros(amountValues, n);
        %mplanredCounter = 1;
        
        %while mplanredCounter < amountValues
        %    for t = 1:n
        %       randomCombination(t) = randperm(2^numcomp,1);
        %    end 
        %   if   ~isempty(mplanred)
        %        if  ~any(isequal(randomCombination, mplanred))
        %            mplanred(mplanredCounter, :) = randomCombination;
        %            mplanredCounter = mplanredCounter + 1;
        %        else
        %           disp("Else"); 
        %        end
        %   end
        %end
        
           
        
        
        
        
        
        
        
        
        
        %gueltigRate = gueltigCounter / ungueltigCounter;
        %disp("Rate der gültigen Werte: " + gueltigRate);

        %csvMatrix = csvread('gueltigRaten.dat');

        %csvQuadruplet = [numcomp, n, skalierung, gueltigRate];
        %csvMatrix(end + 1, :) = csvQuadruplet;
        %csvwrite('gueltigRaten.dat', csvMatrix);
        %type gueltigRaten.dat;
%    end
%end

%% CSV-Abschnitt Ende

%%CSV-Plot

%csvMatrix = csvread('gueltigRaten.dat');
%[X,Y] = meshgrid(400:400:2000,1:20);
%for rangeN = linspace(20,1,20)
%for rangeS = 1:5
%Z(rangeN, rangeS) = csvMatrix(100 + rangeS + (rangeN -1 )* 5, 4);
%end
%end
%surf(X,Y,Z);
%%surf(X,Y,fliplr(flipud(Z)));
%surf(X,Y,fliplr(Z));

%%CSV-Plot-Ende

% Erwartungswert für Gesamtkosten durch spontanen Ausfall für 
% alle Komponenten bestimmen
for ii=1:m
    for nn=1:n    
        EWCGES(ii,nn)=sum(EWC(ii,nn,:));
    end
end
%Kosten der Wartungsfolgen bestimmen
% Cmplanstep=zeros(m,n,m);
% Kosten für jeweilige Wartungskombination bestimmen
intervallC_20180611;
% Matrix erzeugen, in der alle Kombinationen beider Wartungspläne stehen
% for bb=1:m
%     for cc=1:m
%         for dd=1:n % Jede Zelle kombinieren
%             if mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 0
%                 Cmplanstep(cc,dd,bb)=CIS(1)+ MTTR(1)*M + MTTRD(1)*CELd;
%             elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 1
%                 Cmplanstep(cc,dd,bb)=CIS(2)+ MTTR(2)*M + MTTRD(2)*CELd;
%             elseif mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 1
%                 Cmplanstep(cc,dd,bb)=CIS(1) + CIS(2) + (MTTR(1)+MTTR(2))*M + max(MTTRD)*CELd;
%             elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 0
%                 Cmplanstep(cc,dd,bb)=0;
%             elseif mplan(bb,dd,1) <= 1 && mplan(cc,dd,2) == inf
%                 Cmplanstep(cc,dd,bb)=inf;
%             elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= 1
%                 Cmplanstep(cc,dd,bb)=inf;
%             elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= inf
%                 Cmplanstep(cc,dd,bb)=inf;
%             end
%             
%         end
%         
%         CGES(bb,cc)=sum(Cmplanstep(cc,:,bb));
%     end
% end
% for kk=1:(m*m)
%     for dd=1:m
%         for cc=1:m
%             for bb=1:n
%                 if mplan(kk,bb,1) == 1 && mplan(cc,bb,2) == 0 
%                     Cmplanstep(dd,bb)=CIS(1)+ MTTR(1)*M + MTTRD(1)*CELd;
%                 elseif mplan(kk,bb,1) == 0 && mplan(kk,cc,2) == 1 
%                     Cmplanstep(dd,bb)=CIS(2)+ MTTR(2)*M + MTTRD(2)*CELd;
%                 elseif mplan(kk,bb,1) == 1 && mplan(kk,cc,2) == 1 
%                     Cmplanstep(dd,bb)=CIS(1) + CIS(2) + (MTTR(1)+MTTR(2))*M + max(MTTRD(:))*CELd;
%                 elseif mplan(kk,bb,1) == 0 && mplan(kk,cc,2) == 0
%                     Cmplanstep(dd,bb)=0;
%                 end
% 
%         %         CGES(kk,1)=sum(mplan(kk,:,1))*(CWS+CP(1))+sum(A(kk,:,1));
%             end
%     %         Cmplanstep(kk,1)=sum(mplan(kk,:,1))*(CWS+CP(1))+sum(A(kk,:,1));
%         CGES(kk,1)=sum(Cmplanstep(kk,:));
%         end
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
Y01perm=cell2mat(bestInput);
Y0=zeros(length(Y01perm)/numcomp,numcomp);
for i5=1:n
    for i7=1:numcomp
        Y0(i5,i7)=Y01perm(1,(i5-1)*numcomp+i7);
    end
end
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
figure
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
subplot(3,1,3);
bar(X0,Y0, 'stacked');
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
    %patch('XData', [drawvec(i)-1 drawvec(i)+ergAred(2,drawvec(i))-1 drawvec(i)+ergAred(2,drawvec(i))-1 drawvec(i)-1]*skalierung, 'YData', [0 0 numcomp numcomp], 'EdgeColor', 'none', 'FaceColor', 'red', 'FaceAlpha', 0.5);
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


csvMatrix = csvread('Rate_gültigeWerte.dat');

csvQuinttruplet = [numcomp, n, skalierung, gueltigRate, longestRuntime];
csvMatrix(end + 1, :) = csvQuinttruplet;
csvwrite('Rate_gültigeWerte.dat', csvMatrix);
type Rate_gültigeWerte.dat;

csvAuswertungAktiviert = 0;

if csvAuswertungAktiviert == 1
    %% CSV Auswerungsbereich
    csvMatrix = csvread('Rate_gültigeWerte.dat');
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
    barY_n2(2) = csvMatrix(3, 4);
    barY_n2(3) = csvMatrix(4, 4);
    barY_n2(4) = csvMatrix(5, 4);
    barY_n2(5) = csvMatrix(6, 4);
    
    barY_n3(1) = csvMatrix(7, 4);
    barY_n3(2) = csvMatrix(8, 4);
    barY_n3(3) = csvMatrix(9, 4);
    barY_n3(4) = csvMatrix(10,4);
    barY_n3(5) = csvMatrix(11, 4);
    
    barY_n4(1) = csvMatrix(12, 4);
    barY_n4(2) = csvMatrix(13, 4);
    barY_n4(3) = csvMatrix(15, 4);
    barY_n4(4) = 0;
    barY_n4(5) = 0;
    
    barY_n5(1) = csvMatrix(14, 4);
    barY_n5(2) = 0;
    barY_n5(3) = 0;
    barY_n5(4) = 0;
    barY_n5(5) = 0;
    
    
    %Nach Komponenten sortiert
    barY_c3(1) = csvMatrix(2, 4);
    barY_c3(2) = csvMatrix(7, 4);
    barY_c3(3) = csvMatrix(12, 4);
    barY_c3(4) = csvMatrix(14, 4);
    
    barY_c4(1) = csvMatrix(3, 4);
    barY_c4(2) = csvMatrix(8, 4);
    barY_c4(3) = csvMatrix(13, 4);
    barY_c4(4) = 0;
    
    barY_c5(1) = csvMatrix(4, 4);
    barY_c5(2) = csvMatrix(9, 4);
    barY_c5(4) = 0;
    barY_c5(4) = 0;
    
    barY_c6(1) = csvMatrix(5, 4);
    barY_c6(2) = csvMatrix(10, 4);
    barY_c6(4) = 0;
    barY_c6(4) = 0;
    
    
    barY_c7(1) = csvMatrix(6, 4);
    barY_c7(2) = csvMatrix(11, 4);
    barY_c7(4) = 0;
    barY_c7(4) = 0;
    
    
    %Nach Komponenten sortiert und Differenz berechnet
    barY_c3_diff(1) = 0;
    barY_c3_diff(2) = barY_c3(1) - barY_c3(2);
    barY_c3_diff(3) = barY_c3(2) - barY_c3(3);
    barY_c3_diff(4) = barY_c3(3) - barY_c3(4);
    
    barY_c4_diff(1) = 0;
    barY_c4_diff(2) = barY_c4(1) - barY_c4(2);
    barY_c4_diff(3) = barY_c4(2) - barY_c4(3);
    barY_c4_diff(4) = 0;
    
    barY_c5_diff(1) = 0;
    barY_c5_diff(2) = barY_c5(1) - barY_c5(2);
    barY_c5_diff(3) = 0;
    barY_c5_diff(4) = 0;
    
    barY_c6_diff(1) = 0;
    barY_c6_diff(2) = barY_c6(1) - barY_c6(2);
    barY_c6_diff(3) = 0;
    barY_c6_diff(4) = 0;
    
    
    barY_c7_diff(1) = 0;
    barY_c7_diff(2) = barY_c7(1) - barY_c7(2);
    barY_c7_diff(3) = 0;
    barY_c7_diff(4) = 0;
    
    
    
    figure;
    barY_ncomp = [barY_c3; barY_c4; barY_c5; barY_c6; barY_c7];
    titles_ncomp = categorical({'NC=3', 'NC=4', 'NC=5', 'NC=6', 'NC=7'});
    bar(titles_ncomp, barY_ncomp);
    legend('N = 2','N = 3', 'N = 4', 'N = 5');
    
    figure;
    barY_ncomp = [barY_c3_diff; barY_c4_diff; barY_c5_diff; barY_c6_diff; barY_c7_diff];
    titles_ncomp = categorical({'NC=3', 'NC=4', 'NC=5', 'NC=6', 'NC=7'});
    bar(titles_ncomp, barY_ncomp);
    legend('N = 2','N = 3', 'N = 4', 'N = 5');
    
    figure;
    barY_n = [barY_n2; barY_n3; barY_n4; barY_n5];
    titles_n = categorical({'N=2', 'N=3', 'N=4', 'N=5'});
    bar(titles_n, barY_n);
    legend('NC=3', 'NC=4', 'NC=5', 'NC=6', 'NC=7');
end



