input_nonlinear;
%Anzahl Zeitschritte
n=11;
skalierung=1000; %Betriebsstunden pro Zeitschritt
X0=linspace(0,(n-1)*skalierung,n);

% Entscheidungsvariable, wann gewartetet werden soll
% Ist Ma� f�r die Risikobereitschaft des Entscheiders
drisk=1*CP;
%Array f�r Wartungsfolgen (ersch�pfend), j:Anzahl Komponenten
for j=1:length(CIS)    
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
% %Zuf�llige Wartungsfolgen erzeugen
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
for jj=1:2
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
                    EWCEL(ll,mm,jj)=G(ll,mm,jj)*CEL(jj);
%                    end
                elseif mplan(ll,mm,jj) == 0
                    mplan(ll,:,jj)=inf;
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
                mplan(ll,1:mm,jj)=inf;
            end
            % Wenn Erwartungswert der reaktiven Kosten den Schwellwert
            % �berschreitet, dann Zeile unbrauchbar machen
            if EWC(ll,mm,jj)+EWCEL(ll,mm,jj)>=drisk(jj)
                mplan(ll,:,jj)=inf;
            end

        end
    end
    %Kosten durch Stillstand
    A=zeros(m,n,jj);
    A(:,:,jj)=mplan(:,:,jj)*CEL(jj);
end
% Erwartungswert f�r Gesamtkosten durch spontanen Ausfall f�r 
% alle Komponenten bestimmen
for ii=1:m
    for nn=1:n    
        EWCGES(ii,nn)=sum(EWC(ii,nn,:));
    end
end
%Kosten der Wartungsfolgen bestimmen
Cmplanstep=zeros(m,n,m);
CGES=zeros(m,m);

% Matrix erzeugen, in der alle Kombinationen beider Wartungspl�ne stehen
ff=1;
for bb=1:m
    for cc=1:m
        for dd=1:n % Jede Zelle kombinieren
            if mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 0
                Cmplanstep(cc,dd,bb)=CIS(1)+ MTTR(1)*M + MTTRD(1)*CELd;
            elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 1
                Cmplanstep(cc,dd,bb)=CIS(2)+ MTTR(2)*M + MTTRD(2)*CELd;
            elseif mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 1
                Cmplanstep(cc,dd,bb)=CIS(1) + CIS(2) + (MTTR(1)+MTTR(2))*M + max(MTTRD)*CELd;
            elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 0
                Cmplanstep(cc,dd,bb)=0;
            elseif mplan(bb,dd,1) <= 1 && mplan(cc,dd,2) == inf
                Cmplanstep(cc,dd,bb)=inf;
            elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= 1
                Cmplanstep(cc,dd,bb)=inf;
            elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= inf
                Cmplanstep(cc,dd,bb)=inf;
            end
            
        end
        
        CGES(bb,cc)=sum(Cmplanstep(cc,:,bb));
        ff=ff+1;
    end
end
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

ergebnis=min(min(CGES));
[column,row]=find(CGES==ergebnis);
%Umrechnen
% k_1=floor(k/m);
% k_2(1)=k(1)-k_1(1)*m;
% Stillstand durch pr�ventive Wartung
% A=zeros(length(k),1);
% for oo=1:length(k)
%     pp=k(oo);
%     A(oo)=sum(mplan1(pp,:))*MTTRD(1);
% end
Y01perm=permute(mplan(column(1),:,1),[3 1 2]);
Y01=shiftdim(Y01perm,1);
Y02perm=permute(mplan(row(1),:,2),[3 1 2]);
Y02=shiftdim(Y02perm,1);
Y0=zeros(length(Y01),2);
for y0i=1:length(Y01)
    Y0(y0i,1)= Y01(y0i);
    Y0(y0i,2)= Y02(y0i);
end
G01perm=permute(G(column(1),:,1),[3 1 2]);
G01=shiftdim(G01perm,1);
G02perm=permute(G(row(1),:,2),[3 1 2]);
G02=shiftdim(G02perm,1);
% EWC0=shiftdim(EWC(k(1),:,1),1);
figure
subplot(3,1,1);
plot(X0,G01);
if sum(G01) == 0
    limx=1;
elseif sum(G01) > 0
    limx=max(G01);
end
if sum(G02) == 0
    limx2=1;
elseif sum(G02) > 0
    limx2=max(G02);
end


axis([0 n*skalierung 0 limx]);
subplot(3,1,2);
plot(X0,G02);
axis([0 n*skalierung 0 limx2]);
subplot(3,1,3);
bar(X0,Y0, 'stacked');
axis([0 n*skalierung 0 1.5]);


function y = calcG(x,i)
    global T;
    y = 1-exp(-x/T(i));
end
