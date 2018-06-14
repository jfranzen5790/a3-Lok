addpath("FailureExpectations");
input_nonlinear;
%Anzahl Zeitschritte
n=21;
skalierung=1000; %Betriebsstunden pro Zeitschritt
X0=linspace(0,(n-1)*skalierung,n);

% Entscheidungsvariable, wann gewartetet werden soll
% Ist Maß für die Risikobereitschaft des Entscheiders
drisk=1*CP(1);
%Array für Wartungsfolgen (erschöpfend), j:Anzahl Komponenten
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
EWCEL=zeros(m,n,j);
for jj=1:2 %TODO: Hier nicht j statt 2?
    for ll=1:m
        %Zurücksetzen der Variable für jede Zeile
        indexNachWartung = 1;
        for mm=1:n
            %Wenn erster Zeitschritt Wartung, Ausfallwahrscheinlichkeit = 0;
            if mm == 1
                if mplan(ll,mm,jj) == 1
                    if mplan(ll,n,jj) == 0
                        mplan(ll,:,jj)=1;
                    elseif mplan(ll,n,jj) == 1
                    G(ll,mm,jj)=0;
                    EWC(ll,mm,jj)=0;
                    EWCEL(ll,mm,jj)=G(ll,mm,jj)*CEL(jj);
                    end
                elseif mplan(ll,mm,jj) == 0
                    mplan(ll,:,jj)=1;
    %                 G(ll,mm)=calcG(mm*skalierung);
    %                 EWC(ll,mm)=G(ll,mm)*CR(1);
    %                 EWCEL(ll,mm)=G(ll,mm)*CEL(1);
                end
            %Ab zweitem Zeitschritt    
            elseif mm >= 2
                %Für Zeitpunkt 1 überspringen, da kein Vorgänger vorhanden
                %(und dementsprechend auch kein Bedarf zum Resetten)
                if mm-1 > 0
                    %Wenn in einem Zeitschritt zuvor eine Wartung war, wird
                    %der Index auf 1 gesetzt
                    if G(ll, mm-1, jj) == 0
                            indexNachWartung = 1;
                    end
                end
                if mplan(ll,mm,jj)==1
                    G(ll,mm,jj)=0;
                    EWC(ll,mm,jj)=0;
                    EWCEL(ll,mm,jj)=0;
                elseif mplan(ll,mm,jj)==0
                    % Mit dem Index den richtigen Wert nach der letzten Wartung berechnen
                    %TODO: der letzte Wert "1" stellt den Komponentenindex
                    %dar. Also mit j ersetzen, richtig?
                    G(ll,mm,jj)=calcG(indexNachWartung*skalierung, n, skalierung, 1);
                    % G(ll,mm)=G(ll,mm-1)+0.01*skalierung;
                    EWC(ll,mm,jj)=G(ll,mm,jj)*CR(jj);
                    EWCEL(ll,mm,jj)=G(ll,mm,jj)*CEL(jj);
                    indexNachWartung = indexNachWartung + 1;
                end           
            end
            % Wenn Ausfallwahrscheinlichkeit innerhalb der 
            % Simulation > 1, dann komplette Reihe = 1 (0 <= G <= 1)    
    %         if G(ll,mm)>1
    %             mplan1(ll,:)=1;
    %         end
            if EWC(ll,mm,jj)+EWCEL(ll,mm,jj)>=drisk
                mplan(ll,:,jj)=1;
            end

        end
    end
    %Kosten durch Stillstand
    A=zeros(m,n,jj);
    A(:,:,jj)=mplan(:,:,jj)*CEL(jj);
end
%Kosten der Wartungsfolgen bestimmen
CGES=zeros(m,1);
for kk=1:m
        CGES(kk,1)=sum(mplan(kk,:,1))*(CWS+CP(1))+sum(A(kk,:,1));
end

ergebnis=min(CGES);
k=find(CGES==ergebnis);
% Stillstand durch präventive Wartung
% A=zeros(length(k),1);
% for oo=1:length(k)
%     pp=k(oo);
%     A(oo)=sum(mplan1(pp,:))*MTTRD(1);
% end
Y0perm=permute(mplan(k(1),:,2),[3 1 2]);
Y0=shiftdim(Y0perm,1);
G0perm=permute(G(k(1),:,2),[3 1 2]);
G0=shiftdim(G0perm,1);
EWC0=shiftdim(EWC(k(1),:,1),1);
figure
subplot(2,1,1);
plot(X0,G0);
axis([0 n*skalierung 0 max(G0)]);
subplot(2,1,2);
bar(X0,Y0);
axis([0 n*skalierung 0 1.5]);

function y = calcG(x, n, skalierung, componentIndex)
    %TODO: Gamma-Berechnung nach input_nonlinear auslagern und mit
    %Komponenten koppeln
    %Dann werden die Übergabeparameter n und skalierung überflüssig, da der
    %componentIndex dafür ausreicht
    gamma1 = calculate_gamma(n *skalierung + 1, "sinkend");
    y = calculate_FailureExpectation(x, gamma1(x), "kugellager", componentIndex);
    
end
%% noch zu tun:
% 3. mehrere Komponenten