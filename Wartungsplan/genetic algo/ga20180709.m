global skalierung numcomp asize cell_ci CGES n C transmat G EWC A TGESh;

%Anzahl Zeitschritte
n=5;
skalierung=2000; %Betriebsstunden pro Zeitschritt
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
ergwartung=ones(1,n);
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
% lb(1,n+1)=0;
ub=ones(1,n)*asize(1,2);
ub(1,1)=1;
lb(1,n)=asize(1,2);

nonlcon=[];
IntCon=1:n;

avai=0.7;
% options = gaoptimset('PlotFcns',@cost);


options = optimoptions('ga', 'PlotFcn', @gaplotbestf, 'MaxGenerations', 10000, 'MaxStallGenerations', 1000, 'FunctionTolerance', 1e-12);
% options.InitialPopulationRange = [1 8];
[X,fval,exitflag,output,population] = ga(@cost,n,[],[],[],[],lb, ub, nonlcon, IntCon, options);

%% Auswertung des Optimierungsergebnisses
% Ergebnis in Vektor und Wartungsfolge wandeln
Xint=uint8(X);
erg=zeros(numcomp,length(Xint));
for j1=1:length(Xint)-1
    ergvec=cell_ci{1,Xint(j1)};
    for j2=1:numcomp
        erg(j2,j1)=ergvec(1,j2);
    end
end
% Kosten und Ausfallwahrscheinlichkeiten für Zeitschritte bestimmen

        for j3=1:n
            for j4=1:numcomp
                if j3 == 1
                    if erg(j4,j3)== 1
                        ergG(j4,j3)=0;
                        ergA(j4,j3) = MTTRD(1,j4);
                    elseif erg(j4,j3) == 0
                        ergG(j4,j3)=calcG(j3*skalierung,j4);
                        ergA(j4,j3) = 0;
                    end
                elseif j3 >= 2
                    if erg(j4,j3)== 1
                        ergG(j4,j3)=0;
                        ergwartung(1,j4)=j3;
                        ergA(j4,j3) = MTTRD(1,j4);
                    elseif erg(j4,j3) == 0
                        ergG(j4,j3)=calcG((j3-ergwartung(1,j4))*skalierung,j4);
                    end
                    ergEWC(j4,j3)=ergG(j4,j3)*(CEL(1,j4)+CR(1,j4))/CP(j4);
                end
            end
        end
        for j5=1:n
            ergAred(1,j5)=max(ergA(:,j5));
            ergAred(2,j5)=ergAred(1,j5)*8;
        end
        avai=(TGESh-sum(ergAred(2,:)))/TGESh;
        disp(avai);

        
%% Print des Ergebnisses        
X0=linspace(0,(n-1)*skalierung,n);
figure
bar(X0,erg', 'stacked');

% function f=objfcn(x)
%     f = (x(1) + x(2))^2 - x(3);
% end

% %% Outputfkt.
% 
% function stop = outfun()
%     stop = false;
% end


%% Zielfkt. Optimierung
function y=cost(x)
    global A G EWC CEL CR asize cell_ci CP CGES n C numcomp transmat skalierung MTTRD TGESh;    
    %Umwandeln in Vektor
    transmat=zeros(numcomp, n);
    G=zeros(numcomp, n);
    A=zeros(numcomp, n);
    EWC=zeros(numcomp, n);
    trans=cell(1,n);
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
                C(1,i2) = CGES(1,i3);
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
                    G(i5,i4)=calcG(i4*skalierung,i5);
                elseif transmat(i5,i4) == 1
                    G(i5,i4) = 0;
                end
            elseif i4 >= 2
                if transmat(i5,i4) == 0
                    G(i5,i4)=calcG((i4-wartung(1,i5))*skalierung,i5);
%                     G(i5,i4)=calcG((i4)*skalierung,i5);
                elseif transmat(i5,i4) == 1
                    G(i5,i4) = 0;
                    wartung(1,i5)=i4;
                end
            end
            
%             if G(i5,i4) == 0
%                 EWC(i5,i4) = 0;
%             elseif G(i5,i4) > 0
%                 EWC(i5,i4)=G(i5,i4)*(CEL(1,i5)+CR(1,i5));
%             end
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
                if transmat(i6,i8) == 0
                    A(i6,i8)=0;
                elseif transmat(i6,i8) == 1
                    A(i6,i8)=MTTRD(1,i6);
                end
        end
    end
    for j1=1:numcomp
        for j2=1:n
            if EWC(j1,j2)>=CP(1,j1)
                C(1,j2)=10^7;
            end
        end
    end
    Ared=zeros(2,n);
    for j5=1:n
        Ared(1,j5)=max(A(:,j5));
        %Zeitschritte umrechnen
        Ared(2,j5)=ceil(Ared(1,j5)*8/skalierung);
    end
    avail=(TGESh-sum(Ared(1,:))*8)/TGESh;

%     nahe beeinanderliegende Lösung rausfiltern
%     for j6=1:n-max(Ared(2,:))
%         if sum(transmat(:,j6+1:j6+max(Ared(2,:)))) > 0
%             C(1,j6)=10^7;          
%         end
%     end
    y = sum(C(1,:));
end


%% Fkt. für Best. der Ausfallwahrscheinlichkeit zum Zeitpunkt x
function aw = calcG(x,i)
    global T;
    aw = 1-exp(-x/T(i));
end