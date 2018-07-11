%% Wartungsfälle durch Kombination für Zeitschritt definieren
% Alle möglichen Wartungsaktivitäten pro Zeitschritt (ja/nein)
a=[1 0];

if numcomp==2
    combintervall=combvec(combvec(a,a));
    asize=size(combintervall);
    cell_ci=cell(1,asize(1,2));
    for ci=1:asize(1,2)
        cell_ci{1,ci}=[combintervall(1,ci) combintervall(2,ci)];
    end
elseif numcomp==3
    combintervall=combvec(combvec(a,a,a));
    asize=size(combintervall);
    cell_ci=cell(1,asize(1,2));
    for ci=1:asize(1,2)
        cell_ci{1,ci}=[combintervall(1,ci) combintervall(2,ci) combintervall(3,ci)];
    end
elseif numcomp==4
    combintervall=combvec(combvec(a,a,a),a);
    asize=size(combintervall);
    cell_ci=cell(1,asize(1,2));
    for ci=1:asize(1,2)
        cell_ci{1,ci}=[combintervall(1,ci) combintervall(2,ci) combintervall(3,ci) combintervall(4,ci)];
    end
end
    %Kosten für Komponentenwartung
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
calccost;
