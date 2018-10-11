%% Wartungsfälle durch Kombination für Zeitschritt definieren
% Alle möglichen Wartungsaktivitäten pro Zeitschritt
a=[1 0];
combintervall=combvec(a,a,a);
asize=size(combintervall);
cell_ci=cell(1,asize(1,2));
for ci=1:asize(1,2)
    cell_ci{1,ci}=[combintervall(1,ci) combintervall(2,ci) combintervall(3,ci)];
end
    % %Kosten für Komponentenwartung
combCP=zeros(asize(1,1),asize(1,2));
combCEL=zeros(1,asize(1,2));
for a1=1:asize(1,1)
    for a2=1:asize(1,2)
        if combintervall(a1,a2) == 1
            combCP(a1,a2)=CP(a1);
            combCEL(a1,a2)=CEL(a1);
        end
    end
end
% combCGES=zeros(1,asize(1,1));
% for a3=1:asize(1,2)
%     combCGES(1,a3)=sum(combC(:,a3))+max(max(combCEL(:,a3)));
% end
combCELsize=size(combCEL);
combCELmax=zeros(1,combCELsize(1,2));
for c4=1:length(combCEL)
    if sum(combCEL(:,c4))==0
        combCELmax(1,c4)=0;
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
end
combmplan=cell(numcomp*max(q10),n);
combCGES=zeros(numcomp*max(q10),n);
CGES=zeros(numcomp*max(q10),1);
i4=1;
for bb=1:length(mplanred)
    for cc=1:length(mplanred)
        for dd=1:length(mplanred)
            for ee=1:n
                vector=[mplanred(bb,ee,1) mplanred(cc,ee,2) mplanred(dd,ee,3)];
                combmplan{i4,ee}=vector;
                cell_cisize=size(cell_ci);
                for i6=1:cell_cisize(1,2)
                    check=ismember(combmplan{i4,ee},cell_ci{1,i6});
                    if sum(check)==3
                        combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
                    end
                    check=0;
                end
            end
            CGES(i4,1)=sum(combCGES(i4,:));
            i4=i4+1;
        end
    end
end


% e = cellfun('isempty', combmplan);    % Leere Zellen finden
% % combmplan(e) = []; 
% %% Kosten für Kombinationen pro Zeitschritt auf Basis des Gesamtplans definieren
% Cmplanstep=zeros(m,n,m);
% for bb=1:m
%     for cc=1:m
%         for dd=1:n % Jede Zelle kombinieren
%             for ee=1:n
%                 if mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 0 && mplan(ee,dd,3) == 0 
%                     Cmplanstep(ee,cc,dd,bb)=CIS(1)+ MTTR(1)*M + MTTRD(1)*CELd;
%                 elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 1 && mplan(ee,dd,3) == 0
%                     Cmplanstep(cc,dd,bb)=CIS(2)+ MTTR(2)*M + MTTRD(2)*CELd;
%                 elseif mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 1  && mplan(ee,dd,3) == 0
%                     Cmplanstep(cc,dd,bb)=CIS(1) + CIS(2) + (MTTR(1)+MTTR(2))*M + max(MTTRD)*CELd;
%                 elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 0 && mplan(ee,dd,3) == 0
%                     Cmplanstep(cc,dd,bb)=0;
%                 elseif mplan(bb,dd,1) <= 1 && mplan(cc,dd,2) == inf && mplan(ee,dd,3) <= 1
%                     Cmplanstep(cc,dd,bb)=inf;
%                 elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= 1
%                     Cmplanstep(cc,dd,bb)=inf;
%                 elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= inf
%                     Cmplanstep(cc,dd,bb)=inf;
%                 end
%             end
%         end
% 
%         CGES(bb,cc)=sum(Cmplanstep(cc,:,bb));
%     end
% end
