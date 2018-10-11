% Wartungsschritte und -kosten pro Zeitschritt
combmplan=cell(max(q10)^numcomp,n);
combCGES=zeros(max(q10)^numcomp,n);
CGES=zeros(max(q10)^numcomp,1);
i4=1;

if numcomp==3
    for bb=1:length(mplanred)
        for cc=1:length(mplanred)
            for dd=1:length(mplanred)
                for ee=1:n
                    vector=[mplanred(bb,ee,1) mplanred(cc,ee,2) mplanred(dd,ee,3)];
                    combmplan{i4,ee}=vector;
                    cell_cisize=size(cell_ci);
                    for i6=1:cell_cisize(1,2)
                        check=ismember(combmplan{i4,ee},cell_ci{1,i6});
                        if sum(check)==numcomp
                            combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
                        end
                        check=0;
                    end
                end
                % Kosten der Wartungsstrategie = Summe der Kosten für
                % Zeitschritte
                CGES(i4,1)=sum(combCGES(i4,:));
                i4=i4+1;
            end    
        end
    end
elseif numcomp==4
    for bb=1:length(mplanred)
        for cc=1:length(mplanred)
            for dd=1:length(mplanred)
                for icost4=1:length(mplanred)
                    for ee=1:n
                        vector=[mplanred(bb,ee,1) mplanred(cc,ee,2) mplanred(dd,ee,3) mplanred(icost4,ee,4)];
                        combmplan{i4,ee}=vector;
                        cell_cisize=size(cell_ci);
                        for i6=1:cell_cisize(1,2)
%                             check=ismember(combmplan{i4,ee},cell_ci{1,i6});
%                             if sum(check)==numcomp
%                                 combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
%                             end
                            if isequal(combmplan{i4,ee},cell_ci{1,i6})
                                combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
                            end
                        end
                    end
                % Kosten der Wartungsstrategie = Summe der Kosten für
                % Zeitschritte
                CGES(i4,1)=sum(combCGES(i4,:));
                i4=i4+1;
                end
            end    
        end
    end
end


