% Mögliche Wartungsaktivitäten bei einem Zeitschritt
a=[1 0];
% combintervall=combvec(combvec(a,a,a));
% asize=size(combintervall);
% cell_ci=cell(1,asize(1,2));
% for ci=1:asize(1,2)
%     cell_ci{1,ci}=[combintervall(1,ci) combintervall(2,ci) combintervall(3,ci)];
% end    

if numcomp == 1
    combintervall = a;
    asize=size(combintervall);
    cell_ci=cell(1,asize(1,2));
    for ci=1:asize(1,2)
        cell_ci{1,ci}=[combintervall(1,ci)];
    end
    
elseif numcomp==2
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