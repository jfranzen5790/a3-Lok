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
elseif numcomp >= 5
%     combintervall=combvec(combvec(combvec(a,a,a),a),a);

    for i=1:(2^numcomp)
        tempd2b=dec2bin(i-1,numcomp);
        for j=1:length(tempd2b)
            test=tempd2b(j);
            combintervall(j,i)=bin2dec(test);
        end
    end
    combintervall=fliplr(combintervall);
    asize=size(combintervall);
    cell_ci=cell(1,asize(1,2));

    for i=1:asize(1,2)
        if numcomp == 5
        cell_ci{1,i}=[combintervall(1,i) combintervall(2,i) combintervall(3,i) combintervall(4,i) combintervall(5,i)];
        elseif numcomp == 6
        cell_ci{1,i}=[combintervall(1,i) combintervall(2,i) combintervall(3,i) combintervall(4,i) combintervall(5,i) combintervall(6,i)];
        elseif numcomp == 7
        cell_ci{1,i}=[combintervall(1,i) combintervall(2,i) combintervall(3,i) combintervall(4,i) combintervall(5,i) combintervall(6,i) combintervall(7,i)];
        end
    end
end

