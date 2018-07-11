A=cell(length(row),n);
combDT=zeros(length(row),n);
combDTsum=zeros(length(row),1);
for i1=1:length(row)
    for i2=1:n
        A{i1,i2}=combmplan{row(i1),i2};
        for i8=1:length(cell_ci)
            if ismember(A{i1,i2},cell_ci{1,i8})
                combDT(i1,i2)=combMTTRDmax(1,i8);
            end
        end
    end
    combDTsum(i1,1)=sum(combDT(i1,:));
end
minDT=min(combDTsum);
row2=find(minDT==combDTsum);
% Verfügbarkeit bestimmen
av1=(TGESd-combDTsum(row2(1)))/TGESd;