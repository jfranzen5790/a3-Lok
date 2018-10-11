function [transcell,transmat2] = calctransmat(x)
% Zeitstrahl berechnen
    global n numcomp cell_ci;
    transcell=cell(1,n);
    transmat2=zeros(numcomp,n);
    for i=1:n
        transcell{1,i}=cell_ci{1,x(i)};
        transvec=cell_ci{1,x(i)};
        for j=1:numcomp
            transmat2(j,i)=transvec(1,j);        
        end   
    end
end