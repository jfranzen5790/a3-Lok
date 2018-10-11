function [G, EWC] = calcGmat(transmat2, lt,factC)
    global n numcomp CEL CR 
    G=zeros(numcomp,n);
    EWC=zeros(numcomp,n);
    for i=1:numcomp
%         wartung=ones(1,numcomp);
        for j=1:n
            if j==1
                if transmat2(i,j) == 0
                    G(i,j)=calcG(lt(i,j),i);
                elseif transmat2(i,j) == 1
                    G(i,j) = 0;
                end
            elseif j >= 2
                if transmat2(i,j) == 0
                    G(i,j)=calcG(lt(i,j),i);
                elseif transmat2(i,j) == 1
                    G(i,j) = 0;
                end
            end
            if G(i,j) == 0
                EWC(i,j) = 0;
            elseif G(i,j) > 0
                EWC(i,j)=G(i,j)*(CEL(1,i)+CR(1,i))/factC;
            end
        end
    end
end
