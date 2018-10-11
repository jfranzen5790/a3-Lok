function [Amat,Aredmat]= calcdt(transmat2)
    global n numcomp skalierung MTTRD
    Amat=zeros(numcomp,n);
    Aredmat=zeros(2,n);
    for i=1:numcomp
        for j=1:n
            % erster Durchlauf "System neu" darf nicht berücksichtigt
                % werden
                if j == 1
                    Amat(i,j)=0;
                elseif (1 < j) && (j < n)
                    if transmat2(i,j) == 0
                        Amat(i,j)=0;
                    elseif transmat2(i,j) == 1
                        Amat(i,j)=MTTRD(1,i);
                    end
                elseif j == n
                    Amat(i,j)=0;
                end
        end
    end
    
    % Anzahl Zeitschritte Wartung pro Zeitschritt
    for i=1:n
        Aredmat(1,i)=max(Amat(:,i));
        %Zeitschritte umrechnen
        Aredmat(2,i)=ceil(Aredmat(1,i)*8/skalierung);
    end
end
