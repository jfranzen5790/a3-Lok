function C = calcCost(transcell)
    global n cell_ci CGES asize skalierung arbeitszeit;
    C=zeros(1,n);
    for i=1:n
        for j=1:asize(1,2)
            if transcell{1,i} == cell_ci{1,j}
                C(1,i)= CGES(1,j)*1.02;%^(i*(skalierung/(200*arbeitszeit))); %Zinssatz berücksichtigen
                break
            elseif transcell{1,i} ~= cell_ci{1,j}
                C(1,i) = 0;
            end
        end
    end
end