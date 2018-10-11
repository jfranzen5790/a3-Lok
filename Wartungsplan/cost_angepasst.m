function y = cost(x)
%Zielfunktion
    global A G EWC CP n C numcomp decisionfactor skalierung trans transmat TGESh Cmax;    
    
    
    x=round(x);

    [trans, transmat] = calctransmat(x);
    C = calcCost(trans);
    % Stillstandzeiten pro Zeitschritt
    [A, Ared]=calcdt(transmat);
    
    % Lebensdauern bestimmen    
    lifetime=calclifetime(transmat, Ared);    
    
    % Ausfallwahrscheinlichkeiten bestimmten
    [G,EWC] = calcGmat(transmat,lifetime,1);
    

    sumC=sum(C(1,:));
    sumAred=sum(Ared(2,:));
    %Ausfallzeit bestimmen
    techav=calctechav(sumAred);
    
    if techav < 0
        techav = 0;
    end

    y=0.7*(sumC/Cmax)+0.3*(1-techav);
    costCell = {convertCharsToStrings(num2str(x)), y};
    saveCostValues(costCell, 0);   
end