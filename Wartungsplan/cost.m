function y = cost(x)
%Zielfunktion
    global A G EWC CP n C numcomp decisionfactor skalierung trans transmat TGESh Cmax;    
    %Umwandeln in Vektor
%     transmat=zeros(numcomp, n);
%     G=zeros(numcomp, n); % Ausfallwahrscheinlichkeiten
%     A=zeros(numcomp, n); % Ausfallzeiten
%     EWC =zeros(numcomp, n); %Erwartungswert der Audfallkosten
    avhelp=1; % Hilfsvariable, um Ergebnis ungültig zu machen
%     trans=cell(1,n);
    x=round(x);
    
    %Schlechter Ansatz, deaktiviert
    %%checkIfNewValue(x);
    
    %Wartungsplan aus Eingangsvariablen erstellen und Kosten pro Schritt
    %bestimmen
    [trans, transmat] = calctransmat(x);
    C = calcCost(trans);
    % Stillstandzeiten pro Zeitschritt
    [A, Ared]=calcdt(transmat);
    
% Lebensdauern bestimmen    
    lifetime=calclifetime(transmat, Ared);    
    
    % Ausfallwahrscheinlichkeiten bestimmten
    [G,EWC] = calcGmat(transmat,lifetime,1);
    
    
    for i=1:numcomp
        for j=1:n
            if EWC(i,j)>decisionfactor*CP(1,i)
%                 C(1,j)=inf;
                avhelp=0;
            end
        end
    end
    
    gueltig=0; %Hilfsvariable zur Überprüfung, ob die Stillstandszeiten eingehalten werden
    for j5=1:n
        gueltig=0;
        dauer=0;
        if Ared(2,j5)>=1
            dauer=Ared(2,j5);
        end
        if dauer >= 2
            for j6=1:dauer -1
                if j5+j6 <= n
                    gueltig=gueltig+Ared(1,j5+j6);
                end
            end
        end
        if gueltig > 0
            avhelp=0; %Schleife abbrechen, wenn ungültige Reihe
            break;
        end
    end

    if avhelp == 0
        y=1.0;
        costCell = {convertCharsToStrings(num2str(x)), y};
        saveCostValues(costCell, 1);   
    elseif avhelp == 1
        sumC=sum(C(1,:));
        sumAred=sum(Ared(2,:));
        %Ausfallzeit bestimmen
        techav=calctechav(sumAred);
        
        y=0.7*(sumC/Cmax)+0.3*(1-techav);
        costCell = {convertCharsToStrings(num2str(x)), y};
        saveCostValues(costCell, 0);   
    end
end