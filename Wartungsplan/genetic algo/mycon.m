function [c,ceq] = mycon(x)
    %c = calcG(x) - 0.05;
    
    
    global A G EWC CP n C numcomp decisionfactor skalierung trans transmat TGESh Cmax;  
    
    avhelp=-1; % Hilfsvariable, um Ergebnis ungültig zu machen
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
                avhelp=1;
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
            avhelp=1; %Schleife abbrechen, wenn ungültige Reihe
            break;
        end
    end
    
    
    c = avhelp;
    %c = calcG(x) - 0.05;
    
    ceq =  [];% Compute nonlinear equalities at x.
end