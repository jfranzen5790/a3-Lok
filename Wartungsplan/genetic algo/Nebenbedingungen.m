function [c,ceq] = mycon(x)
    %c = calcG(x) - 0.05;
    
    
    global A G EWC CP n C numcomp decisionfactor skalierung arbeitszeit trans transmat TGESh Cmax;  
    
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
    
    
    %1. Nebenbedingung:
    %zustandsabhängige Wartungsvorgabe auf Basis des Erwartungswertes der
    %reaktiven Kosten zum Zeitpunkt t(x)
    praeventive_kosten_procomp = zeros(1, numcomp);
    reaktive_kosten_procomp = zeros(1, numcomp);
    for i=1:numcomp
        for j=1:n
            reaktive_kosten_procomp(i) = max(reaktive_kosten_procomp(i), EWC(i,j));
        end
        praeventive_kosten_procomp(i) = decisionfactor*CP(1,i);
    end
    
    abweichungen_procomp = reaktive_kosten_procomp(:) - praeventive_kosten_procomp(:);
    max_kostenabweichung = max(abweichungen_procomp);


    %2. Nebenbedingung:
    %Stillstand zwischen zwei Wartungen nicht kleiner als 3 Monate
    arbeitstage_pro_monat = 21;
    monate_pro_zeitschritt = skalierung / (arbeitszeit * arbeitstage_pro_monat);
    monate_pro_zeitschritt = round(monate_pro_zeitschritt);
    
    wartezeit_bis_wartung_in_monaten = monate_pro_zeitschritt * n; % Max. mögliche Wartezeit, wenn nie gewartet wird
    
    %Testwert bei ncomp = 7 und n =25
    %Ared = [0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 30 0;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];
    
    for i = 1 : n
        %ausfalldauer_ab_wartung
        if Ared(1,i) > 0 %Prüfen, ob eine Wartung vorliegt
            
            %Vorzeitiges Abbrechen ermöglichen, da in diesem Fall bereits
            %ein Min. erreicht wird
            for j = 1 : Ared(2, i) -1 % Für die Dauer der aktuellen Wartung prüfen (Wichtig für Dauer>1)
                if Ared(1, i + j) > 0
                    %Es liegt bereits während der aktuellen Wartung eine
                    %neue Wartung an -> vollkommen ungültig
                    wartezeit_bis_wartung_in_monaten = 0;
                    break;
                end
            end
            
            %Abbruchkriterium erfüllt
            if wartezeit_bis_wartung_in_monaten == 0
              break;  
            end
            %Abbruchkriterium nicht erfüllt, dann bis zur nächsten Wartung
            %prüfen
            vorheriger_wartungszeitpunkt = i;
            for j = vorheriger_wartungszeitpunkt + 1 : n
                if Ared(1, j) > 0
                    wartezeit_bis_wartung_in_monaten = min(wartezeit_bis_wartung_in_monaten, (j-i - 1) * monate_pro_zeitschritt);
                    break;
                end
            end
        end
    end
    
    min_wartezeit_anforderung_in_monaten = 3;
    wartezeit_grenzwert_abweichung = min_wartezeit_anforderung_in_monaten - wartezeit_bis_wartung_in_monaten;
    
    wartezeit_grenzwert_abweichung = wartezeit_grenzwert_abweichung * 10000;
    
    if((max_kostenabweichung <= 0) && (wartezeit_grenzwert_abweichung <= 0))
        disp("Beide Nebenbedingungen sind erfüllt!");
    elseif(max_kostenabweichung <= 0) && (wartezeit_grenzwert_abweichung ~= 0)
        disp("Nur Nebenbedingungmit EW erfüllt!");
            
    elseif((max_kostenabweichung ~= 0) && (wartezeit_grenzwert_abweichung <= 0))
        disp("Nur Nebenbedingung mit 3 Monaten Wartezeit erfüllt!");
    
    end
                
    
    c = [max_kostenabweichung; wartezeit_grenzwert_abweichung];

    ceq =  [];% Compute nonlinear equalities at x.
end