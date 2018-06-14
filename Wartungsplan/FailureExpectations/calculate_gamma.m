function gamma = calculate_gamma(indexLength, typ)
%CALCULATE_GAMMA: Berechnet beispielhafte Gamma-Werte


gamma_hauptwert = 1;
gamma_varianz = 0.1;

bereich_anfang = gamma_hauptwert - gamma_varianz;
bereich_ende = gamma_hauptwert + gamma_varianz;


gamma_bereich = bereich_anfang :(bereich_ende - bereich_anfang) / (indexLength-1): bereich_ende;

%Nur zum Setzen der Größe, wird später beschrieben
gamma = bereich_anfang :(bereich_ende - bereich_anfang) / (indexLength-1): bereich_ende;
pickerMin = 0;
pickerMax = 1;
picker_softSwitch = 0.25;
iterator = 0.00010;
indexKorrektur = 0;

for i = 1:indexLength
    if typ == "stabil"
        gamma(i) = gamma_bereich(round(0.5 * indexLength));
    elseif typ == "min"
        gamma(i) = gamma_bereich(round(0 * indexLength + 1));
    elseif typ == "max"
        gamma(i) = gamma_bereich(round(1 * indexLength));
    elseif typ == "switchBei25"
        if (mod(i, indexLength * 0.25) < (indexLength * 0.25 / 2))
            gamma(i) = gamma_bereich(round(0 * indexLength + 1));
        else
            gamma(i) = gamma_bereich(round(1 * indexLength));
        end
    elseif typ == "switchBei25_soft"
        if (mod(i, indexLength * 0.5) < (indexLength * 0.5 / 2))
            if picker_softSwitch < 1
                picker_softSwitch = picker_softSwitch + 5* iterator;
                indexKorrektur = 0;
            end
            if picker_softSwitch > 1
                picker_softSwitch = 1;
                indexKorrektur = 0;
            end
            gamma(i) = gamma_bereich(round(picker_softSwitch * indexLength));
        else
            if picker_softSwitch > 0
                picker_softSwitch = picker_softSwitch - 5*iterator;
                indexKorrektur = 1;
            end
            if picker_softSwitch > 1
                picker_softSwitch = 1;
            end
            if picker_softSwitch <= 0
                picker_softSwitch = 0;
                indexKorrektur = 1;
            end
            gamma(i) = gamma_bereich(round(picker_softSwitch * indexLength) + indexKorrektur);
        end
    elseif typ == "steigend"
            if pickerMin < 1
                pickerMin = pickerMin + iterator;
            end
            if pickerMin > 1
                pickerMin = 1;
            end
            gamma(i) = gamma_bereich(round(pickerMin * indexLength));
            
    elseif typ == "sinkend"
            if pickerMax > 0
                pickerMax = pickerMax - iterator;
            end
            if pickerMax < 0
                pickerMax = 0;
            end
            gamma(i) = gamma_bereich(round(pickerMax * indexLength + 1));
    end
end


