function l10h = calculate_lifespan(gamma, lagertyp)
    %CALCULATE_LIFESPAN: Berechnung der Lebenserwartung
    % Returns: L10h (= nominelle Lebensdauer in Betriebsstunden)

    %C = dynamische Tragzahl in Newton
    C = 12200;
    %n = lagerdrehzal
    N_lagerdrehzal = 3000;
    
    %n_lagertyp = Lebensdauerexponent: f�r Rollenlager: p_l = 10/3, f�r Kugellager: p_l = 3
    if lagertyp == "rollenlager"
        n_lagertyp = 10/3;
    elseif lagertyp == "kugellager"
        n_lagertyp = 3;
    end
    %Fr = radiale Komponente der tats�chlichen Lagerbelastung in Newton
    Fr = 45 *100;
    %Fa = axiale Komponente der tats�chlichen Lagerbelastung in Newton
    Fa = 5 *100;
    
    %X = Radiallastfaktor nach Tabelle 1 bzw. 2 in Abh�ngigkeit von Fa / Fr
    x = 0.56;
    
    %y = Axiallastfaktor nach Tabelle 1 bzw. 2 in Abh�ngigkeit von Fa / Fr
    y = 1.3;
    

    %P = dynamische �quivalente >Lagerlast< in Newton
    P = (x * Fr + y * Fa) * gamma;
    
    %L10 = nominelle Lebensdauer f�r eine >Erlebenswahrscheinlichkeit< von 90 % in 106 �berrollungen
    L10 = (C/P)^n_lagertyp;

    
    %L10h = nominelle Lebensdauer in Betriebsstunden
    l10h = (10^6/(60 * N_lagerdrehzal)) * L10;
end

