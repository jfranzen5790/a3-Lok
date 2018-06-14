function FailureExpectation = calculate_FailureExpectation(t, gammaWert, lagertyp, componentIndex)
    %CALCULATE_F_OF_T: Berechnung der Ausfallwahrscheinlichkeit
    %   Berechnet die Ausfallwahrscheinlichkeit für den gegebenen Zeitpunkt und
    %   gegebene Eigenschaften
    
    lambda = 0.9;
    global b;
    beta = b(componentIndex);
    
    FailureExpectation = 1 - exp(log(lambda)*((t/calculate_lifespan(gammaWert, lagertyp)))^beta);
end

