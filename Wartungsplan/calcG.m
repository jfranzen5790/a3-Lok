function aw = calcG(x,i)
% Fkt. für Best. der Ausfallwahrscheinlichkeit zum Zeitpunkt x
    global T;
    aw = 1-exp(-x/T(i));
end
