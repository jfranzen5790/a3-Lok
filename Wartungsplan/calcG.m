function aw = calcG(x,i)
% Fkt. f�r Best. der Ausfallwahrscheinlichkeit zum Zeitpunkt x
    global T;
    aw = 1-exp(-x/T(i));
end
