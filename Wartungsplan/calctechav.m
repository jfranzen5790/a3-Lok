function tav=calctechav(sAred)
%techn. Verfügbarkeit bestimmen
global skalierung TGESh;
    dtcycle=sAred*skalierung;
    tav=((TGESh-dtcycle)/TGESh);
end
