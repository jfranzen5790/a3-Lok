function tav=calctechav(sAred)
%techn. Verf�gbarkeit bestimmen
global skalierung TGESh;
    dtcycle=sAred*skalierung;
    tav=((TGESh-dtcycle)/TGESh);
end
