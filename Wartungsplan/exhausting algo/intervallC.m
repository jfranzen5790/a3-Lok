Cmplanstep=zeros(m,n,m);
for bb=1:m
    for cc=1:m
        for dd=1:n % Jede Zelle kombinieren
            if mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 0
                Cmplanstep(cc,dd,bb)=CIS(1)+ MTTR(1)*M + MTTRD(1)*CELd;
            elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 1
                Cmplanstep(cc,dd,bb)=CIS(2)+ MTTR(2)*M + MTTRD(2)*CELd;
            elseif mplan(bb,dd,1) == 1 && mplan(cc,dd,2) == 1
                Cmplanstep(cc,dd,bb)=CIS(1) + CIS(2) + (MTTR(1)+MTTR(2))*M + max(MTTRD)*CELd;
            elseif mplan(bb,dd,1) == 0 && mplan(cc,dd,2) == 0
                Cmplanstep(cc,dd,bb)=0;
            elseif mplan(bb,dd,1) <= 1 && mplan(cc,dd,2) == inf
                Cmplanstep(cc,dd,bb)=inf;
            elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= 1
                Cmplanstep(cc,dd,bb)=inf;
            elseif mplan(bb,dd,1) == inf && mplan(cc,dd,2) <= inf
                Cmplanstep(cc,dd,bb)=inf;
            end
        end

        CGES(bb,cc)=sum(Cmplanstep(cc,:,bb));
    end
end
