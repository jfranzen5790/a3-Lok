function lt = calclifetime(transmat2, Aredmat)
    global n numcomp skalierung
    lt=zeros(numcomp,n);
    for i = 1:numcomp
        wartung=ones(1,numcomp);
        jlt=1;
        skipi=0;
        while jlt <= n
          
            if jlt==1
                if transmat2(i,jlt) == 0
                    lt(i,jlt)=jlt*skalierung;
                elseif transmat2(i,jlt) == 1
                    lt(i,jlt)=0;
                    wartung(1,i)=jlt;
                end
            elseif jlt > 1
                if sum(transmat2(:,jlt)) == 0
                    if transmat2(i,jlt-skipi) == 0
                        lt(i,jlt)=lt(i,wartung(1,i))+(jlt-wartung(1,i))*skalierung;
                    elseif transmat2(i,jlt-skipi) == 1
                        lt(i,jlt)=(jlt-wartung(1,i))*skalierung;
                    end
                    skipi=0;
        
                elseif sum(transmat2(:,jlt)) > 0
                    skipi=Aredmat(2,jlt);
                    wartung(1,i)=jlt-1+skipi;
                    for k=1:Aredmat(2,jlt)
                        if jlt-1+k<n
                            if transmat2(i,jlt) == 1
                                lt(i,jlt-1+k)=0;
                            elseif transmat2(i,jlt) == 0
                                lt(i,jlt-1+k)=lt(i,jlt-1);
                            end
                        end
                    end
                end
            end
            if skipi == 0
                jlt=jlt+1;
            elseif skipi > 0
                jlt=jlt+skipi;
            end
        end
        jlt=1;
    end
end
