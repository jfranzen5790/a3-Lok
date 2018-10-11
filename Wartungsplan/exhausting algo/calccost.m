global PopSize skalierung numcomp asize cell_ci CGES n Cmax arbeitszeit transmat G EWC A TGESh allCosts costList gueltigCounter ungueltigCounter bestCost bestComb bestInput;

% Wartungsschritte und -kosten pro Zeitschritt
%combmplan=cell(max(q10)^numcomp,n);
%combCGES=zeros(max(q10)^numcomp,n);
%CGES=zeros(max(q10)^numcomp,1);
i4=1;
currentInput = zeros(0);
bestInput = zeros(0);
bestCost = 1;
currentCost = 1;

bestComb = zeros(1, n);
bestComb(:) = 2^numcomp;
currentComb = zeros(1, n);
currentComb(:) = 2^numcomp;

if numcomp==3
    for bb=1:length(mplanred)
        for cc=1:length(mplanred)
            for dd=1:length(mplanred)
                for xx=1:n
                    vector=[mplanred(bb,xx,1) mplanred(cc,xx,2) mplanred(dd,xx,3)];
                    %combmplan{i4,xx}=vector;
                    currentInput{xx}=vector;
                end
                    currentInputIndexVersion = bi2de(cell2mat(currentInput(:))) + 1;
                    currentCost = (cost(currentInputIndexVersion));
                    i4=i4+1;
                    if currentCost < bestCost
                        bestCost = currentCost;
                        bestComb = currentInputIndexVersion;
                        bestInput = currentInput;
                    end
            end    
        end
    end
elseif numcomp==4
    for bb=1:size(mplanred, 1)
        for cc=1:size(mplanred, 1)
            for dd=1:size(mplanred, 1)
                for ee=1:size(mplanred, 1)
                    for xx=1:n
                        %currentInput
                        vector=[mplanred(bb,xx,1) mplanred(cc,xx,2) mplanred(dd,xx,3) mplanred(ee,xx,4)];
                        %combmplan{i4,xx}=vector;
                        currentInput{xx}=vector;
                        %cell_cisize=size(cell_ci);
                        %for i6=1:cell_cisize(1,2)
                            %if isequal(combmplan{i4,ee},cell_ci{1,i6})
                            %    combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
                            %end
                        %end
                    end
                    % Kosten der Wartungsstrategie = Summe der Kosten für
                    % Zeitschritte
                    %CGES(i4,1)=sum(combCGES(i4,:));
                    currentInputIndexVersion = bi2de(cell2mat(currentInput(:))) + 1;
                    currentCost = (cost(currentInputIndexVersion));
                    %[trans, transmat] = calctransmat(currentInputIndexVersion);
                    %C = calcCost(trans);
                    %disp(C);
                    i4=i4+1;
                    if currentCost < bestCost
                        bestCost = currentCost;
                        bestComb = currentInputIndexVersion;
                        bestInput = currentInput;
                    end
                end
            end    
        end
    end

elseif numcomp==5
    for bb=1:size(mplanred, 1)
        for cc=1:size(mplanred, 1)
            for dd=1:size(mplanred, 1)
                for ee=1:size(mplanred, 1)
                    for ff=1:size(mplanred, 1)
                        for xx=1:n
                            %currentInput
                            vector=[mplanred(bb,xx,1) mplanred(cc,xx,2) mplanred(dd,xx,3) mplanred(ee,xx,4) mplanred(ff,xx,5)];
                            %combmplan{i4,xx}=vector;
                            currentInput{xx}=vector;
                            %cell_cisize=size(cell_ci);
                            %for i6=1:cell_cisize(1,2)
                                %if isequal(combmplan{i4,ee},cell_ci{1,i6})
                                %    combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
                                %end
                            %end
                        end
                        % Kosten der Wartungsstrategie = Summe der Kosten für
                        % Zeitschritte
                        %CGES(i4,1)=sum(combCGES(i4,:));
                        currentInputIndexVersion = bi2de(cell2mat(currentInput(:))) + 1;
                        
                        currentCost = (cost(currentInputIndexVersion));
                        %[trans, transmat] = calctransmat(currentInputIndexVersion);
                        %C = calcCost(trans);
                        %disp(C);
                        i4=i4+1;
                        if currentCost < bestCost
                            bestCost = currentCost;
                            bestComb = currentInputIndexVersion;
                            bestInput = currentInput;
                        end
                    end
                end
            end    
        end
    end
elseif numcomp==6
    for bb=1:size(mplanred, 1)
        for cc=1:size(mplanred, 1)
            for dd=1:size(mplanred, 1)
                for ee=1:size(mplanred, 1)
                    for ff=1:size(mplanred, 1)
                        for gg=1:size(mplanred, 1)
                            for xx=1:n
                                %currentInput
                                vector=[mplanred(bb,xx,1) mplanred(cc,xx,2) mplanred(dd,xx,3) mplanred(ee,xx,4) mplanred(ff,xx,5), mplanred(gg,xx,6)];
                                %combmplan{i4,xx}=vector;
                                currentInput{xx}=vector;
                                %cell_cisize=size(cell_ci);
                                %for i6=1:cell_cisize(1,2)
                                    %if isequal(combmplan{i4,ee},cell_ci{1,i6})
                                    %    combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
                                    %end
                                %end
                            end
                            % Kosten der Wartungsstrategie = Summe der Kosten für
                            % Zeitschritte
                            %CGES(i4,1)=sum(combCGES(i4,:));
                            currentInputIndexVersion = bi2de(cell2mat(currentInput(:))) + 1;

                            currentCost = (cost(currentInputIndexVersion));
                            %[trans, transmat] = calctransmat(currentInputIndexVersion);
                            %C = calcCost(trans);
                            %disp(C);
                            i4=i4+1;
                            if currentCost < bestCost
                                bestCost = currentCost;
                                bestComb = currentInputIndexVersion;
                                bestInput = currentInput;
                            end
                        end
                    end
                end
            end    
        end
    end
elseif numcomp==7
    for bb=1:size(mplanred, 1)
        for cc=1:size(mplanred, 1)
            for dd=1:size(mplanred, 1)
                for ee=1:size(mplanred, 1)
                    for ff=1:size(mplanred, 1)
                        for gg=1:size(mplanred, 1)
                            for hh=1:size(mplanred, 1)
                                for xx=1:n
                                    %currentInput
                                    vector=[mplanred(bb,xx,1) mplanred(cc,xx,2) mplanred(dd,xx,3) mplanred(ee,xx,4) mplanred(ff,xx,5), mplanred(gg,xx,6) mplanred(hh,xx,7)];
                                    %combmplan{i4,xx}=vector;
                                    currentInput{xx}=vector;
                                    %cell_cisize=size(cell_ci);
                                    %for i6=1:cell_cisize(1,2)
                                        %if isequal(combmplan{i4,ee},cell_ci{1,i6})
                                        %    combCGES(i4,ee)=sum(combCP(:,i6))+combCELmax(1,i6);
                                        %end
                                    %end
                                end
                                % Kosten der Wartungsstrategie = Summe der Kosten für
                                % Zeitschritte
                                %CGES(i4,1)=sum(combCGES(i4,:));
                                currentInputIndexVersion = bi2de(cell2mat(currentInput(:))) + 1;

                                currentCost = (cost(currentInputIndexVersion));
                                %[trans, transmat] = calctransmat(currentInputIndexVersion);
                                %C = calcCost(trans);
                                %disp(C);
                                i4=i4+1;
                                if currentCost < bestCost
                                    bestCost = currentCost;
                                    bestComb = currentInputIndexVersion;
                                    bestInput = currentInput;
                                end
                            end
                        end
                    end
                end
            end    
        end
    end
end

