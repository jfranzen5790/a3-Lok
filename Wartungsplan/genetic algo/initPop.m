global skalierung n b T numcomp CEL CR CP asize PopSize;

time=0:skalierung:(n-1)*skalierung;
sizetime=size(time);
minT = min(T);
[row, col] = find(T==minT);
minb=b(row,col);

Gip=zeros(1,n);
EWCip=zeros(1,n);
cmin=0; %Hilfsvariable, um Anzahl der Wartungen zu finden
for i = 1:n
    Gip(1,i)=calcG(time(1,i),col);
    EWCip(1,i)=Gip(1,i)*(CEL(1,col)+CR(1,col))/CP(1,col);
    if EWCip(1,i) >= decisionfactor
        cmin=cmin+1;
        for j=1:sizetime(1,2)
            if i+j <= n
                time(1,i+j)=time(1,i+j)-time(1,i);
            end
        end
    end
end

%Matrix mit Initialpopulation erstellen
% IP = randi([1,asize(1,2)],50,n);
IP=ones(PopSize,n)*asize(1,2); 
% IP(:,1)=1; %erster Zeitschritt ist immer Wartung
% IP(:,n)=1; %letzter Zeitschritt ist immer Wartung

IP(1,:)=asize(1,2);
IP(:,1)=1;
% % IP(1,n)=1;

h=round(n/(cmin+1));

for j=1:PopSize
    for i=1:cmin
        
% %         if j == 1
% %             if i == 1
% %                 IP(j,h*i) = 1;
% %             elseif i >= 2
% %                 IP(j,h*i) = 1;
% %             end
% %         elseif j >= 2
            if i == 1
% %                 if j==2
                    % Vektor für Wartungsfolge, in der nur die Komponente mit der
                    %kürzesten Lebensdauer gespeichert wird
                    firstmt=zeros(1,numcomp);
                    for k = 1:numcomp
                        if k == col
                            firstmt(1,k)=1;
                        end
                    end

                    % Finden von firstmt in cell_ci
                    for ii = 1:asize(1,2)
                        if firstmt==cell_ci{1,ii}
                            IP(:,h*i)=119;
                            break
                        end
                    end
% %                 elseif j > 2
% %                     IP(j,h*i)=randi([1 asize(1,2)]);
% %                 end
%                 
%             elseif i >= 2
%                 iprand=rand;
%                 if iprand < 1/3 
%                     IP(j,h*i)=randi([1 asize(1,2)]);
%                 elseif iprand > 1/3 && iprand  <= 2/3
%                     IP(j,h*i)=IP(j,h);
%                 elseif iprand > 2/3
%                     IP(j,h*i)=1;
%                 end
            elseif i >= 1
                iprand=rand;
                if iprand < 1/2 
                    IP(j,h*i)=119;
                elseif iprand > 1/2
                    IP(j,h*i)=1;
                end
            end
            
        end
    end
                

                
                
                
           
                
% %         IP(j,h*i)=1;
%         if i == 1 
%             if j == 1
%                 IP(j,h*i)=1;
%             elseif j >= 2
%                 %Vektor für Wartungsfolge, in der nur die Komponente mit der
%                 %kürzesten Lebensdauer gespeichert wird
%                 firstmt=zeros(1,numcomp);
%                 for k = 1:numcomp
%                     if k == col
%                         firstmt(1,k)=1;
%                     end
%                 end
%             end
%             % Finden von firstmt in cell_ci
%             for ii = 1:asize(1,2)
%                 if firstmt==cell_ci{1,ii}
%                     IP(j,h*i)=ii;
%                     break
%                 end
%             end
%         elseif i >= 2 
%             if j == 1
%                 IP(j,h*i)=1;
%             elseif j >= 2
%                 IP(j,h*i)=randi([1 asize(1,2)]);
%             end
%         elseif i >= 2
%             IP(j,h*i)=IP(j,h*2);
%         end
%     end
% end

%% Fkt. für Best. der Ausfallwahrscheinlichkeit zum Zeitpunkt x
function aw = calcG(x,i)
    global T;
    aw = 1-exp(-x/T(i));
end