close all
clear all
clc


Soc0 = 500;
SocP = 800;
SocM = 200;


ts = 1;

a = 100;
b = 600;
pl = a +(b-a)*rand(1, 24);
pv = [zeros(1, 6) 250*ones(1, 4) 600*ones(1, 4) 250*ones(1, 4) zeros(1, 6)]; 



pb = zeros(1, length(pv));
pg =  zeros(1, length(pv));
% pv = [pv pv pv pv pv pv];
% pl = [pl pl pl pl pl pl];
Soc = zeros(1, length(pv));


%% Algo 1

% for i = 1 : ts : length(pv)
%     pdiff = pv(i) - pl(i);
%     if (i==1)
%         if((pdiff > 0 && Soc0 < SocP) || (pdiff <= 0 && Soc0 > SocM))
%             
%              if (pdiff + Soc0) < SocP
%                 pb(i) = pdiff;
%                 pg(i) = 0;
%             else
%                 pb(i) = SocP - Soc0;
%                 pg(i) = pdiff - pb(i);
%             end
%         else
%             pb(i) = 0;
%             pg(i) = pdiff;
%         end
%         Soc(i) = Soc0 + pb(i);
%     else
%          if((pdiff > 0 && Soc(i-1) < SocP) || (pdiff <= 0 && Soc(i-1) > SocM))
%             pb(i) = pdiff;
%             pg(i) = 0;
%         else
%             pb(i) = 0;
%             pg(i) = pdiff;
%         end
%         Soc(i) = Soc(i-1) + pb(i);
%     end
% end

%% Algo2 



% for i = 1 : ts : length(pv)
%     pdiff = pv(i) - pl(i);
%     if (i==1)
%         if (pdiff > 0 && Soc0 < SocP)
%             if (pdiff + Soc0) < SocP
%                 pb(i) = pdiff;
%                 pg(i) = 0;
%             else
%                 pb(i) = SocP - Soc0;
%                 pg(i) = pdiff - pb(i);
%             end
%             Soc(i) = Soc0 + pb(i);
%         end
%         if (pdiff <= 0 && Soc0 > SocM)
%             if (pdiff + Soc0) > SocM
%                 pb(i) = pdiff;
%                 pg(i) = 0;
%             else
%                 pb(i) = SocM - Soc0;
%                 pg(i) = pdiff - pb(i);
%             end
%             Soc(i) = Soc0 + pb(i);
%         end
%             
%     else
%         if (pdiff > 0 && Soc(i-1) < SocP)
%             if (pdiff + Soc(i-1)) < SocP
%                 pb(i) = pdiff;
%                 pg(i) = 0;
%             else
%                 pb(i) = SocP - Soc(i-1);
%                 pg(i) = pdiff - pb(i);
%             end
%             Soc(i) = Soc(i-1) + pb(i);
%         end
%         if (pdiff <= 0 && Soc(i-1) > SocM)
%             if (pdiff + Soc(i-1)) > SocM
%                 pb(i) = pdiff;
%                 pg(i) = 0;
%             else
%                 pb(i) = SocM - Soc(i-1);
%                 pg(i) = pdiff - pb(i);
%             end
%             Soc(i) = Soc(i-1) + pb(i);
%         end
%     end
% end
%  

%% ALgo 3

for i = 1 : ts : length(pv)
    pdiff = pv(i) - pl(i);
    if (i==1)
        if pdiff > 0
            if Soc0 > SocP
                pb(i) = 0;
                pg(i) = pdiff;
            else
                if (pdiff + Soc0 < SocP)
                    pb(i) = pdiff;
                    pg(i) = 0;
                else
                    pb(i) = SocP - Soc0;
                    pg(i) = pdiff - pb(i);
                end
            end
            Ssoc(i) = Soc(i-1) + pb(i);oc(i) = Soc0 + pb(i);
        end
        if pdiff <= 0
            if Soc0 <= SocM
                pb(i) = 0;
                pg(i) = pdiff;
            else
                if (pdiff + Soc0) > SocM
                    pb(i) = pdiff;
                    pg(i) = 0;
                else
                    pb(i) = SocM - Soc0;
                    pg(i) = pdiff - pb(i);
                end
            end
            Soc(i) = Soc0 + pb(i);
        end   
    else
        if pdiff > 0
            if Soc(i-1) > SocP
                pb(i) = 0;
                pg(i) = pdiff;
            else
                if (pdiff + Soc(i-1) < SocP)
                    pb(i) = pdiff;
                    pg(i) = 0;
                else
                    pb(i) = SocP - Soc(i-1);
                    pg(i) = pdiff - pb(i);
                end
            end
            Soc(i) = Soc(i-1) + pb(i);
        end
        if pdiff <= 0
            if Soc(i-1) <= SocM
                pb(i) = 0;
                pg(i) = pdiff;
            else
                if (pdiff + Soc(i-1)) > SocM
                    pb(i) = pdiff;
                    pg(i) = 0;
                else
                    pb(i) = SocM - Soc(i-1);
                    pg(i) = pdiff - pb(i);
                end
            end
            Soc(i) = Soc(i-1) + pb(i);
        end
    end
end
        


figure
title('PV and PL plot')
plot(pv)
hold on 
plot(pl)
legend('PV', 'PL')  
figure
plot(pb)
hold on 
plot(pg)
title('pb and pg plot')
legend('pb', 'pg')

figure
plot(Soc)
legend('Soc')
title('Soc plot')