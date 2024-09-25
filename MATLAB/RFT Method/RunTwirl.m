%% Run code for torsion and analysis

clear
cla
% numcore = 12
% delete(gcp('nocreate'));
% parpool('local',numcore);
DataFN = 'TwirlWhirlData_2021-10-26_Ribbon_RampUp_Hysteresis';        % File name for multi run data
%%
row = 0;
for q = 1:1
    for w = 1:1
        for run = 1:1
            for f = 1:1
                row = row+1;
                [q,w,run,f]
                cla
                close all
                Gamma = 1;  % ratio of A/C
                sig0 = 1E4;  % sig0 is good at 1E5/Gridnum (when L = 1)
                GridNum = 100;
                sig0 = sig0./GridNum;
                a = 0.1;
                etalist = [1/(4*pi)];
                eta = etalist(1);
%                 AmodList = [0.0,0.2,0.4,0.6,0.8,0.90];
%                 A1 = 1.0 + AmodList(w);     
%                 A2 = 1.0 - AmodList(w);      
                dtlist = [4E-6,1E-8];      % small time step.... units of seconds per time step
                dt = dtlist(1);
                time1 = 0.0001;
                time2 = 1.0;
                zetar = 4*pi*eta*a*a;
%                 Ap = 0.5.*(A1+A2);
%                 om0list = [12.0,9.0].*Ap./(zetar); 
%                 om0 = om0list(1);
%                 om01 = 12.0.*Ap./(zetar);
%                 om02 = (9.0- (f-1).*0.2).*Ap./(zetar);
                TotalTime = time1+time2;              % total time of the simulation (seconds)
                NSteps = 1000;
                THF = 0 ;                % turns thermal force on(1) or off(0)
                RNGmod = 'default';   % either default or shuffle in order to initialize the random number generator. specifying a number sets 'seed' of Mersenne Twister RNG
                
             %% Background 1
%                 A1 = 1.0 + AmodList(2);     
%                 A2 = 1.0 - AmodList(2);  
%                 Ap = 0.5.*(A1+A2);
%                 om01 = 8.8.*Ap./(zetar)
%                 om02 = 8.8.*Ap./(zetar);
%                 myfunc1 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps);
%             %% Background 2
%                 A1 = 1.0 + AmodList(2);    
%                 A2 = 1.0 - AmodList(2);  
%                 Ap = 0.5.*(A1+A2);
%                 om01 = 8.78.*Ap./(zetar)
%                 om02 = 8.78.*Ap./(zetar);
%                 myfunc2 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps);
%             %% Background 3
%                 A1 = 1.0 + AmodList(2);    
%                 A2 = 1.0 - AmodList(2);  
%                 Ap = 0.5.*(A1+A2);
%                 om01 = 8.82.*Ap./(zetar)
%                 om02 = 8.82.*Ap./(zetar);
%                 myfunc3 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps);
%             %% Background 4
%                 A1 = 1.0 + AmodList(2);    
%                 A2 = 1.0 - AmodList(2);  
%                 Ap = 0.5.*(A1+A2);
%                 om01 = 8.85.*Ap./(zetar)
%                 om02 = 8.85.*Ap./(zetar);
%                 myfunc4 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps); 
%             %% Background 5
%                 A1 = 1.0 + AmodList(2);    
%                 A2 = 1.0 - AmodList(2);  
%                 Ap = 0.5.*(A1+A2);
%                 om01 = 8.88.*Ap./(zetar)
%                 om02 = 8.88.*Ap./(zetar);
%                 myfunc5 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps);
%                     
            %% Foreground  
                AmodList = [.9];
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = 4.*Ap./(zetar)
                om02 = 15.*Ap./(zetar);
                [x,y,z,e1,Time,mov,vidtitle,GridNum,info]=TwirlWhirl(Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps);
            
            %     [x,y,z,e1,Time,mov,vidtitle,GridNum,info]=TwirlWhirl(Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps);
            %     F = parfeval(backgroundPool,fcn,n,X1,...,Xm) schedules the function fcn to run in the background. 
                % parfeval runs the function fcn on a background worker. 
                % MATLAB evaluates the function fcn [Y1,...,Yn] = fcn(X1,...,Xm), with m inputs and n outputs.
%             wait(myfunc1)
%             wait(myfunc2)
%             wait(myfunc3)
%             wait(myfunc4)
%             wait(myfunc5)

                
% "done with all"
                
            end 
        end
    end
end

% AggTable3 = table(IterData2(:,1),IterData2(:,2),IterData2(:,3),IterData2(:,4),IterData2(:,5),IterData2(:,6),IterData2(:,7),...
%                   IterData2(:,8),IterData2(:,9),IterData2(:,10),IterData2(:,11),IterData2(:,12),IterData2(:,13),...
%                   IterData2(:,14),IterData2(:,15),IterData2(:,16),IterData2(:,17),...
%                   'VariableNames',{'eta','GammaX','GammY','GammaZ', 'OrigLength','radius','Ac1','Ac2','Cc','Version', ...
%                   'Sigma','GN','dt','GoldSigma','1-min(t2tL/len)','FinalLength','run'});
% writetable(AggTable3,strcat(DataFN,'_data.xlsx'),'Sheet','FlowData');
