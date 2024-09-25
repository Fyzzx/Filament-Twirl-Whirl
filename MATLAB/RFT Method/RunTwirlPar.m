%% Run code for torsion and analysis

clear
cla
numcore = 12
delete(gcp('nocreate'));
parpool('local',numcore);
DataFN = 'TwirlWhirlData_2022-03-29_Visc1_1stOm_Sig3rdO5E5_dt2E-8_VaryA_plectStop';        % File name for multi run data
IterData2 = [];
%%
row = 0;
for q = 1:1
    for w = 1:1
        for run = 1:1
            for f = 1:6
%                 row = row+1;
                [q,w,run,f]
                cla
                close all
                Gamma = 1;  % ratio of A/C
                sig0 = 5E4;  % sig0 is good at 1E5/Gridnum (when L = 1)
                GridNum = 101;
                a = 0.05;
                etalist = [1/(317*pi)];
                eta = etalist(1);
                AmodList = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
%                 A1 = 1.0 + AmodList(w);     
%                 A2 = 1.0 - AmodList(w);      
%                 dtlist = [2E-7,2E-7,2E-7,2E-7,2E-7,2E-7,2E-7,2E-7,2E-7,2E-7];      % small time step.... units of seconds per time step
%                 dtlist = [5E-7,5E-7,5E-7,5E-7,5E-7,5E-7,5E-7,5E-7,5E-7,5E-7];      % small time step.... units of seconds per time step
%                 dtlist = [1E-6,1E-6,1E-6,1E-6,1E-6,1E-6,1E-6,1E-6,1E-6,1E-6];      % small time step.... units of seconds per time step
%                 dtlist = [1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7,1E-7];      % small time step.... units of seconds per time step
%                 dtlist = [1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8,1E-8];
                dtlist = [2E-8,2E-8,2E-8,2E-8,2E-8,2E-8,2E-8,2E-8,2E-8,2E-8];      % small time step.... units of seconds per time step

                dt = dtlist(f);
%                 timelist = [1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5]; % for first omList list
%                 timelist = [0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6];
%                 timelist = [0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05];
%                 timelist = [0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3];
%                 timelist = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];
                timelist = 15*[0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005];
                
                time1 = timelist(f);
                time2 = 0.0000;
                zetar = 4*pi*eta*a*a;
                TotalTime = time1+time2;              % total time of the simulation (seconds)
                NSteps = 500;
                THF = 0 ;                % turns thermal force on(1) or off(0)
                RNGmod = 'default';   % either default or shuffle in order to initialize the random number generator. specifying a number sets 'seed' of Mersenne Twister RNG
                omList = [9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0];         % longest times for 9.3 = 1.5s for cylinder, 1.15 for 0.4 ribbon
%                 omList = [10.2, 10.4, 10.6, 10.8, 11.0, 11.2, 11.4, 11.6];  % longest times for 10.2 = 0.45s
%                 omList = [11.8, 12.0, 12.2, 12.4, 12.6, 12.8, 13.0, 13.2];  % longest times for 11.8 = 0.20s
%                 omList = [13.4, 13.6, 13.8, 14.0, 14.2, 14.4, 14.6, 14.8];  % longest times for 13.4 = 0.10s
            %% Background 1
                A1 = 1.0 + AmodList(f);     
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(1).*Ap./(zetar)
                om02 = omList(1).*Ap./(zetar);
                sig0 = (Ap).*sig0./(GridNum-1);
                myfunc1 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);
            
            %% Background 2
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(2).*Ap./(zetar)
                om02 = omList(2).*Ap./(zetar);
                myfunc2 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);
            
            %% Background 3
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(3).*Ap./(zetar)
                om02 = omList(3).*Ap./(zetar);
                myfunc3 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);
            
            %% Background 4
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(4).*Ap./(zetar)
                om02 = omList(4).*Ap./(zetar);
                myfunc4 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);                    
            
            %% Background 5
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(5).*Ap./(zetar)
                om02 = omList(5).*Ap./(zetar);
                myfunc5 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);    
           
            %% Background 6
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(6).*Ap./(zetar)
                om02 = omList(6).*Ap./(zetar);
                myfunc6 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);  
           
            %% Background 7
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(7).*Ap./(zetar)
                om02 = omList(7).*Ap./(zetar);
                myfunc7 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);  
          
            %% Background 8
                A1 = 1.0 + AmodList(f);    
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(8).*Ap./(zetar)
                om02 = omList(8).*Ap./(zetar);
                myfunc8 = parfeval(@TwirlWhirl,9,Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);  
           
            %% Foreground  
                A1 = 1.0 + AmodList(f);     
                A2 = 1.0 - AmodList(f);  
                Ap = 0.5.*(A1+A2);
                om01 = omList(8).*Ap./(zetar)
                om02 = omList(8).*Ap./(zetar);
%                 sig0 = (Ap).*sig0./(GridNum-1);  %%% REPEAT %%%
                [x,y,z,e1,Time,mov,vidtitle,GridNum,info]=TwirlWhirl(Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps);
            
        %     [x,y,z,e1,Time,mov,vidtitle,GridNum,info]=TwirlWhirl(Gamma,sig0,a,GridNum,dt,A1,A2,eta,TotalTime,THF,RNGmod,time1,time2,om01,om02,NSteps);
        %     F = parfeval(backgroundPool,fcn,n,X1,...,Xm) schedules the function fcn to run in the background. 
                % parfeval runs the function fcn on a background worker. 
                % MATLAB evaluates the function fcn [Y1,...,Yn] = fcn(X1,...,Xm), with m inputs and n outputs.
           
            wait(myfunc1)
            myinfonow = cell2mat(myfunc1.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);
            
            wait(myfunc2)
            myinfonow = cell2mat(myfunc2.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);

            wait(myfunc3)
            myinfonow = cell2mat(myfunc3.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);

            wait(myfunc4)
            myinfonow = cell2mat(myfunc4.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);
% 
            wait(myfunc5)
            myinfonow = cell2mat(myfunc5.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);

            wait(myfunc6)
            myinfonow = cell2mat(myfunc6.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);

            wait(myfunc7)
            myinfonow = cell2mat(myfunc7.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);

            wait(myfunc8)
            myinfonow = cell2mat(myfunc8.OutputArguments(9));
            [row, IterData2] = logData(row,IterData2,myinfonow);

                
"done with all"
                
            end 
        end
    end
end

AggTable3 = table(IterData2(:,1),IterData2(:,2),IterData2(:,3),IterData2(:,4),IterData2(:,5),...
                  IterData2(:,6),IterData2(:,7),IterData2(:,8),IterData2(:,9),IterData2(:,10),...
                  IterData2(:,11),IterData2(:,12),IterData2(:,13),IterData2(:,14),IterData2(:,15),...
                  IterData2(:,16),IterData2(:,17),IterData2(:,18),IterData2(:,19),IterData2(:,20),...
                  IterData2(:,21),IterData2(:,22),IterData2(:,23),IterData2(:,24),IterData2(:,25),...
                  IterData2(:,26),IterData2(:,27),IterData2(:,28),IterData2(:,29),IterData2(:,30),...
                  IterData2(:,31),IterData2(:,32),IterData2(:,33),IterData2(:,34),IterData2(:,35),...
                  IterData2(:,36),IterData2(:,37),IterData2(:,38),IterData2(:,39),...
                  'VariableNames',{'omRatio','Am/Ap','xEndPt','XEndTwirlTime','TimeE3Flip',...
                                   'yZeros', 'zZeros','MaxMaxStretch','AvgMaxStretch','MaxEner',...
                                   'FinalAvgEner','EndEnergy','ETrend','maxTwistEner','maxBendEner',...
                                   'endTwEner','endBendEner','omC','om01','time1',...
                                   'om02','time2','Length','radius','eta',...
                                   'zetar','Gamma','A1','A2','sig0',...
                                   'GridNum','dt','THF','TotalTime','bailout',...
                                   'Plect', 'PlectLoc1', 'PlectLoc2','MostRightNode'});


writetable(AggTable3,strcat(DataFN,'_data.xlsx'),'Sheet','TwirlWhirl');

            


%% function to write info data to IterData
function [row,IterData2] = logData(row,IterData2,myinfonow)
    row = row+1;
    IterData2(row,:) = [myinfonow];
end

