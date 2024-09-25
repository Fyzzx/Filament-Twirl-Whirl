function [x,y,z,e1,Time,mov,vidtitle,GridNum,infodata]=TwirlWhirl(Gamma,sig0,a,GridNum,dt,A1,A2,eta,THF,RNGmod,time1,time2,om01,om02,NSteps)

%% Define parameters
% % Parallel parameters
% numcore =4
% delete(gcp('nocreate'));
% parpool('local',numcore);

%Time parameters
time1 = 0.005;
time2 = 0;
dt = 2*10^-8;  %time step
NSteps = 500;

% Filamaent Parameters
Gamma = 1.0;   % ratio of the bend modulus to the twist modulus
L = 1.0;
GridNum = 101;
s = linspace(0,L,GridNum);
ds = s(2) - s(1);
sig0 = 5E4;
sig0 = sig0.*ds;
a = 0.05; % radius
A1 = 1.5; % bend modulus about e1 direction
A2 = 0.5;


% Fluid parameters
eta = 1/(317*pi);     % 0.001 is water viscosity
THF = 0;         % thermal factor that will dictates whether thre are or are not forces & moments from thermal
RNGmod = 'default';

% spin parameters
om01 = 11.6./eta./(4*pi)./(a.^2)./(L.^2);
om02 = om01;

% % Plect Paraemters
% cRep = 0.01; % constant of repulsion f = (cRep)/delR^2
cRep = 0; % shuts off the attempt to find this
repNum = 0;
plect = 0;
plectLoc = [0,0];
delRMagperp = 1E10.*ones(GridNum,1);
plotbox = 1;
%%
% set(0, 'DefaultFigureRenderer', 'OpenGL');
te3flip = 9999;
L=1;
if om01 == om02
    time2 = 0.0;
end
time3=time2;
time4=time2;
TotalTime = time1+time2+time3+time4;

omList = zeros(NSteps-1,1);
Skip = ceil(TotalTime/dt/NSteps); 
tlist = linspace(TotalTime./NSteps,TotalTime,NSteps-1);
tlist = tlist';

Ap = 0.5.*(A1+A2);
Am = 0.5.*(A1-A2);
Arat = Am/Ap;
s = linspace(0,L,GridNum);
ds = s(2) - s(1);
closenodes = round(2.3*a./ds,0)+1;

zetar = 4*pi*eta.*a.^2;     
zetaperp = 4*pi*eta;      

dragrat = 0.5;      % ratio of the parallel to perpendicular drag coefficents

k1 = 0;         % preferred curvature along the e1 direction
k2 = 0;         % preferred curvature along the e2 direction
tc0 = 0;         % preferred torsion for cell body
rad = a.*(ones(GridNum,1));
t0Half = tc0.*ones(GridNum-1,1);
E1Ang = 0.0;
bailout = 0;
%% Preallocated for speed
E1dat(:,:,:)= zeros(GridNum,NSteps,3);
E2dat(:,:,:)= zeros(GridNum,NSteps,3);
ElEner1 = zeros(NSteps-1,1); ElEner2 = zeros(NSteps-1,1);
ElEner3 = zeros(NSteps-1,1); ElEnerStr = zeros(NSteps-1,1);
Etot = zeros(NSteps-1,1);
%%
dtstr = datestr(now,'yy_mm_dd_HH_MM');
omFactor = (Ap/(L.^2 * zetar));
omRat = om01./omFactor;
%% Thermal
kT = 0.004114;  % room temp kT in pN*um
TFCM = sparse([1 GridNum],[1 GridNum],1./sqrt(2),GridNum,GridNum) ...
    + sparse((2:GridNum-1),(2:GridNum-1),1,GridNum,GridNum);
rng(RNGmod) ; % 'default' with each new run of the code sets RNG to SAME set of random numbers. Each j,k loop iteration get different numbers. 'shuffle' sets the rng based on current time  
mTH = zeros(GridNum,1);

if THF ==1
    TFe1 = sqrt(2*kT*zetaperp*dragrat*ds/dt);
    TFe3 = sqrt(2*kT*zetaperp*ds/dt); % units of force
    TFRot = sqrt(2*kT*zetar/ds/dt); % units of torque/length

else 
    TFe1 = 0.0;
    TFe3 = 0.0;
    TFRot = 0.0;
end
%%  create derivative matrices

L1 = FirstDer6(GridNum,ds);
L2 = SecondDer6(GridNum,ds);
 
O1 =    sparse((3:GridNum-2),(2:GridNum-3),1,GridNum,GridNum) ...
      + sparse((3:GridNum-2),(3:GridNum-2),-2,GridNum,GridNum) ...
      + sparse((3:GridNum-2),(4:GridNum-1),1,GridNum,GridNum) ...
      + sparse(1,2,1,GridNum,GridNum) ...
      + sparse(2,(2:3),[-2 1],GridNum,GridNum) ...
      + sparse(GridNum-1,(GridNum-2:GridNum),[1 -2 0],GridNum,GridNum) ...
      + sparse(GridNum,GridNum-1,1,GridNum,GridNum);
  
O1 = O1./ds; % 2nd spatial derivative

O2 =    sparse((3:GridNum-1),(2:GridNum-2),-0.5,GridNum,GridNum) ...
      + sparse((2:GridNum-2),(3:GridNum-1),0.5,GridNum,GridNum) ...
      + sparse(1,2,0.5,GridNum,GridNum) ...
      + sparse(GridNum,GridNum-1,-0.5,GridNum,GridNum);
  
O3 =   sparse((1:GridNum-1),(1:GridNum-1),1,GridNum,GridNum-1) ...
     + sparse((2:GridNum),(1:GridNum-1),-1,GridNum,GridNum-1);
  
%% Create some matrices

Zero = sparse(GridNum,GridNum);

CV =  sparse([1 GridNum],[1 GridNum],ds./2,GridNum,GridNum) ...
    + sparse((2:GridNum-1),(2:GridNum-1),ds,GridNum,GridNum);

CV2 = [ CV     Zero   Zero;
        Zero   CV     Zero;
        Zero   Zero   CV ];


% define the initial shape of the filament

Om3 = zeros(GridNum,1);

k = 1.8751; % 

x = s';
y = zeros(GridNum,1);
z = 0.02.*( sin(k*s') - sinh(k*s') + ( cos(k) + cosh(k) ).*( cos(k*s') - cosh(k*s') )./( sin(k) + sinh(k) ) );

e1(:,1) = zeros(GridNum,1);
e1(:,2) = ones(GridNum,1);
e1(:,3) = zeros(GridNum,1);

mag = sqrt(e1(:,1).^2 + e1(:,2).^2 + e1(:,3).^2);

e1(:,1) = e1(:,1)./mag;
e1(:,2) = e1(:,2)./mag;
e1(:,3) = e1(:,3)./mag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time(1) = 0;

%%
for j=2:NSteps
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Iter = j;
    x(:,j) = x(:,j-1);
    y(:,j) = y(:,j-1);
    z(:,j) = z(:,j-1);
    Time(j) = Time(j-1);

if Time(j)<= time1
    om0 = om01 + ((Time(j))./time1).*(om02-om01)./2;
else 
    if Time(j)<= time1 + time2
        om0 = om01 + (om02-om01)./2; 
    else
        if Time(j)<= time1 + time2 + time3
            om0 = om01 + (om02-om01)./2 + ((Time(j)-time1-time2)./time3).*(om02-om01)./2;
        else
            om0 = om02;
        end

    end
      
end
omList(j-1,1) = om0;

%%%% Take Full Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if bailout == 0
    for k = 1:Skip
%%        

        Time(j) = Time(j) + dt;

        %% define derivatives
        
        xs = L1*x(:,j);
        ys = L1*y(:,j);
        zs = L1*z(:,j);
        x2s = L2*x(:,j);
        y2s = L2*y(:,j);
        z2s = L2*z(:,j);
        
%% define unit vectors

    mag = sqrt( xs.^2 + ys.^2 + zs.^2 );

    e3(:,1) = xs./mag;
    e3(:,2) = ys./mag;
    e3(:,3) = zs./mag;
    
    e2(:,1) = e3(:,2).*e1(:,3) - e3(:,3).*e1(:,2);
    e2(:,2) = e3(:,3).*e1(:,1) - e3(:,1).*e1(:,3);
    e2(:,3) = e3(:,1).*e1(:,2) - e3(:,2).*e1(:,1);
    
    mag2 = sqrt( e2(:,1).^2 + e2(:,2).^2 + e2(:,3).^2 );
    
    e2(:,1) = e2(:,1)./mag2;
    e2(:,2) = e2(:,2)./mag2;
    e2(:,3) = e2(:,3)./mag2;
    
    Om1(1) = k1(1);
    Om2(1) = k2(1);
    Om1(end) = k1(end);
    Om2(end) = k2(end);

    Ftx = zeros(GridNum,1);
    Fty = zeros(GridNum,1);
    Ftz = zeros(GridNum,1);
    if THF == 1
        randvar = normrnd(0,1, [GridNum*3,1]);
        Fte1 = randvar(1:GridNum,1).*TFe1; % thermal force causing cylinder to translate in e1
        Fte2 = randvar(GridNum+1:2*GridNum,1).*TFe1;
        Fte3 = randvar(2*GridNum+1:end,1).*TFe3; % thermal force causing cylinder to translate in e3

        Fte1 = TFCM*Fte1; % modifying edges by factor of sqrt(2)
        Fte2 = TFCM*Fte2;
        Fte3 = TFCM*Fte3;

        Ftx = Fte1.*e1(:,1) + Fte2.*e2(:,1) + Fte3.*e3(:,1);
        Fty = Fte1.*e1(:,2) + Fte2.*e2(:,2) + Fte3.*e3(:,2);
        Ftz = Fte1.*e1(:,3) + Fte2.*e2(:,3) + Fte3.*e3(:,3); 
    end 
    
%% define curvatures

   Om1 = -( e2(:,1).*x2s + e2(:,2).*y2s + e2(:,3).*z2s );
   Om2 = ( e1(:,1).*x2s + e1(:,2).*y2s + e1(:,3).*z2s );    %Original
   
   MOm = A1.*Om1.*(Om1 - k1) + A2.*Om2.*(Om2 - k2);
   MOm2 = -A1.*Om2.*(Om1 - k1) + A2.*Om1.*(Om2 - k2);
   
%% define twist and matrix for updating e1
% function [SourceE1, MatE1, Om3]=MakeE1Mats(GridNum, e2, ds, e1, zetar, om0, L, Gamma, Ap, CV, dt, MOm2, O3, t0Half, CV2, mTH)
% inputs = GridNum, e2, ds, e1, T0R, T0L, zetar, om0, L, Gamma, Ap, CV, dt,
% MOm2, O3, t0Half, CV2
% outputs SourceE1, MatE1, Om3
% if THF == 1
%     mTH = TFRot.*randn(GridNum,1);
% end

%     E2xMat = sparse((1:GridNum),(1:GridNum),e2(:,1),GridNum,GridNum);
%     E2yMat = sparse((1:GridNum),(1:GridNum),e2(:,2),GridNum,GridNum);
%     E2zMat = sparse((1:GridNum),(1:GridNum),e2(:,3),GridNum,GridNum);
%     
%     TwX = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,1)/ds,GridNum-1,GridNum) ...
%           + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,1)/ds,GridNum-1,GridNum);
%     TwY = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,2)/ds,GridNum-1,GridNum) ...
%           + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,2)/ds,GridNum-1,GridNum);
%     TwZ = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,3)/ds,GridNum-1,GridNum) ...
%           + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,3)/ds,GridNum-1,GridNum);
% 
% Tw = [ TwX, TwY, TwZ ];
% e1Vec = [ e1(:,1); e1(:,2); e1(:,3) ];
%       
% Om3Half = Tw*e1Vec;
% Om3(1) = -zetar.*om0.*L./Gamma./Ap;
% Om3(2:GridNum-1) = 0.5.*(Om3Half(1:GridNum-2) + Om3Half(2:GridNum-1));
% Om3(GridNum) = tc0;
% 
% TwMat = [ E2xMat*O3;
%           E2yMat*O3;
%           E2zMat*O3 ];
% 
% MatE1 = CV2 - (Gamma.*dt./zetar).*(TwMat*Tw);
%       
% SourceE1 = CV2*e1Vec - (Ap.*L.*Gamma.*dt./zetar).*(TwMat*t0Half);
% 
% SourceE1 = SourceE1 + (dt./zetar).*[ CV*((MOm2+mTH).*e2(:,1)); (CV*(MOm2+mTH).*e2(:,2)); (CV*(MOm2+mTH).*e2(:,3)) ];
% end
[SourceE1, MatE1, Om3]=MakeE1Mats(GridNum, e2, ds, e1, zetar, om0, L, Gamma, Ap, CV, dt, MOm2, O3, t0Half, CV2, mTH, tc0);
% build in boundary condition
LeftEnd = [1 GridNum+1 2.*GridNum+1];

MatE1(LeftEnd,:) = 0;
MatE1 = MatE1 + sparse(LeftEnd,LeftEnd,1,3.*GridNum,3.*GridNum);

E1Ang = E1Ang + om0.*dt;

SourceE1(1) = 0;
SourceE1(GridNum+1) = cos(E1Ang);
SourceE1(2.*GridNum+1) = sin(E1Ang);

Ans = MatE1\SourceE1;

e1m(:,1) = Ans(1:GridNum);
e1m(:,2) = Ans(GridNum+1:2.*GridNum);
e1m(:,3) = Ans(2.*GridNum+1:3.*GridNum);
TFnan = isnan(e1m);
if sum(sum(TFnan)) >0
    bailout = 1;
end
% shift MoM
MOmMid = 0.5.*(MOm(1:GridNum-1) + MOm(2:GridNum));
    
%% Self  repulsion
fRep = zeros(GridNum,3);


rnow = [x(:,j) y(:,j) z(:,j)];
    for ii = 5:GridNum % starts at 5 so that the whirling condition doesn't instantly stop the code
    % new shit
    rMe = [x(ii,j) y(ii,j) z(ii,j)];
    delR = (rnow(:,:)-rMe);
%     % delRe1 = delR.*e1;
%     % delRe2 = delR.*e1;
%     % delRMage1 = sqrt(delRe1(:,1).^2 + delRe1(:,2).^2 + delRe1(:,3).^2);
    delRMag = sqrt(delR(:,1).^2 + delR(:,2).^2 + delR(:,3).^2);
%     % delRMage1(delRMage1>2.1*a,:)= 1E2;
    delRMag(delRMag>1.5*a,:) = 1E2;
    backNode = max(ii-closenodes,1);
    frontNode = min(ii+closenodes,GridNum);
    delRMag(backNode:frontNode,:) = 1E2;
    
    if min(delRMag) < 1.5*a
        magU = delRMag(find(delRMag==min(delRMag),1));
        unR =  abs((delR(find(delRMag==min(delRMag),1),:)-delR(ii,:))./magU);
        unRdote3 = abs(unR(1,1).*e3(ii,1) + unR(1,2).*e3(ii,2) + unR(1,3).*e3(ii,3));
        delRMagperp = delRMag./max(1-unRdote3,0.001);
        
    end
    
    if plect == 0 && min(delRMagperp)<1.5*a
        plect = 1;
        repPoint(:,:) = find(delRMag==min(delRMag));

        plectLoc = [ii,repPoint];
        TotalTime = Time(j); % modifies totaltime variable to be new shorter time!!!
        % Bail subroutine
        bailout = 1;
        
    end

    end
%     % delRMage1(backNode:frontNode,:) = 1E2;
% 
% 
%     if min(delRMag)<10 % checks to see if any nodes outside of exclusion zone are close to node ii
%             if plect == 0 
%                 plect = 1; % sets the boolean for data extraction
%             end
%             
%             repWhere = find(delRMag<2.2*a); % find close nodes
%     
%             for kk = 1:length(repWhere) % cycle through close locations and find the actual distance
%                 magU = delRMag(repWhere(kk));
%                 unR =  delR(repWhere(kk),:)./magU;
%                 
%                 fRep(ii,:) = fRep(ii,:)- unR.*(cRep/((delR(repWhere(kk))).^2)); % determines force of reuplsion
%                 
%     % fRep(ii+repWhere(kk),:) = fRep(ii+repWhere(kk),:) + unR.*(cRep/((delR(repWhere(kk))).^2)); % equal and opposite
%             end
%             
%     end
% 
% 
% 
%     end
 

   
   %% solve for tangential force
              
        DelS(1:GridNum-1) = sqrt(   ( x(2:GridNum,j) - x(1:GridNum-1,j) ).^2 ... 
                                  + ( y(2:GridNum,j) - y(1:GridNum-1,j) ).^2 ...
                                  + ( z(2:GridNum,j) - z(1:GridNum-1,j) ).^2 );
                              
%        q1 = sig0.*(1 - ds./DelS');
%        q = sig0.*(DelS'./ds - 1);
       q1 = sig0.*(DelS'./ds - 1);
       q3 = 5E5.*(ds).^2*sig0.*(DelS'./ds - 1).^3;
       
%        q(q>0) = 1*q(q>0);
%        if max(abs(DelS./ds-1))>1.6
%            keyboard
%        end
% e3 at the half grid
                          
    e3m(:,1) = ( x(2:GridNum,j) - x(1:GridNum-1,j) )./DelS';
    e3m(:,2) = ( y(2:GridNum,j) - y(1:GridNum-1,j) )./DelS';
    e3m(:,3) = ( z(2:GridNum,j) - z(1:GridNum-1,j) )./DelS';
         
%% integrate position full time step

Mxx = zetaperp.*(CV*sparse((1:GridNum),(1:GridNum),1-dragrat.*e3(:,1).^2,GridNum,GridNum));
Mxy = zetaperp.*(CV*sparse((1:GridNum),(1:GridNum),-dragrat.*e3(:,1).*e3(:,2),GridNum,GridNum)); 
Mxz = zetaperp.*(CV*sparse((1:GridNum),(1:GridNum),-dragrat.*e3(:,1).*e3(:,3),GridNum,GridNum));

Myy = zetaperp.*(CV*sparse((1:GridNum),(1:GridNum),1-dragrat.*e3(:,2).^2,GridNum,GridNum));
Myz = zetaperp.*(CV*sparse((1:GridNum),(1:GridNum),-dragrat.*e3(:,2).*e3(:,3),GridNum,GridNum));

Mzz = zetaperp.*(CV*sparse((1:GridNum),(1:GridNum),1-dragrat.*e3(:,3).^2,GridNum,GridNum));
% 
% Bxx = sparse((1:GridNum),(1:GridNum),Ap+Am.*(e1(:,1).^2 - e2(:,1).^2),GridNum,GridNum);
% Bxy = sparse((1:GridNum),(1:GridNum),Am.*(e1(:,1).*e1(:,2) - e2(:,1).*e2(:,2)),GridNum,GridNum); 
% Bxz = sparse((1:GridNum),(1:GridNum),Am.*(e1(:,1).*e1(:,3) - e2(:,1).*e2(:,3)),GridNum,GridNum);
% 
% Byy = sparse((1:GridNum),(1:GridNum),Ap+Am.*(e1(:,2).^2 - e2(:,2).^2),GridNum,GridNum);
% Byz = sparse((1:GridNum),(1:GridNum),Am.*(e1(:,2).*e1(:,3) - e2(:,2).*e2(:,3)),GridNum,GridNum);
% 
% Bzz = sparse((1:GridNum),(1:GridNum),Ap+Am.*(e1(:,3).^2 - e2(:,3).^2),GridNum,GridNum);

% TanForce = (   sparse((1:GridNum-1),(1:GridNum-1),(q-MOmMid)./ds,GridNum,GridNum) ...
%              - sparse((1:GridNum-1),(2:GridNum),(q-MOmMid)./ds,GridNum,GridNum) ...
%              - sparse((2:GridNum),(1:GridNum-1),(q-MOmMid)./ds,GridNum,GridNum) ...
%              + sparse((2:GridNum),(2:GridNum),(q-MOmMid)./ds,GridNum,GridNum) );

TanForce = (   sparse((1:GridNum-1),(1:GridNum-1),(q1+q3-MOmMid)./ds,GridNum,GridNum) ...
             - sparse((1:GridNum-1),(2:GridNum),(q1+q3-MOmMid)./ds,GridNum,GridNum) ...
             - sparse((2:GridNum),(1:GridNum-1),(q1+q3-MOmMid)./ds,GridNum,GridNum) ...
             + sparse((2:GridNum),(2:GridNum),(q1+q3-MOmMid)./ds,GridNum,GridNum) );
             
% BigMat = [ Mxx+(dt).*(O1*Bxx*L2)    Mxy+(dt).*(O1*Bxy*L2)     Mxz+(dt).*(O1*Bxz*L2);
%            Mxy+(dt).*(O1*Bxy*L2)    Myy+(dt).*(O1*Byy*L2)     Myz+(dt).*(O1*Byz*L2);
%            Mxz+(dt).*(O1*Bxz*L2)    Myz+(dt).*(O1*Byz*L2)     Mzz+(dt).*(O1*Bzz*L2) ];
%% For 2 moduli - PR
if A2 >= A1
    BigMat = [Mxx+(A2.*dt).*(O1*L2)  Mxy                      Mxz;
              Mxy                    Myy+(A2.*dt).*(O1*L2)    Myz;
              Mxz                    Myz                      Mzz+(A2.*dt).*(O1*L2) ];

    hrssX = (O1*(Om1.*e2(:,1))); hrssY = (O1*(Om1.*e2(:,2))); hrssZ = (O1*(Om1.*e2(:,3)));
else
    BigMat = [Mxx+(A1.*dt).*(O1*L2)   Mxy                      Mxz;
              Mxy                     Myy+(A1.*dt).*(O1*L2)    Myz;
              Mxz                     Myz                      Mzz+(A1.*dt).*(O1*L2) ];
    hrssX = (O1*(Om2.*e1(:,1))); hrssY = (O1*(Om2.*e1(:,2))); hrssZ = (O1*(Om2.*e1(:,3)));
end
%%
% define Source vectors

hx = Gamma.*( Om3-tc0 ).*(ys.*z2s - zs.*y2s); 
hy = Gamma.*( Om3-tc0 ).*(zs.*x2s - xs.*z2s); 
hz = Gamma.*( Om3-tc0 ).*(xs.*y2s - ys.*x2s);

Fm(:,1) = -A1.*k1.*e2(:,1) + A2.*k2.*e1(:,1);
Fm(:,2) = -A1.*k1.*e2(:,2) + A2.*k2.*e1(:,2);
Fm(:,3) = -A1.*k1.*e2(:,3) + A2.*k2.*e1(:,3);
%% CW version
%  XRHS =   Mxx*x(:,j) + Mxy*y(:,j) + Mxz*z(:,j) ...
%         + dt.*( O1*Fm(:,1) + O2*hx - TanForce*x(:,j) + Ftx);
%     
%  YRHS =   Mxy*x(:,j) + Myy*y(:,j) + Myz*z(:,j) ...
%         + dt.*( O1*Fm(:,2) + O2*hy - TanForce*y(:,j) + Fty);
%     
%  ZRHS =   Mxz*x(:,j) + Myz*y(:,j) + Mzz*z(:,j) ...
%         + dt.*( O1*Fm(:,3) + O2*hz - TanForce*z(:,j) + Ftz);
%% PR version
XRHS =   Mxx*x(:,j) + Mxy*y(:,j) + Mxz*z(:,j)  ...
    + dt.*( -O1*(Fm(:,1)) + O2*hx + Ftx - TanForce*(x(:,j))) ...
    + dt.*(A1 - A2).*hrssX ...
    + dt.*fRep(:,1);

YRHS =   Mxy*x(:,j) + Myy*y(:,j) + Myz*z(:,j) ...
    + dt.*( -O1*(Fm(:,2)) + O2*hy + Fty  - TanForce*(y(:,j))) ...
    + dt.*(A1 - A2).*hrssY ...
    + dt.*fRep(:,2);

ZRHS =   Mxz*x(:,j) + Myz*y(:,j) + Mzz*z(:,j)  ...
    + dt.*( -O1*(Fm(:,3)) + O2*hz + Ftz - TanForce*(z(:,j))) ...
    + dt.*(A1 - A2).*hrssZ ...
    + dt.*fRep(:,3); 
        
 Source = [ XRHS;
            YRHS;
            ZRHS ];
        
% build in boundary conditions

BigMat([LeftEnd, LeftEnd+1],:) = 0;
BigMat =  BigMat + sparse(LeftEnd,LeftEnd,1,3.*GridNum,3.*GridNum);
BigMat(2,1:GridNum) = L1(1,:);
BigMat(GridNum+2,GridNum+1:2.*GridNum) = L1(1,:);
BigMat(2.*GridNum+2,2.*GridNum+1:3.*GridNum) = L1(1,:);

Source(LeftEnd) = 0;
Source(2) = 1;
Source(GridNum+2) = 0;
Source(2.*GridNum+2) = 0;

Answer = BigMat\Source;
if bailout == 0
    x(:,j) = Answer(1:GridNum);
    y(:,j) = Answer(GridNum+1:2*GridNum);
    z(:,j) = Answer(2*GridNum+1:3*GridNum);
end
%% Update e1 based on a temporal rotation matrix

        xs2 = L1*x(:,j);
        ys2 = L1*y(:,j);
        zs2 = L1*z(:,j);
        
        mag = sqrt( xs2.^2 + ys2.^2 + zs2.^2 );

    e32(:,1) = xs2./mag;
    e32(:,2) = ys2./mag;
    e32(:,3) = zs2./mag;
    
    % calculating e1 dotted into the new e3

    e1dote3 = e1m(:,1).*e32(:,1) + e1m(:,2).*e32(:,2) + e1m(:,3).*e32(:,3);
    % e2dote3 = e2(:,1).*e32(:,1) + e2(:,2).*e32(:,2) + e2(:,3).*e32(:,3);
    
    % costheta is e3 dotted into new e3
    costheta = e3(:,1).*e32(:,1) + e3(:,2).*e32(:,2) + e3(:,3).*e32(:,3);
     
     e1b(:,1) =  e1m(:,1) - e1dote3.*(e3(:,1)+e32(:,1))./(1+costheta);
            
     e1b(:,2) =  e1m(:,2) - e1dote3.*(e3(:,2)+e32(:,2))./(1+costheta);
            
     e1b(:,3) =  e1m(:,3) - e1dote3.*(e3(:,3)+e32(:,3))./(1+costheta);
            
        mag = sqrt(e1b(:,1).^2+e1b(:,2).^2+e1b(:,3).^2);
        
        e1b(:,1) = e1b(:,1)./mag;
        e1b(:,2) = e1b(:,2)./mag;
        e1b(:,3) = e1b(:,3)./mag;
        
        e1 = e1b;
     % end of k loop           
    end
end
%%
% collecting data for output     
    E1dat(:,j,1) = e1(:,1); E1dat(:,j,2) = e1(:,2); E1dat(:,j,3) = e1(:,3);
    E2dat(:,j,1) = e2(:,1); E2dat(:,j,2) = e2(:,2); E2dat(:,j,3) = e2(:,3);
    

    %% Plotting
    
% quiver3 plot
% figure(2); quiver3(x(:,j),y(:,j),z(:,j),XRHS,YRHS,ZRHS);hold on; plot3(x(:,j),y(:,j),z(:,j)); hold off
% figure(2); 
% quiver3(x(:,j),y(:,j),z(:,j),e2(:,1),e2(:,2),e2(:,3));hold on; plot3(x(:,j),y(:,j),z(:,j)); hold off

if j == 2 || j == 1
Hfig = figure(1); % nothing on first iter... did pop up on 2nd
clf(1);
end

if j>=3
clf(1);
% Hfig.WindowState = 'minimized'; % %minimizes plot window so it doesn't interupt typing 
end

r(:,1) = x(:,j);
r(:,2) = y(:,j);
r(:,3) = z(:,j);
[X,Y,Z] = sphere(60);

% [X,Y,Z] = sphere(60)/;
if A1 == A2
    tubeplot2(r,Om3,rad,20,[0.3 0.7 0.4],[0.3 0.7 0.4],1);axis equal;axis([-0.4.*L L+0.1 -0.8.*L 0.8.*L -0.8.*L 0.8.*L]);view([2 15])
    hold on
    h1 = patch(surf2patch(rad(GridNum).*X+r(GridNum,1),rad(GridNum).*Y+r(GridNum,2),rad(GridNum).*Z+r(GridNum,3),Z)); % patch last point with sphere
    set(h1,'FaceColor',[0.3 0.7 0.4],'EdgeColor','none'); 
    h2 = patch(surf2patch(rad(1).*X+r(1,1),rad(1).*Y+r(1,2),rad(1).*Z+r(1,3),Z)); % patch 1st point with sphere
    set(h2,'FaceColor',[0.8 0.3 0.3],'EdgeColor','none');box off;axis off; 

else
    ribbonplot(r,e1,e2,1.0.*a,0.2.*a,[0.8 0.3 0.4],1);axis equal;axis([-0.4.*L L+0.1 -0.8.*L 0.8.*L -0.8.*L 0.8.*L]);view([2 15])
end
hold on

box off;axis off;
lightangle(40,40);
set(gcf,'Color',[1 1 1]);
% quiver3(e2(:,1),e2(:,2),e2(:,3),x(:,j),y())
if plotbox == 1    
    view([45 15])
    ySc = 0.70;
    zSc = 0.70; 
    xLHlimit = -0.4;
    xRHlimit = 1.1;
    box1p1 = [xLHlimit -L.*ySc -L.*zSc]; % -,-,-
    box1p2 = [xLHlimit +L.*ySc -L.*zSc]; % -,+,-
    box1p3 = [xLHlimit +L.*ySc +L.*zSc]; % -,+,+
    box1p4 = [xLHlimit -L.*ySc +L.*zSc]; % -,-,+
    box1plot = [box1p1;box1p2;box1p3;box1p4];
    box2p1 = [xRHlimit -L.*ySc -L.*zSc];
    box2p2 = [xRHlimit +L.*ySc -L.*zSc];
    box2p3 = [xRHlimit +L.*ySc +L.*zSc];
    box2p4 = [xRHlimit -L.*ySc +L.*zSc];
    box2plot = [box2p1;box2p2;box2p3;box2p4];
    
    for oo = 1:4
        plist1 = mod(oo-1,4)+1;
        plist2 = mod(oo,4)+1;
        linenow = [box1plot(plist1,:);box1plot(plist2,:)];
        plot3(linenow(:,1), linenow(:,2), linenow(:,3), 'LineWidth',2.0,'color',[0,0,0]+0.4)
        linenow = [box2plot(plist1,:);box2plot(plist2,:)];
        plot3(linenow(:,1), linenow(:,2), linenow(:,3), 'LineWidth',2.0,'color',[0,0,0]+0.4)
        linenow = [box1plot(oo,:);box2plot(oo,:)];
        plot3(linenow(:,1), linenow(:,2), linenow(:,3), 'LineWidth',2.0,'color',[0,0,0]+0.4)
    
    end
    box1p1 = [xLHlimit -L.*ySc -L.*zSc]; % -,-,-
    box1p2 = [xLHlimit +L.*ySc -L.*zSc]; % -,+,-
    box1p3 = [xLHlimit +L.*ySc +L.*zSc]; % -,+,+
    box1p4 = [xLHlimit -L.*ySc +L.*zSc]; % -,-,+
    box1plot = [box1p1;box1p2;box1p3;box1p4];
    box2p1 = [xRHlimit -L.*ySc -L.*zSc];
    box2p2 = [xRHlimit +L.*ySc -L.*zSc];
    box2p3 = [xRHlimit +L.*ySc +L.*zSc];
    box2p4 = [xRHlimit -L.*ySc +L.*zSc];
    box2plot = [box2p1;box2p2;box2p3;box2p4];

    patch([box1p1(1) box1p2(1) box1p3(1) box1p4(1)],[box1p1(2) box1p2(2) box1p3(2) box1p4(2)],[box1p1(3) box1p2(3) box1p3(3) box1p4(3)],[0.99,0.99,0.99])
    patch([box1p1(1) box1p2(1) box2p2(1) box2p1(1)],[box1p1(2) box1p2(2) box2p2(2) box2p1(2)],[box1p1(3) box1p2(3) box2p2(3) box2p1(3)],[0.99,0.99,0.99])
    patch([box1p2(1) box1p3(1) box2p3(1) box2p2(1)],[box1p2(2) box1p3(2) box2p3(2) box2p2(2)],[box1p2(3) box1p3(3) box2p3(3) box2p2(3)],[0.99,0.99,0.99])
    
    plot3(x(:,j),y(:,j),zeros(GridNum,1)-L.*zSc, 'LineWidth',3.0,'color',[0,0,0.5]+0.2);
    plot3(zeros(GridNum,1)+xLHlimit,y(:,j),z(:,j), 'LineWidth',3.0,'color',[0,0,0.5]+0.2);
    plot3(x(:,j),zeros(GridNum,1)+L.*ySc,z(:,j), 'LineWidth',3.0,'color',[0,0,0.5]+0.2);
    
end
% if j == NSteps
% 
%     %% save last plots
% figure(89);
% for pp = 1:10
% plot3(x(:,j-pp),y(:,j-pp),z(:,j-pp));hold on;
% end
% hold off
vidtitle = strcat('TwirlWhirl-Ribbon_',num2str(omRat, '%6.3g'),'omRat_','Arat_',num2str(Arat, '%6.3g'),'_',num2str(dt, '%6.3g'),'dt_',dtstr,'THF_',num2str(THF, '%6.3g'),'_eta', num2str(eta*1000, '%6.3g'));
% figtitle = strcat(vidtitle,'.fig');
% savefig(figtitle)
% 
% 
% end

%%
ElEner1(j-1,1) = sum(A1.*ds./2.*(Om1(2:GridNum)-k1).^2);
ElEner2(j-1,1) = sum(A2.*ds./2.*(Om2(2:GridNum)-k2).^2);
ElEner3(j-1,1) = sum(Ap.*ds.*Gamma./2.*(Om3(2:GridNum)-tc0).^2);
ElEnerStr(j-1,1) = sum((sig0./2./ds).*(DelS(2:end)-ds).^2);
Etot(j-1,1) = ElEner1(j-1,1)+ElEner2(j-1,1)+ElEner3(j-1,1)+ElEnerStr(j-1,1);
title(strcat('Elastic Energy = ',num2str(round(Etot(j-1,1),4, 'significant'),'%3.2f'),...
    ', Time:',num2str(Time(j),'%2.4f'),...
    ', \omega_{o}:',num2str(om0,'%3.2f'),...
    ', \omega_{c}:', num2str(8.9*omFactor,'%3.2f')),...
    'FontSize',12);

% making plot large for graphics to be clearer
% Hfig = figure(1);
hold off
set(Hfig, 'Position',  [100, 100, 600, 400]);

mov(:,j-1) = getframe(Hfig);
%% Number of Times Crossed Zero
yShift = y(3:GridNum).*y(2:GridNum-1);
zShift = z(3:GridNum).*z(2:GridNum-1);
yCros = find(yShift<0);
zCros = find(zShift<0);
if isempty(yCros)
    yZeros = 0;
else
    yZeros = length(yCros);
end
if isempty(zCros)
    zZeros = 0;
else
    zZeros = length(zCros);
end

maxDelS(j-1,1) = max(abs(DelS'./ds-1));
%%
display = ["max 1-Dels/ds", "step", "om0","Length","Last X", "Time", "Energy", "e3x";...
num2str(max(abs(1-DelS./ds))), strcat(string(Iter),'/',string(NSteps)),om0,sum(DelS) ,x(end,j), Time(j), Etot(j-1), e3(end,1)]  

if te3flip == 9999 && e3(end,1)<0
    te3flip = Time(j);
end


end %end of j loop
%% Make a video 

writerObj = VideoWriter(strcat(vidtitle));
writerObj.FrameRate = 20; %FPS
%Open Video Writer
open(writerObj);
for i=1:length(mov)
    %convert image to frame of movie
    frame = mov(i);
    writeVideo(writerObj,frame);
end
%close the video writer
close(writerObj);


%% Create tables of data
eBend = ElEner1 + ElEner2;

Tx = table(x);
Ty = table(y);
Tz = table(z);
TOm = table(Om1,Om2,Om3);

% TabTFtx = table(TFTx);
E1xTab = table(E1dat(:,:,1));
E1yTab = table(E1dat(:,:,2));
E1zTab = table(E1dat(:,:,3));
E2xTab = table(E2dat(:,:,1));
E2yTab = table(E2dat(:,:,2));
E2zTab = table(E2dat(:,:,3));

omC = 8.9.*(Ap)./((L.^2).*zetar);

xLast = L-x(end,2:NSteps);
xLShift = x(end,2:end-1).*x(end,3:end);

flipLoc = find(xLShift<0);

if isempty(flipLoc)
    tW1 = 999999;
    twirlWhen = 0;      
else
    twirlWhen = Time(flipLoc);
    tW1 = twirlWhen(1);
end
tTW = table(twirlWhen);

xEndPt = x(end,end);
timesteps = 51;
    timeRun = Time(min(NSteps,timesteps)) - Time(1);

lastAvgEner(1,1) = mean(Etot(max(1,NSteps-timesteps):NSteps-1));
lastEner(1,1) = Etot(end,1);
enerTrend = (lastEner - lastAvgEner)./ timeRun;
maxEner = max(max(Etot));
maxmaxDelS = max(maxDelS);
avgmaxDelS = mean(maxDelS);
maxTwEner = max(ElEner3);
maxBendEner = max(eBend);  
endTwEner = ElEner3(end);
endBendEner = eBend(end); 
AmOvAp = Am./Ap;
farRight = find(x(:,j) == max(x(:,j)));

infodata = [omRat, AmOvAp, xEndPt, tW1, te3flip,...
            yZeros, zZeros,maxmaxDelS, avgmaxDelS, maxEner,...
            lastAvgEner, lastEner,enerTrend, maxTwEner, maxBendEner,...
            endTwEner, endBendEner,omC,om01,time1,...
            om02,time2,L, a, eta,...
            zetar, Gamma,A1, A2, sig0,...
            GridNum,dt,THF, TotalTime,bailout,...
            plect,plectLoc(1),plectLoc(2),farRight];

 EnergyTable = table(tlist,ElEner1,ElEner2,eBend,ElEner3,ElEnerStr,Etot,omList,xLast',...
'VariableNames', {'Time','E1','E2','E1+E2','E3','ETan','ETotal','omList','xDelta'});
     

 InfoTable = table(infodata(:,1),infodata(:,2),infodata(:,3),infodata(:,4),infodata(:,5),...
                  infodata(:,6),infodata(:,7),infodata(:,8),infodata(:,9),infodata(:,10),...
                  infodata(:,11),infodata(:,12),infodata(:,13),infodata(:,14),infodata(:,15),...
                  infodata(:,16),infodata(:,17),infodata(:,18),infodata(:,19),infodata(:,20),...
                  infodata(:,21),infodata(:,22),infodata(:,23),infodata(:,24),infodata(:,25),...
                  infodata(:,26),infodata(:,27),infodata(:,28),infodata(:,29),infodata(:,30),...
                  infodata(:,31),infodata(:,32),infodata(:,33),infodata(:,34),infodata(:,35),...
                  infodata(:,36),infodata(:,37),infodata(:,38),infodata(:,39),...
                  'VariableNames',{'omRatio','Am/Ap','xEndPt','XEndTwirlTime','TimeE3Flip',...
                                   'yZeros', 'zZeros','MaxMaxStretch','AvgMaxStretch','MaxEner',...
                                   'FinalAvgEner','EndEnergy','ETrend','maxTwistEner','maxBendEner',...
                                   'endTwEner','endBendEner','omC','om01','time1',...
                                   'om02','time2','Length','radius','eta',...
                                   'zetar','Gamma','A1','A2','sig0',...
                                   'GridNum','dt','THF','TotalTime','bailout',...
                                   'Plect', 'PlectLoc1', 'PlectLoc2','MostRightNode'});


%% Write tables to excel file
writetable(EnergyTable,strcat(vidtitle,'_data.xlsx'),'Sheet','AllEnergy')
writetable(InfoTable,strcat(vidtitle,'_data.xlsx'),'Sheet','Info');
writetable(TOm,strcat(vidtitle,'_data.xlsx'),'Sheet','Omegas');
% writetable(Tx,strcat(vidtitle,'_data.xlsx'),'Sheet','x');
% writetable(Ty,strcat(vidtitle,'_data.xlsx'),'Sheet','y');
% writetable(Tz,strcat(vidtitle,'_data.xlsx'),'Sheet','z');
% writetable(E1xTab,strcat(vidtitle,'_data.xlsx'),'Sheet','E1x');
% writetable(E1yTab,strcat(vidtitle,'_data.xlsx'),'Sheet','E1y');
% writetable(E1zTab,strcat(vidtitle,'_data.xlsx'),'Sheet','E1z');
% % writetable(E2xTab,strcat(vidtitle,'_data.xlsx'),'Sheet','E2x');
% % writetable(E2yTab,strcat(vidtitle,'_data.xlsx'),'Sheet','E2y');
% % writetable(E2zTab,strcat(vidtitle,'_data.xlsx'),'Sheet','E2z');



end

