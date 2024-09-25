function [SourceE1, MatE1, Om3]=MakeE1Mats(GridNum, e2, ds, e1, zetar, om0, L, Gamma, Ap, CV, dt, MOm2, O3, t0Half, CV2, mTH, tc0)
% inputs = GridNum, e2, ds, e1, T0R, T0L, zetar, om0, L, Gamma, Ap, CV, dt,
% MOm2, O3, t0Half, CV2
% outputs SourceE1, MatE1, Om3
% if THF == 1
%     mTH = TFRot.*randn(GridNum,1);
% end

    E2xMat = sparse((1:GridNum),(1:GridNum),e2(:,1),GridNum,GridNum);
    E2yMat = sparse((1:GridNum),(1:GridNum),e2(:,2),GridNum,GridNum);
    E2zMat = sparse((1:GridNum),(1:GridNum),e2(:,3),GridNum,GridNum);
    
    TwX = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,1)/ds,GridNum-1,GridNum) ...
          + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,1)/ds,GridNum-1,GridNum);
    TwY = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,2)/ds,GridNum-1,GridNum) ...
          + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,2)/ds,GridNum-1,GridNum);
    TwZ = - sparse((1:GridNum-1),(1:GridNum-1),0.5.*e2(2:GridNum,3)/ds,GridNum-1,GridNum) ...
          + sparse((1:GridNum-1),(2:GridNum),0.5.*e2(1:GridNum-1,3)/ds,GridNum-1,GridNum);

Tw = [ TwX, TwY, TwZ ];
e1Vec = [ e1(:,1); e1(:,2); e1(:,3) ];
      
Om3Half = Tw*e1Vec;
Om3(1) = -zetar.*om0.*L./Gamma./Ap;
Om3(2:GridNum-1) = 0.5.*(Om3Half(1:GridNum-2) + Om3Half(2:GridNum-1));
Om3(GridNum) = tc0;

TwMat = [ E2xMat*O3;
          E2yMat*O3;
          E2zMat*O3 ];

MatE1 = CV2 - (Gamma.*dt./zetar).*(TwMat*Tw);
      
SourceE1 = CV2*e1Vec - (Ap.*L.*Gamma.*dt./zetar).*(TwMat*t0Half);

SourceE1 = SourceE1 + (dt./zetar).*[ CV*((MOm2+mTH).*e2(:,1)); (CV*(MOm2+mTH).*e2(:,2)); (CV*(MOm2+mTH).*e2(:,3)) ];
end