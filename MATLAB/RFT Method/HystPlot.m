% x = zeros(20,1);
% y = zeros(20,1);
% z = zeros(20,1);
AMS = [0,.2,.4,.6,.8,.9];
for ii = 1:length(AMS)
plot3(x(y==AMS(ii)),y(y==AMS(ii)),z(y==AMS(ii)),'.')
hold on
grid on
end

% mesh(x,y,z)
% or surf(x,y,z,'FaceAlpha',0.5)