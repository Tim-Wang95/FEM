close all; clear; clc;
%% Input matrix value

% constant values
E = 100; A1 = 1; A2 = 2; L = 10;
alphaT = 0.001; deltaT = 30;

% local stiffness
K(:,:,1) = E*A1/(sqrt(2)*L)*localstiffness(pi/4);
K(:,:,2) = E*A1/(sqrt(2)*L)*localstiffness(3*pi/4);
K(:,:,3) = E*A2/L*localstiffness(0);
K(:,:,4) = E*A2/L*localstiffness(0); 
K(:,:,5) = E*A2/L*localstiffness(0);

% element freedom table
EFT(:,1) = [1 2 7 8]; EFT(:,2) = [3 4 9 10]; EFT(:,3) = [5 6 7 8];
EFT(:,4) = [7 8 9 10]; EFT(:,5) = [9 10 11 12];

% global stiffness
dof = 12; K_ = zeros(dof);
for m=1:length(K(1,1,:))
   for i=1:4
       for j=1:4
           K_(EFT(i,m),EFT(j,m)) = K_(EFT(i,m),EFT(j,m))+K(i,j,m); 
       end
   end
end

%% Calculation

% for node 7-10
f = [0;-5;5*sqrt(2);-5-5*sqrt(2)]; k = K_(7:10,7:10);
u = inv(k)*f;
U = [[0;0;0;-2;0;0];u;0;0];

% thermal forces
F_T = zeros(12,1); m_T = [-1;0;1;0];
for t=3:5
    F_T(EFT(:,t)) = F_T(EFT(:,t))+E*A2*alphaT*deltaT*m_T;
end
    
% total reaction forces
K_*U-F_T

%% Local matrix function
function Ke = localstiffness(theta)
    Ke(1,1:4) = [cos(theta)^2 sin(theta)*cos(theta) -cos(theta)^2 -sin(theta)*cos(theta)];
    Ke(2,1:4) = [sin(theta)*cos(theta) sin(theta)^2 -sin(theta)*cos(theta) -sin(theta)^2];
    Ke(3,1:4) = [-cos(theta)^2 -sin(theta)*cos(theta) cos(theta)^2 sin(theta)*cos(theta)];
    Ke(4,1:4) = [-sin(theta)*cos(theta) -sin(theta)^2 sin(theta)*cos(theta) sin(theta)^2];
end


