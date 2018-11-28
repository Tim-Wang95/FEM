close all; clear; clc;
%% Input matrix value

% constant values
E = 100; A1 = 1; A2 = 2; L = 10;
alphaT = 0.001; deltaT = 30;

% assign values
Elast = [E;E;E;E;E]; area = [A1;A1;A2;A2;A2]; leng = [(sqrt(2)*L);(sqrt(2)*L);L;L;L];
ang = [pi/4;3*pi/4;0;0;0];

% local stiffness
K = zeros(4,4,5);
for mem=1:5
   K(:,:,mem) = Elast(mem)*area(mem)/leng(mem)*localstiffness(ang(mem)); 
end

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

% compute the member forces
u = zeros(4,1,5); f = zeros(5,1); T = zeros(4,4,5);
for mem=1:5
    T(:,:,mem) = transformM(ang(mem)); u(:,:,mem) = T(:,:,mem)*U(EFT(:,mem));
    f(mem) = Elast(mem)*area(mem)/leng(mem)*(u(3,1,mem)-u(1,1,mem));
end


%% Functions and Definitions
% local matrix function
function Ke = localstiffness(theta)
    Ke(1,1:4) = [cos(theta)^2 sin(theta)*cos(theta) -cos(theta)^2 -sin(theta)*cos(theta)];
    Ke(2,1:4) = [sin(theta)*cos(theta) sin(theta)^2 -sin(theta)*cos(theta) -sin(theta)^2];
    Ke(3,1:4) = [-cos(theta)^2 -sin(theta)*cos(theta) cos(theta)^2 sin(theta)*cos(theta)];
    Ke(4,1:4) = [-sin(theta)*cos(theta) -sin(theta)^2 sin(theta)*cos(theta) sin(theta)^2];
end

% coord. transform matrix
function T = transformM(theta)
    T(1,1:4) = [cos(theta) sin(theta) 0 0];
    T(2,1:4) = [-sin(theta) cos(theta) 0 0];
    T(3,1:4) = [0 0 cos(theta) sin(theta)];
    T(4,1:4) = [0 0 -sin(theta) cos(theta)];
end
