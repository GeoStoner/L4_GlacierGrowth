% Model of glacier development and stabilization. Ryan Stoner Feb 21, 2016
% For Geology modeling seminar
clear all
%% Initialize



zmax = 2550;            % m, maximum altitude
ELA = 2500;             % m
s = 0.2;               % slope
dx = 1;               % m
xmax = 200;           % m
x = 0:dx:xmax;          % m
zbas = zmax - s*x;          % m
z = zbas;         % m

gamma = 0.01;           % m/yr

dt = 1;                 % yr
tmax = 2000;              % yr
t = 0:dt:tmax;            % yr

N = 3;

icedens = 917;          % kg/m^3
g = 9.8;                % m/s^2
W = 100;
h = zeros(1,length(x));   % m, initial thickness of ice
usl = 0.1;              % m/yr, sliding velocity of ice
A = 2.1*10^-16;         % yr^-1,Pa^-3
%% Loop
imax = length(t);
for i=1:imax

b = gamma*((zmax-zbas)-ELA)
dzdx = abs(diff(z)/dx);

dHdx = diff(h)/dx;


hedge = h(1:length(x)-1)+0.5*dHdx;

Q = usl*hedge+ A*(icedens*g*abs(dzdx)).^3.*(hedge.^5)/5;

Q = [0 Q 0];
dHdt = b - 1/W*(diff(Q)/dt);



% for i =1:imax

% B = f(z) % not a function of bedrock

% H
% z = 

h = h+dHdt*dt;
z = zbas +h;
%hold on
%plot(x,z)
% end


end
%hold off