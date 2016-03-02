% Model of glacier development and stabilization. Ryan Stoner Feb 29, 2016
% For Geology modeling seminar
clear
%% Initialize

figure(1)
clf



zmax = 2550;            % m, maximum altitude
ELAave = 2200;          % m
%ELA = ELAave;           % m
s = 0.05;               % slope
dx = 100;               % m
xmax = 20000;           % m
x = 0:dx:xmax;          % m
zbas = zmax* 5.^(-x/50000);          % m
z = zbas;         % m
AMP = 60;               % m, amplitude of change

gamma = 0.01;           % m/yr

dt = 0.002;                 % yr
tmax = 150;              % yr
t = 0:dt:tmax;            % yr
nplots = 20;
tplot = tmax/nplots;

N = 3;

icedens = 917;          % kg/m^3
g = 9.8;                % m/s^2
W = 100;
h = zeros(1,length(x));   % m, initial thickness of ice
usl = 0.01;              % m/yr, sliding velocity of ice
A = 2.1*10^-16;         % yr^-1,Pa^-3

nframe = 0;
%% Loop
% First calculate the ELA adjusted to sinusoidal variations
% Calculate the mass balance
% Calculate slope of ice
% Interpolate hedge values so that matrix dimensions will match
% Calculate flux and pad edges
% Calculate change in thickness of ice
% Recalculate thickness and elevation of ice
% Plot ice, basement, add labels including one for time

imax = length(t);
for i=1:imax
ELA = ELAave+60*sin(t(i));
b = gamma*(z-ELA);
dzdx = abs(diff(z)/dx);

dHdx = diff(h)/dx;
hedge = h(1:length(x)-1)+0.5*dHdx;

Q = (usl*hedge)+ A*(icedens*g*abs(dzdx)).^3.*(hedge.^5)/5;

Q = [0 Q 0];
dHdt = b - 1/W*(diff(Q)/dt);



% for i =1:imax

% B = f(z) % not a function of bedrock

% H
% z = 

h = h+dHdt*dt;
h = max(h,0);
z = zbas + h;

if(mod(t(i),tplot)==0)
nframe = nframe+1;
    figure(1)
plot(x/1000,zbas,'k')
hold on
plot(x/1000,z,'r')

%plot(x/1000,ELA*ones(1,length(z)),'b')
trounded = round(t(i),0);
strtime = num2str(trounded);

txt = uicontrol('Style','text',...
      'Position',[430 0 80 20],...
      'Fontsize',12,...
      'FontWeight','Bold',... 
      'HorizontalAlignment','left',...
      'String',['Time (yr): ' strtime]);

xlabel('distance (km)')
ylabel('elevation (m)')
title('Inception and Growth of Glacier over Time')

pause(0.1)

end

% end


end
hold off