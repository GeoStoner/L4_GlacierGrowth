% Model of glacier development and stabilization. Ryan Stoner Feb 29, 2016
% For Geology modeling seminar
clear
%% Initialize


% Load .txt file of pleistocene data. 
load pleist_del18O.txt;

tmax = 50;                                  % kyr
tmaxind = find(pleist_del18O(:,1)<=tmax);   % find times less than tmax
pltime = flipud(pleist_del18O(tmaxind,1));  % kyr, flipped past to present
pldel18 = (pleist_del18O(tmaxind,2));       % del180 values,Past to present

dt = 0.002;                                 % yr
t = (tmax*1000:-dt:0)';                     % yr
nplots = 1000;                              % number of plots
tplot = tmax*1000/nplots;                   % yr, plot interval

% Interpolate del18O values to smaller times and then subtract mean
pldel18sm = interp1(pltime,pldel18,t,'spline');

pldel18meaned = pldel18sm - mean(pldel18sm);
f = 2/(max(pldel18meaned)-min(pldel18meaned));

% Try to have pldel18norm values range from -1 to 1.
pldel18norm = pldel18meaned * f;

% Plot of del18O values with interpolated values
% figure(1)
% plot(pltime,pldel18,'ro',t,pldel18sm,'m:.');
% xlim([0 tmax]);

zmax = 2550;                        % m, maximum altitude
ELAave = 2200;                      % m
s = 0.05;                           % slope
dx = 100;                           % m
xmax = 20000;                       % m
x = 0:dx:xmax;                      % m
zbas = zmax* 5.^(-x/50000);         % m
z = zbas;                           % m
AMP = 60;                           % m, amplitude of change
gamma = 0.01;                       % m/yr


N = 3;                              % power of strain vs. stress for ice

icedens = 917;                      % kg/m^3
g = 9.8;                            % m/s^2
W = 100;                            % m, width of glacier
h = zeros(1,length(x));             % m, initial thickness of ice
usl = 0.01;                         % m/yr, sliding velocity of ice
A = 2.1*10^-16;                     % yr^-1,Pa^-3


nframe = 0;

%% Loop
% First calculate the ELA adjusted to del18O values
% Calculate the mass balance
% Calculate slope of ice
% Interpolate hedge values so that matrix dimensions will match
% Calculate flux and pad edges
% Calculate change in thickness of ice
% Recalculate thickness and elevation of ice
% Plot ice, basement, add labels including one for time

imax = length(t);
for i=1:imax
    
ELA = ELAave+AMP*pldel18norm(i);
b = gamma*(z-ELA);
dzdx = abs(diff(z)/dx);

dHdx = diff(h)/dx;
hedge = h(1:length(x)-1)+0.5*dHdx;

Q = (usl*hedge)+ A*(icedens*g*abs(dzdx)).^3.*(hedge.^5)/5;
Q = [0 Q 0];
dHdt = b - 1/W*(diff(Q)/dt);

h = h+dHdt*dt;
h = max(h,0);
z = zbas + h;

if(mod(t(i),tplot)<=0.000001)
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
      'String',['Time: ' strtime]);

xlabel('distance (km)')
ylabel('elevation (m)')
title('Inception and Growth of Glacier over Time')

pause(0.01)
hold off
end

end
