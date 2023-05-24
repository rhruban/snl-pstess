% make PST play-in data from PST simulation data
% Ron Hruban 2023

%% Housekeeping
close all;
clear all;
clc;

%% setup
load('PSTdata_d2asbeghp.mat');

t = g.t; % time vector
mac_num = 1; % machine to be validated number
genb = g.bus.bus_int(mac_num); % Bus Number of machine to be validated
pmub = g.bus.bus_int(10); % "PMU location bus"

% !!! check line reactance in the next line !!!( it is not saved in g)
disp('Did you verify line impedance?');
line_imp=1i*0.0167; %check line impedance
I = (g.bus.bus_v(genb,:) - g.bus.bus_v(pmub,:))./(line_imp); % from gen to pmu bus
S = g.bus.bus_v(pmub,:).*conj(I); % into PMU bus
P = real(S);
Q = imag(S);

%% Plot overall PMU data

figure

subplot(511)
plot(t(2:end)./600,1./diff(t),'k','linewidth',2)
ylabel('f_s (sps)')

subplot(512)
plot(t./600,g.freq.bus_freq(genb,:),'k','linewidth',2);
ylabel('f_d (Hz)')

subplot(513)
plot(t./600,abs(g.bus.bus_v(pmub,:)),'k','linewidth',2);
hold on;
plot(t./600,abs(g.bus.bus_v(genb,:)),'r--','linewidth',2);
hold off;
ylabel('V (kV)')

subplot(514)
plot(t./600,P,'k','linewidth',2)
hold on;
plot(t./600,g.mac.pelect(1,:),'r--','LineWidth',2)
hold off;
ylabel('P (MW)')

subplot(515)
plot(t./600,Q,'k','linewidth',2)
hold on;
plot(t./600,g.mac.pelect(mac_num,:),'r--','LineWidth',2)
hold off;
legend('PMU','Gen','Location','best')
ylabel('Q (MVAR)')
xlabel('Time (min.)')

x = get(gcf,'Position');
set(gcf,'Position',[x(1) 0.5*x(2) x(3) 1.5*x(4)]);

%% Build V
% Initial conditions at PMU
disp('PMU initial conditions')
P0 = P(1)
Q0 = Q(1)
V0 = abs(g.bus.bus_v(pmub,1))
ang0 = angle(g.bus.bus_v(pmub,1))/pi*180

%% Initial bus matrix
y=[-P0;-Q0];
ymhos=-1/(line_imp);
ymag=abs(ymhos);
yang=angle(ymhos)/pi*180;
f=@(x) [V0*x(1)*ymag*cosd(ang0-x(2)-yang)+V0*V0*ymag*cosd(ang0-ang0-yang+180);...
        V0*x(1)*ymag*sind(ang0-x(2)-yang)+V0*V0*ymag*sind(ang0-ang0-yang+180)];

x0=[1;0];
J=@(x) [V0*ymag*cosd(ang0-x(2)-yang),V0*x(1)*ymag*sind(ang0-x(2)-yang);...
        V0*ymag*sind(ang0-x(2)-yang),-1*V0*x(1)*ymag*cosd(ang0-x(2)-yang)];

x=x0;
tol=1e-8;
delx=1;
count = 0;
myiter = 10000;
while count<=myiter 
    JJ=J(x);
    yhat=f(x);
    delx=JJ\(y-yhat);
    x=x+delx;
    count=count+1;
    if norm(delx)<tol
        break
    end
end
count
if count >= myiter
    error('Newton Raphson Did Not Converge')
end

V20 = x(1)
ang20 = wrapTo180(x(2))
if V20 < 0
    V20 = -V20
    ang20 = wrapTo180(ang20+180)
end
P20 = V0*V20*ymag*cosd(ang20-ang0-yang)+V20*V20*ymag*cosd(ang20-ang20-yang+180)
Q20 = V0*V20*ymag*sind(ang20-ang0-yang)+V20*V20*ymag*sind(ang20-ang20-yang+180)


vtemp = g.bus.bus_v(pmub,:);
ttemp = t;

if 0>(sum(isnan(vtemp))+sum(isinf(vtemp))+sum(isnan(ttemp))+sum(isinf(ttemp)))
    error('NAN or INF in v or t')
end

t = t(1):1/600:t(end);
v = interp1(ttemp,vtemp,t,'spline','extrap');
p = interp1(ttemp,P,t,'spline','extrap'); % into PMU
q = interp1(ttemp,Q,t,'spline','extrap'); % into PMU
i = interp1(ttemp,I,t,'spline','extrap'); % into PMU



mytitle = ['d_PST_Val_playin.mat']

%save(mytitle,'P0',"Q0","V0","ang0","t","v","p","q","i","Q20","P20","V20","ang20");

figure()
plot(g.t,abs(g.bus.bus_v(pmub,:)));
hold on;
plot(t,abs(v),'--');
hold off;