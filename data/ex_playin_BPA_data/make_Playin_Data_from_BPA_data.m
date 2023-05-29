%% Make playin data from bpa data
% Ron Hruban - May 2023
% Note BPA Data has weird negative sines on q and i that need to be checked
% and fixed
% This is setup for Centralia G1

%% Housekeeping
close all;
clear all;
clc;

%% setup
load('centralia_data.mat','D');

% pu settings
Vbase = 500e3/sqrt(3);
Sbase = 870;
Ibase = Sbase*1e5/3/Vbase;
fbase = 60;
kGen = 1;

% desired data range
tRange = [20387.3 20417];
% tRange = [13494.2 13519.2];

nR = find(D.t{kGen}>=tRange(1) & D.t{kGen}<=tRange(2));
k = kGen;

%% Plot overall PMU data
% D.Vt{k}=D.Vt{k}*exp(1i*pi);
% D.I{k}=D.I{k}*exp(1i*pi);

% f = angle(D.Vt{k}(2:end)./D.Vt{k}(1:end-1))./(2*pi*diff(D.t{k}));
% f = [0; f];
fi = angle(D.I{k}(2:end)./D.I{k}(1:end-1))./(2*pi*diff(D.t{k}));
fi = [0; fi];
Ip=real(D.I{k});
Iq=imag(D.I{k});

figure

subplot(511)
plot(D.t{k}(2:end)./60,1./diff(D.t{k}),'k','linewidth',2)
ylabel('f_s (sps)')

subplot(512)
plot(D.t{k}./60,D.f{k},'k',D.t{k}([nR(1) nR(end)])./60,D.f{k}([nR(1) nR(end)]),'ro','linewidth',2,'markersize',10);
ylabel('f_d (Hz)')

subplot(513)
plot(D.t{k}./60,abs(D.Vt{k})./1e3,'k',D.t{k}([nR(1) nR(end)])./60,abs(D.Vt{k}([nR(1) nR(end)]))./1e3,'ro','linewidth',2,'markersize',10)
ylabel('|V| (kV)')

subplot(514)
plot(D.t{k}./60,D.P{k},'k',D.t{k}([nR(1) nR(end)])./60,D.P{k}([nR(1) nR(end)]),'ro','linewidth',2,'markersize',10)
ylabel('P (MW)')

subplot(515)
plot(D.t{k}./60,D.Q{k},'k',D.t{k}([nR(1) nR(end)])./60,D.Q{k}([nR(1) nR(end)]),'ro','linewidth',2,'markersize',10)
ylabel('Q (MVAR)')
xlabel('Time (min.)')


x = get(gcf,'Position');
set(gcf,'Position',[x(1) 0.5*x(2) x(3) 1.5*x(4)]);

%% Plot
figure();
subplot(311);
hold on;
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),abs(D.Vt{k}(nR(1) :nR(end)))/Vbase,'k','LineWidth',2);
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),abs(D.I{k}(nR(1) :nR(end)))/Ibase,'r','LineWidth',2);
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),abs(-D.I{k}(nR(1) :nR(end)))/Ibase,'b--','LineWidth',2);
title('Use this plot to determine the direction of I (I data)');
legend('Voltage','Current','-Current','location','best');
ylabel('Magnitude (PU)')
hold off;

subplot(312);
hold on;
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),angle(D.Vt{k}(nR(1) :nR(end)))*180/pi,'k','LineWidth',2);
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),angle(D.I{k}(nR(1) :nR(end)))*180/pi,'r','LineWidth',2);
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),angle(-D.I{k}(nR(1) :nR(end)))*180/pi,'b--','LineWidth',2);
legend('Voltage','Current','-Current','location','best');

ylabel('Angle (Degrees)')
hold off;

subplot(313);
hold on;
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),wrapTo180((angle(D.Vt{k}(nR(1) :nR(end)))-angle(D.I{k}(nR(1) :nR(end))))*180/pi),'r','LineWidth',2);
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),wrapTo180((angle(D.Vt{k}(nR(1) :nR(end)))-angle(-D.I{k}(nR(1) :nR(end))))*180/pi),'b','LineWidth',2);
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),wrapTo180(angle(D.P{k}(nR(1) :nR(end))/Sbase+1i*D.Q{k}(nR(1) :nR(end))/Sbase)*180/pi),'k','LineWidth',2);
plot(D.t{k}(nR(1) :nR(end))-D.t{k}(nR(1)),wrapTo180(angle(D.P{k}(nR(1) :nR(end))/Sbase+1i*-D.Q{k}(nR(1) :nR(end))/Sbase)*180/pi),'g--','LineWidth',2);
legend('V & I data','V & -I data','P & Q data','P & -Q data')
ylabel('Power Factor Angle (Degrees)');
xlabel('Time (sec.)');


% D.Vt{k}=D.Vt{k}*exp(1i*pi);

%% Build V
P0 = D.P{k}(nR(1))/Sbase
Q0 = -D.Q{k}(nR(1))/Sbase %negate for BPA data based on plots
f0 = D.f{k}(nR(1))/fbase
V0 = abs(D.Vt{1}(nR(1))/Vbase)
ang0 = angle(D.Vt{1}(nR(1))/Vbase)/pi*180

%% Initial bus matrix
y=[-P0;-Q0];
ymhos=-1/(0.00348+ 1i*0.179742);
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
    error("Whoa buddy, that Newton Raphson didn't have enough iterations to converge");
end



if x(1) < 0
    V20 = -x(1)
    ang20 = wrapTo180(wrapTo180(x(2))+180)
else
    V20 = x(1)
    ang20 = wrapTo180(x(2))

end
P20 = V0*V20*ymag*cosd(ang20-ang0-yang)+V20*V20*ymag*cosd(ang20-ang20-yang+180)
Q20 = V0*V20*ymag*sind(ang20-ang0-yang)+V20*V20*ymag*sind(ang20-ang20-yang+180)


vtemp = D.Vt{1}(nR)/Vbase;
ttemp = D.t{k}(nR);

if 0>(sum(isnan(vtemp))+sum(isinf(vtemp))+sum(isnan(ttemp))+sum(isinf(ttemp)))
    error('NAN or INF in v or t')
end

t = D.t{k}(nR(1)):1/600:D.t{k}(nR(end));
v = interp1(ttemp,vtemp,t,'linear','extrap');
p = interp1(ttemp,D.P{k}(nR)/Sbase,t,'spline','extrap');
q = interp1(ttemp,-D.Q{k}(nR)/Sbase,t,'spline','extrap');
i = interp1(ttemp,-D.I{k}(nR)/Ibase,t,'spline','extrap');

t = t-D.t{k}(nR(1));

mytitle = ['d_centralia_play_' num2str(nR(1)) 'to' num2str(nR(end)) '.mat']

%% !!!!! This is here because angle wrap issue with pst sim!!!!!!!
ang20=ang20-ang0;
v=v*exp(1i*-ang0/180*pi);
ang0=ang0-ang0;

%save(mytitle,'P0',"Q0","f0","V0","ang0","t","v","p","q","i","Q20","P20","V20","ang20");