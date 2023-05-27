%% This code compares the play-in simulation output with the pmu data
% Ron Hruban 2023
close all;
clear all;
clc;


%% Load data
pmu = load("d_centralia_play_1223239to1225021.mat"); % generated from PST simulation
pst = load("playindata_centralia_1223239to1225021.mat"); % playin


%% plot PMU comparison
figure();
subplot(211);
hold on;
plot(pmu.t,-pmu.p,'LineWidth',2,'DisplayName','PMU');
plot(pst.g.t,pst.g.mac.pelect(1,:),'LineWidth',2,'DisplayName','Play-in');
ylabel('P')
hold off;

subplot(212);
hold on;
plot(pmu.t,-pmu.q,'LineWidth',2,'DisplayName','PMU');
plot(pst.g.t,pst.g.mac.qelect(1,:),'LineWidth',2,'DisplayName','Play-in');
legend('Location','best');
ylabel('Q')
xlabel('time (s)')
hold off;


figure();
subplot(211);
hold on;
plot(pmu.t,abs(pmu.v),'LineWidth',2,'DisplayName','PMU');
plot(pst.g.t,abs(pst.g.bus.bus_v(1,:)),'LineWidth',2,'DisplayName','Play-in');
ylabel('V mag (pu)')
hold off;

subplot(212);
hold on;
plot(pmu.t,unwrap(angle(pmu.v)),'LineWidth',2,'DisplayName','PMU');
plot(pst.g.t,unwrap(angle(pst.g.bus.bus_v(1,:))),'LineWidth',2,'DisplayName','Play-in');
legend('Location','best');
xlabel('time (s)')
ylabel('V angle (rad)')
hold off;