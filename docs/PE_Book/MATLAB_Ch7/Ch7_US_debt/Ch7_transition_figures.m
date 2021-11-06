%   Ch7_transition_figures.m
%
%   presents the figure for the transition in Chapter 7.4 in
%   Heer(2018), Public Economics - The Macroeconomic Perspective
%
%   first run Ch7_debt_transition.m for policy=0,1,2,3
%
%

load('Ch7_US_debt_transition.mat');

    nt1=30;
    periodst = linspace(2010,2010+5*nt1-5,nt1);
    periodst = periodst';

    figure
    subplot(3,2,1);
    plot(periodst,kt_policy0(1:nt1),periodst,kt_policy1(1:nt1),periodst,kt_policy2(1:nt1),periodst,kt_policy3(1:nt1));
    xlabel('year');
	title('capital stock');
    subplot(3,2,2);
    plot(periodst,nt_policy0(1:nt1),periodst,nt_policy1(1:nt1),periodst,nt_policy2(1:nt1),periodst,nt_policy3(1:nt1));
	title('labor');
    xlabel('year');
    subplot(3,2,3);
    plot(periodst,taunt_policy0(1:nt1)-taupt_policy0(1:nt1),periodst,taunt_policy1(1:nt1)-taupt_policy1(1:nt1),periodst,taunt_policy2(1:nt1)-taupt_policy2(1:nt1),periodst,taunt_policy3(1:nt1)-taupt_policy3(1:nt1));
	title('\tau^n');
    xlabel('year');
    subplot(3,2,4);
    plot(periodst,taupt_policy0(1:nt1),periodst,taupt_policy1(1:nt1),periodst,taupt_policy2(1:nt1),periodst,taupt_policy3(1:nt1));
	title('\tau^p');
    xlabel('year');
    subplot(3,2,5);
    plot(periodst,trt_policy0(1:nt1),periodst,trt_policy1(1:nt1),periodst,trt_policy2(1:nt1),periodst,trt_policy3(1:nt1));
	title('government transfers');
    xlabel('year');
    subplot(3,2,6);
    plot(periodst,bt_policy0(1:nt1),periodst,bt_policy1(1:nt1),periodst,bt_policy2(1:nt1),periodst,bt_policy3(1:nt1));
	title('government debt');
    xlabel('year');
    legend('Policy 1','Policy 2','Policy 3','Policy 4');