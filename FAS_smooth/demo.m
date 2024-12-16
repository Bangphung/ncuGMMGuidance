clear; clc; close all;

% input.mat has input signal (x) and sampling interval (dt)
load input.mat;

%% EXAMPLE 1 (use debug option is 'True' and provide figure name and path):
[hp_freq, lp_freq] = cornerFreqs(x,dt,'plot_name','Demo','plot_path','./figures/','debug','True');

%% EXAMPLE 2 (without plotting):
%[hp_freq, lp_freq] = cornerFreqs(x,dt);