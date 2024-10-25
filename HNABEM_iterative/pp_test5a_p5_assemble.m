% reading in data and rearranging, p = 5

clear all
addpath('test5a_p5_7_HNA_pi_4/')

load('test5a_HNA_pmax5_overlap2_quaddof_5_10')

%relabelling
G1_data_HNA_5_10 = G1_data_HNA;
G2_data_HNA_5_10 = G2_data_HNA;
phi1_HNA_5_10 = phi1_HNA;
phi2_HNA_5_10 = phi2_HNA;

% clearing up
clear G1_data_HNA G2_data_HNA phi1_HNA phi2_HNA

%loading second load of data
load('test5a_HNA_pmax5_overlap2_quaddof_20_80.mat')
G1_data_HNA_20_80 = G1_data_HNA;
G2_data_HNA_20_80 = G2_data_HNA;
phi1_HNA_20_80 = phi1_HNA;
phi2_HNA_20_80 = phi2_HNA;

% clearing up
clear G1_data_HNA G2_data_HNA phi1_HNA phi2_HNA

%relabel
G1_data_HNA = {}; 
G1_data_HNA = {G1_data_HNA_5_10{1},...
    G1_data_HNA_5_10{2}, G1_data_HNA_20_80{1},...
    G1_data_HNA_20_80{2}, G1_data_HNA_20_80{3}}

G2_data_HNA = {}; 
G2_data_HNA = {G2_data_HNA_5_10{1}, ...
    G2_data_HNA_5_10{2}, G2_data_HNA_20_80{1},...
    G2_data_HNA_20_80{2}, G2_data_HNA_20_80{3}}

phi1_HNA = {}; 
phi1_HNA = {phi1_HNA_5_10{1}, phi1_HNA_5_10{2}, ...
    phi1_HNA_20_80{1}, phi1_HNA_20_80{2}, phi1_HNA_20_80{3}}

phi2_HNA = {}; 
phi2_HNA = {phi2_HNA_5_10{1}, phi2_HNA_5_10{2}, ...
    phi2_HNA_20_80{1}, phi2_HNA_20_80{2}, phi2_HNA_20_80{3}}

 save('test5a_HNA_pmax5_overlap2_dof5_80', 'G1_data_HNA', ...
     'G2_data_HNA', 'phi1_HNA', 'phi2_HNA')


