% reading in data and rearranging p = 8

clear all

load('test5a_HNA_pmax8_overlap2_quaddof_5_20')

%relabelling
G1_data_HNA_5_20 = G1_data_HNA;
G2_data_HNA_5_20 = G2_data_HNA;
phi1_HNA_5_20 = phi1_HNA;
phi2_HNA_5_20 = phi2_HNA;

% clearing up
clear G1_data_HNA G2_data_HNA phi1_HNA phi2_HNA

%loading second load of data
load('test5a_HNA_pmax8_overlap2_quaddof_40_80.mat')
G1_data_HNA_40_80 = G1_data_HNA;
G2_data_HNA_40_80 = G2_data_HNA;
phi1_HNA_40_80 = phi1_HNA;
phi2_HNA_40_80 = phi2_HNA;

% clearing up
clear G1_data_HNA G2_data_HNA phi1_HNA phi2_HNA

%relabel
G1_data_HNA = {}; 
G1_data_HNA = {G1_data_HNA_5_20{1},...
    G1_data_HNA_5_20{2}, G1_data_HNA_5_20{3},...
    G1_data_HNA_40_80{1}, G1_data_HNA_40_80{2}}

G2_data_HNA = {}; 
G2_data_HNA = {G2_data_HNA_5_20{1}, ...
    G2_data_HNA_5_20{2}, G2_data_HNA_5_20{3},...
    G2_data_HNA_40_80{1}, G2_data_HNA_40_80{2}}

phi1_HNA = {}; 
phi1_HNA = {phi1_HNA_5_20{1}, phi1_HNA_5_20{2}, ...
    phi1_HNA_5_20{3}, phi1_HNA_40_80{1}, phi1_HNA_40_80{2}}

phi2_HNA = {}; 
phi2_HNA = {phi2_HNA_5_20{1}, phi2_HNA_5_20{2}, ...
    phi2_HNA_5_20{3}, phi2_HNA_40_80{1}, phi2_HNA_40_80{2}}

 save('test5a_HNA_pmax8_overlap2_dof5_80', 'G1_data_HNA', ...
     'G2_data_HNA', 'phi1_HNA', 'phi2_HNA')


