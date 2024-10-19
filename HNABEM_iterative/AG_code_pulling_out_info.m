function [v_N1, GOA1, colMatrix1, colRHS1, col_points1,...
v_N2, GOA2, colMatrix2, colRHS2, col_points2, VHNA1, VHNA2, Xstruct1, Xstruct2] ...
    = AG_code_pulling_out_info(pMax, cL, sigmaGrad, nLayers, OverSample, ...
    Gamma1, Gamma2, kwave, uinc )
% In this function we are computing all the things we need to from AG code
% in one place.
% Inputs:
%
% Outputs:
%
%

% construct the HNA basis (single mesh):
VHNA1 = HNAoverlappingMesh(Gamma1, pMax, kwave, nLayers, sigmaGrad);
% DOFs1 = length(VHNA1.el); %get total #DOFs

VHNA2 = HNAoverlappingMesh(Gamma2, pMax, kwave, nLayers, sigmaGrad);
% DOFs2 = length(VHNA2.el); %get total #DOFs

% construct the single layer potential 'operator' ---------------------------
S1=singleLayer(kwave,Gamma1);
S2=singleLayer(kwave,Gamma2);

%solve (and time)
[v_N1, GOA1, colMatrix1, colRHS1, T1, Xstruct1, OpPsiX1, fX1] = ColHNA(S1, VHNA1, uinc, Gamma1,'oversample', OverSample);%, 'progress');
[v_N2, GOA2, colMatrix2, colRHS2, T2, Xstruct2, OpPsiX2, fX2] = ColHNA(S2, VHNA2, uinc, Gamma2,'oversample', OverSample);%, 'progress');

for j = 1:length(colRHS1)
    col_points1(j, 1) = Xstruct1(1, j).x;
    
end

for j = 1:length(colRHS2)
    col_points2(j, 1) = Xstruct2(1, j).x;
end

