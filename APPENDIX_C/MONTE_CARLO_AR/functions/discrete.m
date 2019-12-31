function s=discrete(pst)

% PURPOSE: This function draws a discrete random number s from the 
% vector of probabilities pst.
%******************************************************************
% USAGE: s=discrete(pst)
% *****************************************************************
% INPUT:
% pst = vector of states probabilities
% OUTPUT:
% s = state draw
%******************************************************************
% Written by
% Davide Pettenuzzo
% Phd Boconi University - Milan
% Italy
% email: davide.pettenuzzo@uni-bocconi.it

states=(1:rows(pst))';
lb=lag(cumsum(pst));  %lower bound
ub=cumsum(pst);       %upper bound

%check condition
u=rand(1,1);
s=states .* (u > lb & u <= cumsum(pst));
s=selif(s,s>0);