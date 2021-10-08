%% Booth function implementation
%% Author:  Pietro Mosca
%% Email:   pietromosca1994@gmail.com
%% Date:    20.09.2020
%% Reference: https://en.wikipedia.org/wiki/Test_functions_for_optimization

function [y]=booth(x)
    y=(x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
end
