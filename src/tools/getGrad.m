%% Author: Pietro Mosca
%% Email: pietromosca1994@gmail.com
%% Date: 22.09.2020

%% Description:
% Function which computes numerically the gradient of a n-dimensional
% function

%% Function Arguments
% f: function handle
% x0: target point

%% Parameters
% epsilon: small number

function [grad, history]=getGrad(f, x0)
    %% algorithm parameters
    epsilon=10^-4;

    %% variables initialization
    grad=zeros(size(x0));
    history.feval=0;

    %% gradient computation
    for i=1:length(grad)
        x_hi=x0;
        x_lo=x0;
        x_hi(i)=x_hi(i)+epsilon/2;
        x_lo(i)=x_lo(i)-epsilon/2;

        grad(i)=(feval(f, x_hi)-feval(f, x_lo))/epsilon;
        history.feval=history.feval+2;
    end
end