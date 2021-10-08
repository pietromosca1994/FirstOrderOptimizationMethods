%% Author: Pietro Mosca
%% Email: pietromosca1994@gmail.com
%% Date: 22.09.2020

%% Description:
% function which brackets the target function f in a interval in which the
% target function has values of opposite sign. If the condition f(a)*f(b)
% is met, the existence of a root in the interval is guaranteed.

%% Function Arguments
% f: target function handler
% a: low bracket (initial guess)
% b: high bracket (initial guess)
% k: expansion factor

%% Parameters
% maxiter: maximum number of iterations

function [a,b,fa,fb,history]=BracketSignChange(f, a, b, k, maxiter)
    
    %% algorithm
    % check that a<b and in case swap
    if a>b
        % swap a and b
        temp=a;
        a=b;
        b=temp;
    end
    
    % log history
    i=1;
    history.a(i)=a;
    history.b(i)=b;
    history.feval=0;
    
    % initilise center and half width
    center=(a+b)/2;
    half_width=(b-a)/2;

    % repeat untile f(a) and f(b) have different signs
    for i=1:maxiter
        fa=f(a);
        fb=f(b);
        history.feval=history.feval+2;
        
        if fa*fb>=0
            half_width=half_width*k;
            a=center-half_width;
            b=center+half_width;
            
            % log history
            history.a(i)=a;
            history.b(i)=b;
            history.fa(i)=fa;
            history.fb(i)=fb;
        else
            break;       
        end
    end
    
    if not(fa*fb<0)
        error('Maximum number of iterations exceeded');
    end
    
    history.steps=i;
end
