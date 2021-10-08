%% Author: Pietro Mosca
%% Email: pietromosca1994@gmail.com
%% Date: 22.09.2020

%% Description:
% Function which implements gradient descent using different paraemters
% from different authors
% - Fletcher-Reeves
% - Polak-Ribier
% - Hestenes-Stiefel
% - Dai-Yuan

%% Function Arguments
% f: target function handler
% x0: starting point
% b= 'FR' for Fletcher-Reeves coefficient
% b= 'PR' for Polak-Ribier coefficient 
% b= 'HS' for Hestenes-Stiefel coefficient
% b= 'DY' for Dai-Yuan coefficient

%% Parameters
% ytol: y tollerance for convergence
% maxiter: maximum number of iterations for the algorithm to stop

function [x, history]=ConjugateGradient(f, x0, method, ytol, maxiter)
    %% algorithm parameters initialization
    ytol=10^-6;
    maxiter=10000;
    history.feval=0;
    
    %% variables initialization
    i=1;
    
    %% algorithm
    history.x(i, :)=x0;
    history.f(i)=f(x0);
    
    % subscript 0 indicates cached values
    % initialization
    [grad0, g_history]=getGrad(f, x0);
    % descent direction
    grad0=-grad0;
    history.feval=history.feval+g_history.feval;
    % s initialization
    s0=grad0;
    % argmin(f(x+alpha*s)) is the same as s*grad(f(x+alpha*s))=0
    falpha=@(alpha)(x0*getGrad(f,x0+alpha*grad0)');
    % line search
    [alpha0, ls_history]=LineSearchBrent(falpha, 0, 1, 0);
    history.feval=ls_history.feval*2; % univariate gradient evaluation (2 function evaluations)
    % update position
    x=x0+alpha0*s0;
    
    fx0=f(x0);
    fx=f(x);
    e=abs(fx-fx0);
    
    x0=x;
    
    % convergence criteria
    while e>=ytol && i<maxiter
        i=i+1;
        % compute the steepes direction
        [grad, g_history]=getGrad(f, x);
        % descent direction
        grad=-grad;
        history.feval=history.feval+g_history.feval;
        
        % compute b accordingly to the authors
        if strcmp(method, 'FR') %
            b=(grad*grad')/(grad0*grad0');
        elseif strcmp(method, 'PR')
            b=(grad*(grad-grad0)')/(grad0*grad0');
        elseif strcmp(method, 'HS')
            b=(grad*(grad-grad0)')/(-s0*(grad-grad0)');
        elseif strcmp(method, 'DY')
            b=(grad*grad')/(-s0*(grad-grad0)');
        end
        
        % update the conjugate direction
        s=grad+b*s0;
        
        % perform a line search optimize alpha (stepsize)
        falpha=@(alpha)(s*getGrad(f,x+alpha*s)');
        [alpha, ls_history]=LineSearchBrent(falpha,0,alpha0,0);
        history.feval=history.feval+ls_history.feval*2;
        
        % update position
        x=x0+alpha*s;
        
        fx0=f(x0);
        fx=f(x);
        
        e=abs(fx0-fx);
        
        % update parameters
        s0=s;
        grad0=grad;
        x0=x;
        
        % log
        history.x(i, :)=x;
        history.f(i)=fx;
        
    end
    
    if i==maxiter
        error('Reached maximum number of iterations')
    end
    
    history.steps=i;
end
