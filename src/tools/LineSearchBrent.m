%% Author: Pietro Mosca
%% Email: pietromosca1994@gmail.com
%% Date: 22.09.2020

%% Description:
% function which find's the root of a target function combining the
% bisection method (reliable bracketing method), the secant method and 
% inverse quadratic interpolation (fast open methods) accordingly to 
% Brent's Method (Dekker's method 1969)

%% Function Arguments
% f: target function handler
% a: low bracket (initial guess)
% b: high bracket (initial guess)

%% Parameters
% ytol: y tollerance to stop the search
% maxiter: maximum number of iterations

function [b, history]=LineSearchBrent(f, a, b, verbose)
    
    %% algorithm paramters initialization 
    ytol=2*10^-3;
    maxiter=1000;
    
    %% algorithm
    % bracket the function ensuring f(a)*f(b)<0
    [a,b,fa,fb,bsc_history]=BracketSignChange(f,a,b,2,10000);
    history.feval=bsc_history.feval; % increment function evaluation counter
    
    % b is the current guess for the root of f ( |f(b)| is the closest point
    % to 0)
    if abs(fa)<abs(fb)
        % swap a b to have |f(a)|>|f(b)|
        temp=a;
        a=b;
        b=temp;
        
        temp=fa;
        fa=fb;
        fb=temp;
    end
   
    % end of the interval that contains the sign change
    c=a;
    fc=fa;
    
    % bisection flag
    mflag=1;
    % cache
    d=c;
    i=1;
    
    % log root
    history.a(i)=a;
    history.b(i)=b;   
    
    % repeat until convergence
    while abs(fb)>ytol && i<maxiter % convergence condition
        % update step number
        i=i+1;

        if abs(fa-fc)>ytol && (fb-fc)>ytol
            % inverse quadratic interpolation
            % possibility to speed up the lagrange polynomial using RST
            % implementation https://mathworld.wolfram.com/BrentsMethod.html
            s=(a*fb*fc)/((fa-fb)*(fa-fc))+(b*fa*fc)/((fb-fa)*(fb-fc))+(c*fa*fb)/((fc-fa)*(fc-fb));
        else
            % secant method
            s=b-fb*(b-a)/(fb-fa);
        end
        
        % Brent's conditions
        delta=abs(2*eps*abs(fb));
        min1=abs(s-b);
        min2=abs(b-c);
        min3=abs(c-d);
        
        if (s<(3*a+b)/4 && s>b) || ... %|| (s>b && s<(3*a+b)) || ... & s is not between (3*a+b)/4 and b
           mflag==1 && min1>=min2/2 || ...  
           mflag==0 && min1>=min3/2 || ...
           mflag==1 && min2<delta || ...
           mflag==0 && min3<delta
       
            % bisection method
            s=(a+b)/2;
            mflag=1;
        else
            mflag=0;
        end
        
        fs=f(s);
        history.feval=history.feval+1;
                
        d=c;
        c=b;        
        
        if fa*fs<0
            b=s;
            fb=fs;
        else
            a=s;
            fa=fs;
        end
        
        if abs(fa)<abs(fb)
        % swap a b to have |f(a)|>|f(b)|
            temp=a;
            a=b;
            b=temp;
            
            temp=fa;
            fa=fb;
            fb=temp;
        end
        
        % log history
        history.b(i)=b;
        history.a(i)=a;
      
        % verbose mode
        if verbose==1
            disp(['i:' num2str(i)]);
            disp(['a:' num2str(a)]);
            disp(['b:' num2str(b)]);
            disp(['s:' num2str(s)]);
            disp(['fa:' num2str(fa)]);
            disp(['fb:' num2str(fb)]);
            disp(['fs:' num2str(fs)]);           
            disp(['c:' num2str(c)]);
            disp('_________________________________');
        end
    
    end
    
    if i==maxiter
        error('Reached maximum number of iterations')
    end
    
    history.steps=i; 
end