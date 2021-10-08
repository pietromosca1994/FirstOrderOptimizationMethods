%% Rosenbrock function implementation
%% Author:  Pietro Mosca
%% Email:   pietromosca1994@gmail.com
%% Date:    20.09.2020
%% Reference: https://en.wikipedia.org/wiki/Rosenbrock_function

function [y]=rosenbrock(x)
    d=length(x);
    y=0;
    for i=1:(d-1)
        xi=x(i);
        xnext=x(i+1);
        new=100*(xnext-xi^2)^2+(xi-1)^2;
        y=y+new;
    end
  %y=sum(100*(x(2:end)-x(1:end-1).^2).^2+(1-x(1:end-1).^2));
  %y=(1-x(1))^2+2*(x(2)-x(1)^2)^2;
  %y=sum(100*(x(2:end)-x(1:end-1).^2).^2+(x(1,end-1)-1).^2);
end