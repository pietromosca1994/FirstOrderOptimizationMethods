%% Adam Optimizer Implementation
%% Author:  Pietro Mosca
%% Email:   pietromosca1994@gmail.com
%% Date:    20.09.2020
%% Reference: https://arxiv.org/pdf/1412.6980.pdf

function [x_min, history]=Adam(f, x0)
  %% algorithms parameters initialization 
  % a: stepsize
  params.a=0.01;
  % b1, b2 exponential decay rates
  params.b1=0.9;
  params.b2=0.999;
  % epsilon: convergence threshold
  params.epsilon=10^-8;
  
  %% variables initialization
  x_min=x0;
  % f_min: function minimum
  f_min=f(x_min);
  % m: first momentum
  m=zeros(size(x0));
  % m_exp: expected first momentum
  m_exp=0;
  % v: second momentum;
  v=zeros(size(x0));
  % v_exp: expected second momentum
  v_exp=0;
  % grad: gradient
  grad=zeros(size(x0));
  % t: step
  t=0;
  
  %% optimization
  history.x_min(1,:)=x_min;
  history.f_min(1)=f_min;
  err=1;
  
  while err>params.epsilon && t<10000
    t=t+1;
    % get gradients w.r.t. stochastic objective at timestep t
    grad=getGrad(f, x_min);
    % update biased first moment estimate
    m=params.b1*m+(1-params.b1)*grad;
    % update biased second raw moment estimate
    v=params.b2*v+(1-params.b2)*grad.^2;
    % compute bias-corrected first moment estimate
    m_exp=m./(1-params.b1^t);
    % compute bias-corrected second raw moment estimate
    v_exp=v./(1-params.b2^t);
    % update parameters
    x_min=x_min-params.a*m_exp/(sqrt(v_exp)+params.epsilon);
    % update for gradient descent
    % x_min=x_min-params.a*grad;
    
    % compute f_min
    f_min=f(x_min);
    % log history
    history.x_min(t+1,:)=x_min;
    history.f_min(t+1)=f_min;
    history.t=t;
    err=abs(history.f_min(t+1)-history.f_min(t));
   end
  
  disp(['x_min: ', num2str(x_min)]);
  disp(['f_min: ', num2str(f_min)]);
  disp(['t: ', num2str(t)]);
 
end
