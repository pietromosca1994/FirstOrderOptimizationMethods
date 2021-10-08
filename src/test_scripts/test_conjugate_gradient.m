clear all;
clc;

project_path='C:\Users\PIM\Desktop\Documents\GitHub_Repo\Conjugate_Gradient';

run('init_lib_FirstOrderOptimizationMethods.m');

% algorithm parameters definition
method='PR';
% 'FR' for Fletcher-Reeves coefficient
% 'PR' for Polak-Ribier coefficient 
% 'HS' for Hestenes-Stiefel coefficient
% 'DY' for Dai-Yuan coefficient

ytol=10^-1;
maxiter=10000;

opt_bound=[-2, 2;  % dimension 1 
           -2, 2]; % dimension 2
       
f=@rosenbrock;

% optimization
x0=opt_bound(1,:)+rand*(opt_bound(2,:)-opt_bound(1,:)); % starting point initilization
[x_min, history]=ConjugateGradient(f, x0, 'PR', ytol, maxiter);

% results display 
disp(['x_min: ', num2str(x_min)])
disp(['f_min :' , num2str(history.f(end))]);
disp(['steps :', num2str(history.steps)]);
disp(['feval :', num2str(history.feval)]);

% plotting
res=100;
x1=linspace(min(min(history.x(:,1)), opt_bound(1,1)), max(max(history.x(:,1)),opt_bound(1,2)), res);
x2=linspace(min(min(history.x(:,2)), opt_bound(2,1)), max(max(history.x(:,2)), opt_bound(2,2)), res);
y=zeros(numel(x1), numel(x2));

for i=1:numel(x1)
  for j=1:numel(x2)
    y(i,j)=f([x1(i), x2(j)]);
  end
end

figure()
subplot(1,2,1)
colormap ("default");
surf(x1, x2, y')
shading interp;
xlabel('x1');
ylabel('x2');
zlabel('y');
hold on;
plot3(history.x(:,1), history.x(:,2), history.f, '-o', 'MarkerFaceColor', 'r', 'Color','r');
%scatter3(history.x(:,1), history.x(:,2), history.f, 'r');
set(gcf,'color','w');

subplot(1,2,2)
plot(history.f);


