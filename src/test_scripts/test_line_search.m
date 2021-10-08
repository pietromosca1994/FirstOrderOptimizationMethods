clear all;
clc;

project_path='C:\Users\PIM\Desktop\Documents\GitHub_Repo\Conjugate_Gradient';

% %function to evaluate
% f=@(x)(atan(x));
% %a lower bracket  
% a=0.5*pi;
% %upper bracket
% b=0.7*pi;

% % function to evaluate
% f=@(x)((x+3).*(x-1).^2);
% % a lower bracket  
% a=-4;
% % upper bracket
% b=4/3;

% function to evaluate
f=@(x)(x-4);
% a lower bracket  
a=-4;
% upper bracket
b=4/3;


x=linspace(-2*pi,2*pi,500);
y=f(x);

[root, history]=LineSearchBrent(f,a,b,0);
disp(['root: ', num2str(root)]);
disp(['steps: ', num2str(history.steps)]);
disp(['feval: ', num2str(history.feval)]);

figure()
subplot(1,3,1)
plot(x,y)
xlabel('x');
ylabel('y');
hold on 
scatter(history.b, f(history.b), 'r');

subplot(1,3,2);
plot(history.a);
xlabel('n');
ylabel('a');

subplot(1,3,3);
plot(history.b);
xlabel('n');
ylabel('b');

