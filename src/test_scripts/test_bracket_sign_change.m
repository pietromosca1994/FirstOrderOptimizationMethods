clear all;
clc;

project_path='C:\Users\PIM\Desktop\Documents\GitHub_Repo\Conjugate_Gradient';

% function to evaluate
f=@(x)(atan(x));
% a lower bracket  
a=0.5*pi;
% upper bracket
b=0.7*pi;
% k=2 usual setting for expansion factor k
k=2;

x=linspace(-2*pi,2*pi,500);
y=f(x);

[a, b, history]=BracketSignChange(f,a,b,k);

figure()
subplot(1,3,1)
plot(x,y)
xlabel('x');
ylabel('y');
hold on 
scatter(history.a, f(history.a), 'r');
scatter(history.b, f(history.b), 'b');

subplot(1,3,2);
plot(history.a);
xlabel('n');
ylabel('a');

subplot(1,3,3);
plot(history.b);
xlabel('n');
ylabel('b');

disp(['a: ', num2str(a)]);
disp(['b: ', num2str(b)]);
disp(['steps: ', num2str(history.steps)]);

