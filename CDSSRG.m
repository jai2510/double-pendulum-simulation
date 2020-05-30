%% Analytical Mechanics

syms theta1(t) x(t) theta2(t) %  d.o.fs
%syms L(t)

theta1_dot = diff(theta1,t);
x_dot = diff(x,t);
theta2_dot = diff(theta2,t);

% Parameters
k = 1;l1 =0.1 ; l2 = 0.5 ; m1 = 0.05 ; m2=0.01 ;g=9.8;

% initial conditions
theta1_0 = pi/18 ; 
x_0 = 0 ;
theta2_0 = pi/9 ; 
theta1_dot_0=0;
x_dot_0=0;
theta2_dot_0=0;
 y0 = [theta1_0,x_0,theta2_0,theta1_dot_0,x_dot_0,theta2_dot_0];
%y0 = [theta1_0,x_0,theta2_0]%theta1_dot_0,x_dot_0,theta2_dot_0];

%
KE(t) = 0.5*m1*((l1*theta1_dot)^2+(x_dot)^2+(x*theta1_dot)^2+2*l1*x*(theta1_dot^2)) + 0.5*m2*((l1*theta1_dot)^2+(x_dot)^2+(x*theta1_dot)^2+(l2*theta1_dot)^2+2*(l1*x*(theta1_dot)^2+l1*l2*theta1_dot*theta2_dot*cos(theta1-theta2)+x_dot*l2*theta2_dot*sin(theta1-theta2)+x*l2*theta1_dot*theta2_dot*cos(theta1-theta2)));

PE(t) = 0.5*k*x*x - m1*g*(l1+x)*cos(theta1) - m2*g*(l2*cos(theta2) + (l1+x)*cos(theta1));

L = KE - PE;
%
f1=diff(L,t)./diff(theta1_dot,t);
f2=diff(L,t)./diff(x_dot,t);
f3=diff(L,t)./diff(theta2_dot,t);
%
% Equations
eqn1 = diff(f1,t) == diff(KE-PE,t)./diff(theta1,t);
eqn2 = diff(f2,t) == diff(KE-PE,t)./diff(x,t);
eqn3 = diff(f3,t) == diff(KE-PE,t)./diff(theta2,t);
%
vars = [theta1(t),x(t),theta2(t)];

V = odeToVectorField([eqn1,eqn2,eqn3]);

M = matlabFunction(V,'vars', {'t','Y'});

% interval
t0=0;
tf=100;
interval = [t0 tf];

% solutions

ySol = ode45(M,interval,y0);

tValues = linspace(interval(1),interval(2),10000);
yValues = deval(ySol,tValues,1);

% plotting
plot(tValues,yValues)





