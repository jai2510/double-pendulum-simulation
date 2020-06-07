
syms theta1(t) x(t) theta2(t) %  d.o.fs

k = 1;l1 =0.1 ; l2 = 0.5 ; m1 = 0.05 ; m2=0.01 ;g=9.8;

% initial conditions
theta1_0 = pi/18 ; 
x_0 = 0 ;
theta2_0 = pi/9 ; 
theta1_dot_0=0;
x_dot_0=0;
theta2_dot_0=0;
y0 = [theta1_0,x_0,theta2_0,theta1_dot_0,x_dot_0,theta2_dot_0];

eqn1 = diff((m1+m2)*(l1*l1 + x*x)*diff(theta1,1)+2*(m1+m2)*l1*x*diff(theta1,1)+m2*l2*cos(theta1-theta2)*diff(theta2,1)*(l1+x),1) == -m2*l2*(l1+x)*sin(theta1-theta2)*diff(theta1,1)*diff(theta2,1)+m2*l2*diff(theta2,1)*diff(x,1)*cos(theta1-theta2)-(m1+m2)*g*(l1+x)*sin(theta1); 
eqn2 = diff((m1+m2)*diff(x,1)+m2*l2*sin(theta1-theta2)*diff(theta2,1),1) == (m1+m2)*(l1+x)*(diff(theta1,1))^2 + m2*g*cos(theta1) - k*x;
eqn3 = diff(m2*l2*l2*diff(theta2,1)+m2*(l1+x)*l2*diff(theta1,1)*cos(theta1-theta2)+m2*l2*diff(x,1)*sin(theta1-theta2),1) == m1*l2*(l1+x)*sin(theta1-theta2)*diff(theta1,1)*diff(theta2,1)-m2*g*l2*sin(theta2);

vars = [theta1(t),x(t),theta2(t)];

V = odeToVectorField([eqn1,eqn2,eqn3]);

% interval
t0=0;
tf=40;
interval = [t0 tf];

% solutions

M = matlabFunction(V,'vars', {'t','Y'});

ySol = ode45(M,interval,y0);

tValues = linspace(interval(1),interval(2),200);

theta1_values= deval(ySol(1),tValues,1);
x_values = deval(ySol(1),tValues,2);
theta2_values = deval(ySol(1),tValues,3);
theta1_dot_values = deval(ySol(1),tValues,4);
x_dot_values = deval(ySol(1),tValues,5);
theta2_dot_values = deval(ySol(1),tValues,6);

% plotting
%hold on
plot(tValues,theta1_values)
%plot(tValues,x_values)
%plot(tValues,theta2_values)
%plot(tValues,theta1_dot_values)
%plot(tValues,x_dot_values)
%plot(tValues,theta2_dot_values)

syms theta1 x theta2 theta1_dot x_dot theta2_dot;

KE(theta1,x,theta2,theta1_dot,x_dot,theta2_dot) = 0.5*(m1)*((diff(((l1)+x)*sin(theta1),t))^2 + (diff(((l1)+x)*cos(theta1),t))^2) + 0.5*(m2)*((diff(((l1)+x)*sin(theta1)+(l2)*sin(theta2),t))^2 + (diff(((l1)+x)*cos(theta1)+(l2)*cos(theta2),t)))^2;
PE(theta1,x,theta2,theta1_dot,x_dot,theta2_dot) = -(m1)*g*((l1)+x)*cos(theta1) - (m2)*g*(((l1)+x)*cos(theta1)+l2*cos(theta2)) + 0.5*k*x*x ;

a=subs(KE+PE,{theta1 x theta2 theta1_dot x_dot theta2_dot},{theta1_values x_values theta2_values theta1_dot_values x_dot_values theta2_dot_values});

%plot(tValues,a)
axis([0 40 -0.5 0.5])
hold off
