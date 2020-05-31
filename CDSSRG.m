syms theta1(t) x(t) theta2(t) %  d.o.fs

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


eqn1 = diff(theta1,2)==(-m1*g*((l1)+x)*sin(theta1)-m2*g*((l1)+x)*sin(theta1)-2*(m1)*x*(theta1_dot)*x_dot-2*(m1)*(l1)*x_dot*theta1_dot-2*(m2)*x*theta1_dot*x_dot-(m2)*(l1)*(l2)*(diff(theta2_dot,1)*cos(theta1-theta2)+(theta2_dot^2)*sin(theta1-theta2))-(m2)*(l2)*(x*cos(theta1-theta2)*diff(theta2,2)+x*sin(theta1-theta2)*(theta2_dot)^2))./((m1)*(l1)*(l1)+(m1)*x*x+2*(m1)*(l1)*(x)+(m2)*(l1)*(l1)+(m2)*x*x+2*(m2)*(l1)*x);
eqn2 = diff(x,2)==(m1*(x+l1)*theta1_dot^2 +m2*(x+l1)*theta1_dot^2 + m2*l2*theta1_dot*theta2_dot*cos(theta1-theta2)-k*x+(m1+m2)*g*cos(theta1)-m2*l2*(diff(theta2_dot,1)*sin(theta1-theta2)+theta2_dot*(theta1_dot-theta2_dot)*cos(theta1-theta2)))/(m1+m2);
eqn3 = diff(theta2,2)==(-m2*l2*x_dot*theta2_dot*cos(theta1-theta2)-m2*l1*l2*(diff(theta1,2)*cos(theta1-theta2)-sin(theta1-theta2)*theta1_dot^2)-m2*l2*(diff(x_dot,1)*sin(theta1-theta2)+x_dot*cos(theta1-theta2)*(theta1_dot-theta2_dot))-m2*l2*(x_dot*theta1_dot*cos(theta1-theta2)+x*diff(theta1,2)*cos(theta1-theta2)-x*theta1_dot^2*sin(theta1-theta2)))./(m2*l2*l2);

vars = [theta1(t),x(t),theta2(t)];

V = odeToVectorField([eqn1,eqn2,eqn3]);


% interval
t0=0;
tf=2000;
interval = [t0 tf];

% solutions

M = matlabFunction(V,'vars', {'t','Y'});

ySol = ode45(M,interval,y0);

tValues = linspace(interval(1),interval(2),200);
yValues1 = deval(ySol(1),tValues,1);
yValues2 = deval(ySol(1),tValues,2);
yValues3 = deval(ySol(1),tValues,3);
yValues4 = deval(ySol(1),tValues,4);
yValues5 = deval(ySol(1),tValues,5);
yValues6 = deval(ySol(1),tValues,6);

% uncomment whichever ones you want to plot, from below.
% plotting
plot(tValues,yValues1)
%plot(tValues,yValues2)
%plot(tValues,yValues3)
