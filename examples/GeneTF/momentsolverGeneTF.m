function [average,var,outm] = momentsolverGeneTF(Tmax)

%specifies the degree of the closure, i.e. the size of the approximate moment system
maxdeg = 2;

%specifies the closure type
closureMethod = 'dm';

%reads the specified model
net=readNet('GeneTF.net');

%computes both the open moment system and the closed approximation
%also creates a function funname that can be passed to an ode solver
mdyn=closureDynamics(net,maxdeg,closureMethod,'funname');

%time horizon
Ts = 1;
t = 0:Ts:Tmax;

%sets the intial condition of the ode 
x0=mdyn.Mu.x0;

tic;

options = odeset('RelTol',1e-13);
%uses ode23s to solve the approximate system
[t,mu]=ode23s(@(t,x)funname(x),t,x0,options);

%converts uncentered moments to centered moments
[average,stddev] = getCMoments(net,mdyn,mu,{'X';'Y';'Z'});

outm = toc;

vari = stddev.^2;

color = 'g';

figure(1)
hold on
plot(t,average(:,1),color)

figure(2)
hold on
plot(t,average(:,2),color)

figure(3)
hold on
plot(t,average(:,3),color)

figure(4)
hold on
plot(t,vari(:,1),color)

figure(5)
hold on
plot(t,vari(:,2),color)

figure(6)
hold on
plot(t,vari(:,3),color)

% 
% [handles,average,stddev]=plotCMoments(net,mdyn,t,mu,symExpression_1 ,style_1,symExpression_2 ,style_2,...)

end

