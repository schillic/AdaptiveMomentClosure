function pdf=getDistribution(net,mdyn,mu,symExpression,x,distribution,symBoolean)
% pdf=getDistribution(net,mdyn,mu,symExpression,support,type)
% 
% Given a vector mu of uncentered moments, this function estimates the
% distribution of a given expression involving the populations of the
% different species. Internally, this function uses getCMoments() to
% compute means and standard deviations. 
% 
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% mdyn is a structure characterizing the moment dynamics, typically
%      computed using closureDynamics(), momentDynamics(), or momentClosure(). 
%      Details can be obtained with 'help mdyn'
%
% mu : Time series of uncentered moments with one time instant per row
%      and one moment per columns. Such vectors are typically obtained
%      from an ode solver that computes the solution to some
%      approximate moment dynamics. The meaning of each column is
%      specified by mdyn.Mu
%
% mu is a vector containing the value of Mu or [Mu;Phi] at a given
%    time instant, typically computed using ode23s().
%
% symExpression is a symbolic expression whose distribution will be estimated
%
%      This function uses getCMoments() to compute the mean and
%      standard deviation from which the distribution is
%      estimated. The reader is referred to the documentation of
%      getCMoments() for details on the computation of the mean and
%      standard deviation.   
%
% x is a vector containing the points at which the distribution will
%   be computed.
%
% distribution is a character string specifying which type of
%              distribution should be assumed. This parameter can take
%              the following values:
%   . 'normal' when a normal distribution should be assumed.
%   . 'lognormal' when a lognormal distribution should be assumed.
%   . 'binomial' when a binomial distribution should be assumed.
%   . 'mix_lognormal' when it should be assumed that the distribution
%              is lognormal, when conditioned to any one of the two
%              values for a given boolean-valued expression
%              'symBoolean'. This will result in a convex combination
%              (mixture) of two lognormal distributions. 
%   . 'mix_normal' when it should be assumed that the distribution
%              is lognormal, when conditioned to any one of the two
%              values for a given boolean-valued expression
%              'symBoolean'. This will result in a convex combination
%              (mixture) of two normal distributions. 
%
% symBoolean (optional) is a boolean-valued symbolic expression needed
%            for any of the 'mixture' distributions discussed above.
%
% Output
% ------
%
% pdf is a vector containing the values of the probability density
%     function (pdf) at the points in x. 
%
% Copyright (C) 2006  Joao Hespanha

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%
% Modification history:
%
% Created May 8, 2007 
%
% January 10, 2007
% updated inline documentation

if prod(size(mu))~=max(size(mu))
  error('getDistribution: uncentered input should be a column vector and not a [%d,%d] matrix\n', ...
	size(mu,1),size(mu,2));
end

% convert to row vector
mu=mu(:)'; 
  
switch lower(distribution)
 case 'normal',
  [average,stddev]=getCMoments(net,mdyn,mu,symExpression);
  pdf=1/stddev/sqrt(2*pi)*exp(-(x-average).^2/2/stddev^2);
 
 case 'lognormal',
  [average,stddev]=getCMoments(net,mdyn,mu,symExpression);
  si=sqrt(log(1+stddev^2/average^2));
  mu=log(average)-si/2;
  pdf=exp(-((log(x)-mu).^2)/2/si^2)./x/si/sqrt(2*pi);

 case 'binomial',
  [average,stddev]=getCMoments(net,mdyn,mu,symExpression);
  p=1-stddev^2/average;
  n=average/p;
  pdf=gamma(n+1)./gamma(x+1)./gamma(n-x+1).*p.^x.*(1-p).^(n-x);

 case 'mix_normal',
  symExpression=sym(symExpression);
  symBoolean=sym(symBoolean);

  p=getCMoments(net,mdyn,mu,symBoolean);
  [average1,stddev1]=getCMoments(net,mdyn,mu,symExpression*symBoolean);
  average=average1/p;stddev=stddev1/p;
  pdf1=1/stddev/sqrt(2*pi)*exp(-(x-average).^2/2/stddev^2);
  [average0,stddev0]=getCMoments(net,mdyn,mu,symExpression*(1-symBoolean));
  average=average1/(1-p);stddev=stddev1/(1-p);
  pdf0=1/stddev/sqrt(2*pi)*exp(-(x-average).^2/2/stddev^2);
  pdf=pdf1*p+pdf0*(1-p);
  
 case 'mix_lognormal',
  symExpression=sym(symExpression);
  symBoolean=sym(symBoolean);

  p=getCMoments(net,mdyn,mu,symBoolean);
  [average1,stddev1]=getCMoments(net,mdyn,mu,symExpression*symBoolean);
  average=average1/p;stddev=stddev1/p;
  si=sqrt(log(1+stddev^2/average^2));
  mu=log(average)-si/2;
  pdf1=exp(-((log(x)-mu).^2)/2/si^2)./x/si/sqrt(2*pi);
  
  [average0,stddev0]=getCMoments(net,mdyn,mu,symExpression*(1-symBoolean));
  average=average0/(1-p);stddev=stddev0/(1-p);
  si=sqrt(log(1+stddev^2/average^2));
  mu=log(average)-si/2;
  pdf0=exp(-((log(x)-mu).^2)/2/si^2)./x/si/sqrt(2*pi);
  pdf=pdf1*p+pdf0*(1-p);
  plot(x,(1-p)*pdf0,'-',x,p*pdf1,'--')
  legend('pdf0','pdf1')
  
 otherwise
  error('queryDistribution: unknown distribution type ''%s''\n', ...
	distribution);
end
