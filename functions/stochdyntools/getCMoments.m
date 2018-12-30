function [average,stddev]=getCMoments(net,mdyn,mu,symExpression,assumezero)
% [average,stddev]=getCMoments(net,mdyn,mu,symExpression,assumezero)
% 
% Given a time series mu of uncentered moments, typically computed
% using ode23s(), this function computes the means and standard
% deviations of a given vector of expressions involving the populations
% of the different species.
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
% mu is a matrix containing a time series of Mu or [Mu;Phi], typically
%    computed using ode23s(). Each row of mu corresponds to a
%    different time instant and each column to the different moment. 
%
% symExpression is an array of symbolic expressions whose means and
%    standard deviations will be computed. The function returns NaN if
%    the computation of a mean and/or standard deviation cannot be
%    computed exactly using the moments available (except when
%    assumezero is nonzero, see below).
%
%    The replacement rules in net.substitutionrule will be applied to
%    the symbolic expression (and to its square in order to compute
%    the variance) before taking expectations. The boolean nature 
%    of the variables declared as such, will be taken into account.
% 
% assumezero is an optional boolean variable that, when nonzero, 
%    instructs the function to assume zero all moments that are not
%    available to compute the mean and/or standard deviations.
%
% Attention: In the current implementation of this function, the
%    deterministic nature of the variables declared as such will be
%    ignore. This means, e.g., that NaN will be generated if one asks
%    for the expected value of X^2 and mu only contains the expected
%    value of X, even if X was declared deterministic.
%
% Outputs
% -------
%
% average is a time series of the expected value of the vector
%         symExpression, with one time instant per row and one element
%         of the symExpression per column. 
%
% stddev (optional) is a time series of the standard deviation of the
%        vector symExpression, with one time instant per row and one
%        element of the symExpression per column.
%
% Attention: For populations with low stochasticity (i.e., for which
%    standard deviations are much smaller than the means), the moment
%    closure approximation errors can lead to negative variances. Such
%    populations should be declared as deterministic in the .net file.
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
% June 27, 2007
% mdyn passed as a single structure
%
% December 26, 2007
% made compatible with 'lna' moment closure
%
% January 10, 2007
% updated inline documentation
%
% July 18, 2009
% added assumezero input

if nargin<5
    assumezero=0;
end

mpnb=symengine;

if isfield(mdyn,'Phi')  
  % LNA  
  if length(mdyn.Mu.mon)+length(mdyn.Phi.mon)~=size(mu,2)
    error('getCMoments: mdyn.Mu.mon has %d elements \& mdyn.Phi.mon has %d elements, whose sum is different than the number of columns %d of mu\n',length(mdyn.Mu.mon),length(mdyn.Phi.mon),size(mu,2));
  end
  mu=mu(:,1:length(mdyn.Mu.mon));
else
  if length(mdyn.Mu.mon)~=size(mu,2)
    error('getCMoments: mdyn.Mu.mon has %d elements, which is different than the number of columns %d of mu\n', length(mdyn.Mu.mon),size(mu,2));
  end
end

nSpecies=length(net.species);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create symbolic species vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symSpecies=[];
boolSpecies=zeros(nSpecies,1);
determSpecies=zeros(nSpecies,1);
stochSpecies=zeros(nSpecies,1);
constSpecies=zeros(nSpecies,1);
for thisSpecies=1:nSpecies
  symSpecies=[symSpecies,sym(net.species(thisSpecies).id)];
    
  if strcmp(net.species(thisSpecies).type,'boolean')
    boolSpecies(thisSpecies)=1;    
  elseif strcmp(net.species(thisSpecies).type,'deterministic')
    determSpecies(thisSpecies)=1;    
  elseif strcmp(net.species(thisSpecies).type,'stochastic')
    stochSpecies(thisSpecies)=1;    
  elseif strcmp(net.species(thisSpecies).type,'constant')
    constSpecies(thisSpecies)=1;    
  else
    error('unsupported species type ''%s''',net.species(thisSpecies).type)  
  end  
end

symExpression=sym(symExpression);

average=zeros(size(mu,1),length(symExpression));

for i=1:length(symExpression)
%  applyRules(net,symExpression(i))
  [coef,mon,idx]=monomialindexes(applyRules(net,symExpression(i)),symSpecies,mpnb);
  
  %% Boolean simplification: x^n -> x  
  k=find(idx & (ones(size(idx,1),1)*boolSpecies'));
  idx(k)=1;
  
  %% Constant simplification: E[x^n] -> initialAmount^n  
  k=find(idx .* (ones(size(idx,1),1)*constSpecies'));
  [id,sp]=ind2sub(size(idx),k);
%coef  
  for j=1:length(id)
    e=zeros(1,nSpecies);
    e(sp(j))=1;    
    % add to coeficient
    coef(id(j))=coef(id(j))*net.species(sp(j)).initialAmount^idx(id(j),sp(j));
    % remove from monomials    
    idx(id(j),sp(j))=0;
  end
%coef  
%  idx
  
  for k=1:length(mon)  % loop over monomials
    if coef(k)~=0
      if any(idx(k,:))
	mon(k)=prod(symSpecies.^idx(k,:));
	found=0;
	for j=1:length(mdyn.Mu.mon) % find monomial in  mdyn.Mu.mon
	  if mdyn.Mu.mon(j)-mon(k)==0
	    try
	      average(:,i)=average(:,i)+double(coef(k))*mu(:,j);
	    catch
	      error('getCMoments: coefficient ''%s'' in expression ''%s'' is not a double\n',char(coef(k)),char(symExpression(i)));
	    end
	    found=1;
	    break;
	  end
	end
	if ~found
          if ~assumezero
              average(:,i)=average(:,i)+NaN;
              warning('getCMoments: monomial ''%s'' needed to compute the mean of ''%s''\n',...
                      char(mon(k)),char(symExpression(i)));
          else
              warning('getCMoments: monomial ''%s'' needed to compute the mean of ''%s'' - assumedzero\n',...
                      char(mon(k)),char(symExpression(i)));

          end
	end
      else % monomial is just 1
	try
	  average(:,i)=average(:,i)+double(coef(k));
	catch
	  error('getCMoments: coefficient ''%s'' in expression ''%s'' is not a double\n',char(coef(k)),char(symExpression(i)));
	end
      end
    end
  end
end

if nargout<2
  return
end

stddev=zeros(size(mu,1),length(symExpression));

for i=1:length(symExpression)
%  symExpression(i)^2
%  applyRules(net,symExpression(i)^2)
  [coef,mon,idx]=monomialindexes(applyRules(net,symExpression(i)^2),symSpecies,mpnb);

  %% Boolean simplification: x^n -> x  
  k=find(idx & (ones(size(idx,1),1)*boolSpecies'));
  idx(k)=1;  
  
  %% Constant simplification: E[x^n] -> initialAmount^n  
  k=find(idx .* (ones(size(idx,1),1)*constSpecies'));
  [id,sp]=ind2sub(size(idx),k);
%coef  
  for j=1:length(id)
    e=zeros(1,nSpecies);
    e(sp(j))=1;    
    % add to coeficient
    coef(id(j))=coef(id(j))*net.species(sp(j)).initialAmount^idx(id(j),sp(j));
    % remove from monomials    
    idx(id(j),sp(j))=0;
  end
%coef  
%  idx

for k=1:length(mon)  % loop over monomials
    if coef(k)~=0
      if any(idx(k,:))
	mon(k)=prod(symSpecies.^idx(k,:));
	found=0;
	for j=1:length(mdyn.Mu.mon) % find monomial in  mdyn.Mu.mon
	  if mdyn.Mu.mon(j)-mon(k)==0
	    try
	      stddev(:,i)=stddev(:,i)+double(coef(k))*mu(:,j);
	    catch
	      error('getCMoments: coefficient ''%s'' in expression ''%s'' is not a double\n',char(coef(k)),char(symExpression(i)));
	    end
	    found=1;
	    break;
	  end
	end
	if ~found
          if ~assumezero
             stddev(:,i)=stddev(:,i)+NaN;
             warning('getCMoments: monomial ''%s'' needed to compute the standard deviation of ''%s''\n',...
                     char(mon(k)),char(symExpression(i)));
          else
              warning('getCMoments: monomial ''%s'' needed to compute the standard deviation of ''%s'' - assumed zero\n',...
		  char(mon(k)),char(symExpression(i)));
          end
	end
      else % monomial is just 1
	try
	  stddev(:,i)=stddev(:,i)+double(coef(k));
	catch
	  error('getCMoments: coefficient ''%s'' in expression ''%s'' is not a double\n',char(coef(k)),char(symExpression(i)));
	end
      end
    end
  end
end

stddev=sqrt(stddev-average.^2);

%delete(mpnb)  % should be used, but apparently creates errors
