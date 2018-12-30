function centeredMu=centerMoments(mdyn,mu)
% centeredMu=centerMoments(mdyn,mu)
% 
% This function takes a vector of uncentered moments and produces a vector
% of centered moments.
% 
% Inputs
% ------
%
% mdyn : structure describing the moment dynamics, as returned by 
%        closureDynamics() or momentDynamics().
%
% mu : Time series of uncentered moments with one time instant per row
%      and one moment per columns. Such vectors are typically obtained
%      from an ode solver that computes the solution to some
%      approximate moment dynamics. The meaning of each column is
%      specified by mdyn.Mu
%
% Outputs
% -------
%
% centeredMu : Time series of centered moments.
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
% Created October 27, 2006 
%
% Modified June 27, 2007
% mdyn passed as a single structure

mpnb=symengine;

nSpecies=size(mdyn.Mu.ndx,2);

symSpecies=[];
symMeans=[];
meansMu=zeros(size(mu,1),nSpecies);
for thisSpecies=1:nSpecies
  symSpecies=[symSpecies;sym(sprintf('S%d',thisSpecies))];
  symMeans=[symMeans;sym(sprintf('M%d',thisSpecies))];
  
  % find mean
  idx=zeros(1,nSpecies);
  idx(thisSpecies)=1;
  kk=find(all(mdyn.Mu.ndx==ones(size(mdyn.Mu.ndx,1),1)*idx,2));
  if isempty(kk)
    error('centerMoments: failed to find mean of ''%s'' in monomial vector',char(symSpecies(thisSpecies)));
  end
  meansMu(:,thisSpecies)=mu(:,kk);
end
%symSpecies
%symMeans

centeredMu=zeros(size(mu));

for thisMoment=1:size(mu,2)
  if sum(mdyn.Mu.ndx(thisMoment,:))==1
    % do not center first-order moments
    centeredMu(:,thisMoment)=mu(:,thisMoment);
  else
    centered=prod((symSpecies-symMeans).^(mdyn.Mu.ndx(thisMoment,:)'));
    [coef,mon,idx]=monomialindexes(centered,symSpecies,mpnb);
    
    for i=1:length(coef)
      [ccoef,mmon,iidx]=monomialindexes(coef(i),symMeans,mpnb);
      m=zeros(size(mu,1),1);
      for j=1:length(ccoef)
	m=m+double(ccoef(j))*prod(meansMu.^(ones(size(meansMu,1),1)*iidx(j,:)),2);
      end
      
      % find the monomial 
      if any(idx(i,:))
	k=find(all(mdyn.Mu.ndx==ones(size(mdyn.Mu.ndx,1),1)*idx(i,:),2));
	if isempty(k)
	  error('centerMoments: failed to find monomial ''%s'' in monomial vector',char(mon(i)));
	end
	m=m.*mu(:,k);
      end
      
      centeredMu(:,thisMoment)=centeredMu(:,thisMoment)+m;
    end
  end
end
  
%delete(mpnb)  % should be used, but apparently creates errors
