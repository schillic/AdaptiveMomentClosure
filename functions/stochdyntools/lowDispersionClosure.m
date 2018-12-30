function approxBarMu=lowDispersionClosure(net,mdyn)
% Auxiliary m-script used by momentClosure.m to compute the moment closure
% functions obtained by setting high-order centered moments to zero.
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
% Created May 12, 2007
%
% June 27, 2007
% Get moment dynamics as a single structure
%
% Dec 31, 2007
% Do computations with unnormalized moments, instead of normalized
% ones, i.e., do not divide by the mean 
%
% Jan 10, 2008
% Fixed to take into account constant species

mpnb=symengine;

verbose=0;

nMoments=size(mdyn.Mu.ndx,1);
nBarMu=size(mdyn.barMu.ndx,1);

nSpecies=length(net.species);
nReactions=length(net.reaction);

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

%% Create symbolic species vector
Phi.mon=symSpecies.';
Phi.sym=expname(net,eye(nSpecies));

% find moments that need to be determined
% . moments that appear in BarMu
% . moments that appear in expansions of BarMu
% but not in mdyn.Mu.ndx

toclose=mdyn.barMu.ndx;
for thisMon=1:size(mdyn.barMu.ndx,1)
  if any(mdyn.barMu.ndx(thisMon,:))
    % compute centered moment
    dmu=prod((Phi.mon-Phi.sym).^(mdyn.barMu.ndx(thisMon,:)'));
    [coef,mon,idx]=monomialindexes(dmu,Phi.mon,mpnb);
    toclose=union(toclose,idx,'rows');
  end    
end

for thisMon=size(toclose,1):-1:1
  % monomial in state?    
  k=find(all(mdyn.Mu.ndx==ones(nMoments,1)*toclose(thisMon,:),2));
  if ~any(toclose(thisMon,:)) || ~isempty(k)
    toclose(thisMon,:)=[];
  end    
end
%toclose

% compute centered moments that must be set to zero
equation='';
unknown='';
for thisMon=1:size(toclose,1)
  % compute centered moment
  dmu=prod((Phi.mon-Phi.sym).^(toclose(thisMon,:)'));
  [coef,mon,idx]=monomialindexes(dmu,Phi.mon,mpnb);
  
  % take expected value
  dmu=0;
  for i=1:size(idx,1)
    if any(idx(i,:))
      dmu=dmu+coef(i)*expname(net,idx(i,:));
    else
      dmu=dmu+coef(i); 
    end    
  end  
  equation=[equation,',',char(dmu)];
  unknown=[unknown,',',char(expname(net,toclose(thisMon,:)))];
end  
if length(equation)>1
  equation(1)=[];
  unknown(1)=[];
end  

% find moments
sol=solve(equation,unknown);

% extract moment to output
approxBarMu=sym(zeros(size(mdyn.barMu.ndx,1),1));

for thisMon=1:size(mdyn.barMu.ndx,1)
  if any(mdyn.barMu.ndx(thisMon,:))
    if verbose  
      fprintf('lowDispersionClosure: closing %s ~ ',char(expname(net,mdyn.barMu.ndx(thisMon,:))))
    end
    
    if isstruct(sol)
      approxBarMu(thisMon)=expand(getfield(sol,char(expname(net,mdyn.barMu.ndx(thisMon,:)))));
    else % single unknown?
      approxBarMu(thisMon)=expand(sol);
    end  
    
    if verbose  
      fprintf('%s\n',char(approxBarMu(thisMon)))
    end    
  else
    approxBarMu(thisMon)=1;
  end    
end

%mdyn.barMu.ndx
%approxBarMu

%delete(mpnb)  % should be used, but apparently create errors
