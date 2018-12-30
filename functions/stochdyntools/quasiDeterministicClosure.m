function mdyn=quasiDeterministicClosure(net,mdyn)
% Auxiliary m-script used by momentClosure.m to compute the moment closure
% functions obtained by setting high-order quasi-deterministic moments
% to zero.
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
% Created Dec 31, 2007
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

%% find first order moments 
first=find(sum(mdyn.Mu.ndx,2)==1);

% Deterministic dynamics (Phi)
mdyn.Phi.mon=symSpecies.';
[mdyn.Phi.sym,texrules]=expname(net,eye(nSpecies),'phi','\phi_{','}');
mdyn.texrules=[texrules;mdyn.texrules];

mdyn.Phi.sym(find(constSpecies))=mdyn.Phi.mon(find(constSpecies));

mdyn1.barMu.ndx=[mdyn.Mu.ndx;mdyn.barMu.ndx];
barPhi=zeroVarianceClosure(net,mdyn1);
mdyn.dotPhi=sym(zeros(nSpecies,1));
mdyn.dotPhi(find(~constSpecies))=[mdyn.dotMu.A(first,:),mdyn.dotMu.B(first,:)]*barPhi;
mdyn.dotPhi=expand(mysubs(mdyn.dotPhi,mdyn.Mu.sym(first),mdyn.Phi.sym,0,mpnb));

% find moments that need to be determined
% . moments that appear in BarMu
% . moments that appear in expansions of BarMu
% but not in mdyn.Mu.ndx

toclose=mdyn.barMu.ndx;
for thisMon=1:size(mdyn.barMu.ndx,1)
  if any(mdyn.barMu.ndx(thisMon,:))
    % compute quasi-determinsitic centered moment
    dmu=prod((mdyn.Phi.mon-mdyn.Phi.sym).^(mdyn.barMu.ndx(thisMon,:)'));
    [coef,mon,idx]=monomialindexes(dmu,mdyn.Phi.mon,mpnb);
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
  dmu=prod((mdyn.Phi.mon-mdyn.Phi.sym).^(toclose(thisMon,:)'));
  
  [coef,mon,idx]=monomialindexes(dmu,mdyn.Phi.mon,mpnb);
  
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
mdyn.approxDotMu.barMu=sym(zeros(size(mdyn.barMu.ndx,1),1));

for thisMon=1:size(mdyn.barMu.ndx,1)
  if any(mdyn.barMu.ndx(thisMon,:))
    if verbose  
      fprintf('lowDispersionClosure: closing %s ~ ',char(expname(net,mdyn.barMu.ndx(thisMon,:))))
    end
    
    if isstruct(sol)
      mdyn.approxDotMu.barMu(thisMon)=expand(getfield(sol,char(expname(net,mdyn.barMu.ndx(thisMon,:)))));
    else % single unknown?
      mdyn.approxDotMu.barMu(thisMon)=expand(sol);
    end  
    
    if verbose  
      fprintf('%s\n',char(approxDotMu.barMu(thisMon)))
    end    
  else
    mdyn.approxDotMu.barMu(thisMon)=1;
  end    
end

%mdyn.barMu.ndx
%approxDotMu.barMu

% compute approxDotMu
if isempty(mdyn.dotMu.B)
  mdyn.approxDotMu.sym=expand(mdyn.dotMu.A*mdyn.Mu.sym);
else
  mdyn.approxDotMu.sym=expand(mdyn.dotMu.A*mdyn.Mu.sym...
      +mdyn.dotMu.B*mdyn.approxDotMu.barMu);
end

% compute approxDotMu.A, approxDotMu.C
mdyn.approxDotMu.A=jacobian(mdyn.approxDotMu.sym,mdyn.Mu.sym);
mdyn.approxDotMu.c=expand(mdyn.approxDotMu.sym-mdyn.approxDotMu.A*mdyn.Mu.sym);

% remove constants from Phi vector

mdyn.Phi.sym(find(constSpecies))=[];
mdyn.Phi.mon(find(constSpecies))=[];
mdyn.dotPhi(find(constSpecies))=[];

if length(mdyn.Phi.mon)~=length(mdyn.Mu.mon(first)) ...
      || any(mdyn.Phi.mon~=mdyn.Mu.mon(first))
  error('Mu vector without 1st order moments in the right order')
end

mdyn.Phi.x0=mdyn.Mu.x0(first);

%delete(mpnb)  % should be used, but apparently create errors
