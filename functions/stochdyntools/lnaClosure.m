function mdyn=lnaClosure(net,mdyn,volume)
% Auxiliary m-script used by momentClosure.m to compute Van Kampen's
% linear noise approximation
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
% Created December 25, 2007 

mpnb=symengine;

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

%% find first order moments 
first=find(sum(mdyn.Mu.ndx,2)==1);

% Deterministic dynamics (Phi)
mdyn.Phi.mon=symSpecies.';
[mdyn.Phi.sym,texrules1]=expname(net,eye(nSpecies),'phi','\phi_{','}');

mdyn1.barMu.ndx=[mdyn.Mu.ndx;mdyn.barMu.ndx];
barPhi=zeroVarianceClosure(net,mdyn1);
mdyn.dotPhi=sym(zeros(nSpecies,1));
mdyn.dotPhi(find(~constSpecies))=[mdyn.dotMu.A(first,:),mdyn.dotMu.B(first,:)]*barPhi;
mdyn.dotPhi=expand(mysubs(mdyn.dotPhi/volume,mdyn.Mu.sym(first),volume*mdyn.Phi.sym,0,mpnb));
mdyn.dotPhi=limit(mdyn.dotPhi,volume,Inf);

% Stochastic Xi
fprintf('\n   ...lnaClosure: computing dotXi(Mu)... working...');

nMoments=size(mdyn.Mu.ndx,1);
mdyn.Xi.mon=sym(zeros(nMoments,1));
[mdyn.Xi.sym,texrules2]=expname(net,mdyn.Mu.ndx,'xi','\xi_{','}');
mdyn.texrules=[texrules2;texrules1;mdyn.texrules];

mdyn.Xi.mu2xi=sym(zeros(nMoments,1));
dotXi=sym(zeros(nMoments,1));
for i=1:nMoments
  mdyn.Xi.mon(i)=prod((mdyn.Phi.mon/volume^(1/2)-mdyn.Phi.sym*volume^(1/2)).^(mdyn.Mu.ndx(i,:)'));

  [coef,mon,idx]=monomialindexes(mdyn.Xi.mon(i),mdyn.Phi.mon,mpnb);
  for thisMon=1:size(idx,1)
    if sum(idx(thisMon,:))==0
      mdyn.Xi.mu2xi(i)=mdyn.Xi.mu2xi(i)+coef(thisMon);
      % derivatives with respect to Phi's  
      for j=1:nSpecies    
	dotXi(i)=dotXi(i)+diff(coef(thisMon),mdyn.Phi.sym(j))*mdyn.dotPhi(j);
      end    
    else
      % monomial in Mu
      k=find(all(mdyn.Mu.ndx==ones(nMoments,1)*idx(thisMon,:),2));
      if isempty(k)
	error('monomial ''%s'' not found in Mu -- this should not happend',char(mon(thisMon)))      
      else
	mdyn.Xi.mu2xi(i)=mdyn.Xi.mu2xi(i)+coef(thisMon)*mdyn.Mu.sym(k);
	% derivatives with respect to Mu's
	dotXi(i)=dotXi(i)+coef(thisMon)*mdyn.dotMu.sym(k);
	% derivatives with respect to Phi's  
	for j=1:nSpecies      
	  dotXi(i)=dotXi(i)+mdyn.Mu.sym(k)*diff(coef(thisMon),mdyn.Phi.sym(j))*mdyn.dotPhi(j);
	end    
      end    
    end      
  end
end

% Stochastic barXi
fprintf('\n   ...lnaClosure: computing barXi...     working...');

nBarMoments=size(mdyn.barMu.ndx,1);
barXi.mon=sym(zeros(nBarMoments,1));
barXi.sym=expname(net,mdyn.barMu.ndx,'xi');
barXi.mu2xi=sym(zeros(nBarMoments,1));
for i=1:nBarMoments
  barXi.mon(i)=prod((mdyn.Phi.mon/volume^(1/2)-mdyn.Phi.sym*volume^(1/2)).^(mdyn.barMu.ndx(i,:)'));
  
  [coef,mon,idx]=monomialindexes(barXi.mon(i),mdyn.Phi.mon,mpnb);
  for thisMon=1:size(idx,1)
    if sum(idx(thisMon,:))==0
      barXi.mu2xi(i)=barXi.mu2xi(i)+coef(thisMon);
    else      
      % monomial in Mu
      k=find(all(mdyn.Mu.ndx==ones(nMoments,1)*idx(thisMon,:),2));
      if isempty(k)
	% monomial in barMu
	k=find(all(mdyn.barMu.ndx==ones(nBarMoments,1)*idx(thisMon,:),2));
	if isempty(k)
	  error('monomial ''%s'' not found in Mu nor barMu -- this should not happend',char(mon(thisMon)))      
	else
	  barXi.mu2xi(i)=barXi.mu2xi(i)+coef(thisMon)*mdyn.barMu.sym(k);
	end	  
      else
	barXi.mu2xi(i)=barXi.mu2xi(i)+coef(thisMon)*mdyn.Mu.sym(k);
      end    
    end      
  end
end

% change variables back to Xi
fprintf('\n   ...lnaClosure: computing dotXi(Xi)... working...');

equations=[mdyn.Xi.mu2xi-mdyn.Xi.sym;barXi.mu2xi-barXi.sym];
k=find(equations~=0);             % only equations not already zero
equations=sym2cell(equations(k));
k=find(sum(mdyn.barMu.ndx,2)>0);  % only "real" variables
variables=sym2cell([mdyn.Mu.sym;mdyn.barMu.sym(k)]);
solution=solve(equations{:},variables{:});

muEq=subssolve(dotXi,variables,solution,mpnb);  % single solution
mdyn.Xi.xi2mu=subssolve(mdyn.Mu.sym,variables,solution,mpnb);
mdyn.Xi.xi2mu=mdyn.Xi.xi2mu{1};
mdyn.Xi.xi2barMu=subssolve(mdyn.barMu.sym,variables,solution,mpnb);
mdyn.Xi.xi2barMu=mdyn.Xi.xi2barMu{1};

mdyn.dotXi=expand(muEq{1});
%mdyn.approxDotXi.sym=mdyn.dotXi;
%mdyn.approxDotXi.sym(5:14)=limit(mdyn.dotXi(5:14),volume,Inf);
mdyn.approxDotXi.sym=limit(mdyn.dotXi,volume,Inf);
% compute approxDotXi.A, approxDotXi.C
mdyn.approxDotXi.A=jacobian(mdyn.approxDotXi.sym,mdyn.Xi.sym);
mdyn.approxDotXi.c=expand(mdyn.approxDotXi.sym-mdyn.approxDotXi.A*mdyn.Xi.sym);

%% find Mu dynamics
fprintf('\n   ...lnaClosure: computing dotMu...     working...');


mdyn.approxDotMu.sym=sym(zeros(nMoments,1));
Mu=subssolve(mdyn.Mu.sym,variables,solution,mpnb);
Mu=Mu{1};  

for i=1:nMoments
  for j=1:nSpecies
    mdyn.approxDotMu.sym(i)=mdyn.approxDotMu.sym(i)+diff(Mu(i),mdyn.Phi.sym(j))*mdyn.dotPhi(j);
  end    

  for j=1:nMoments  
    mdyn.approxDotMu.sym(i)=mdyn.approxDotMu.sym(i)+diff(Mu(i),mdyn.Xi.sym(j))*mdyn.approxDotXi.sym(j);
  end
end

% change variables back to Mu
mdyn.approxDotMu.sym=expand(mysubs(mdyn.approxDotMu.sym,mdyn.Xi.sym,mdyn.Xi.mu2xi,0,mpnb));

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

mdyn.Phi.mon=mdyn.Phi.mon/volume;
mdyn.Phi.x0=mdyn.Mu.x0(first)/double(subsParameters(net,volume,mpnb));

%delete(mpnb)  % should be used, but apparently creates errors

