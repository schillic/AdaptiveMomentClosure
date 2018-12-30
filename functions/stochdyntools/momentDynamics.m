function mdyn=momentDynamics(net,maxdeg)
% mdyn=momentDynamics(net,maxdeg)
%
% This function computes the exact (open) moment dynamics
%           \dot Mu = A Mu + B barMu       ($)
%     where 
%       Mu    - column vectors containing the moments of interest
%       A,B   - appropriate matrices
%       barMu - column vector with higher-order moments
%
% The (exact) time-derivative of the uncentered moment
%   E[psi(X1,X2,...,Xn)],  \psi(X1,X2,\dots,Xn)= X1^m1 X2^m2 ... Xn^mn
% can be obtained by computing
%   d/dt E[psi(X1,X2,...,Xn)] = E[(L psi)((X1,X2,...,Xk)]   (&)
% where $(L psi)(X1,X2,...,Xk) is an expression obtained by applying
% to psi(.) the generator of the Markov process that describes the 
% populations. 
%
% Regardless of the method used, the following is assumed
% in computing the expected value in the left-hand side of (&)
% (1) All replacement rules specified in the 'substitutions' section
%     of the .net file are applied to '(L psi)((X1,X2,...,Xk)' 
%     *before* taking the expected value in the right-hand side of (&).
% (2) For every species X1 that was declared 'deterministic' it is
%     assumed that 
%        E[ X1 f(X2,...Xk) ] = E[ X1 ] E[ f(X2,...Xk) ]
%     for any function $f(.)$ of the remaining species.
% (3) For every species X1 that was declared 'boolean,' it is assumed that
%        X1^n = X1,  for all n >=  1. 
%     This introduces no error as long as the number of molecules of
%     X1 is indeed restricted to the set {0,1}. 
%
% The dynamics in ($) do not provide a closed-system of equations
% because one can only use it to determine Mu if barMu is
% given. However, moment closure techniques can be used to turn ($)
% into a closed system of ODEs [cf. function momentClosure()]
%
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% maxdeg specified which moment to compute. It can be:
%      1) A scalar integer that specifies the largest degree for the
%      uncentered moments in Mu.  
%
%      2) A scalar with fractional part, whose floor() specifies the
%      largest degree for the uncentered moments in Mu, but boolean
%      variables are not taken into account for the degree of a
%      moment. In this case, more moments are included in Mu for the
%      same value of 'maxdeg'.
%
%      3) A matrix with one row per moment and one column per
%      species that species the powers of the monomials to be
%      included in the moment vector 
%      (see .ndx fields in mduyn structure)
%
% Output
% ------
%
% mdyn is a structure characterizing the exact moment dynamics in ($)
%      and the approximate moment dynamics in ($$) or ($$$). 
%      Details can be obtained with 'help mdyn'
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
% Created October 13, 2006 
%
% May 5, 2007
% Added processing of rules
%
% May 7, 2007
% Added initial condition as output
%
% May 29, 2007
% Added dotMu output
%
% June 27, 2007
% Returns everything in a single structure
% Allowed for .5 in maxdeg
%
% June 28, 2007
% Returns texrules
%
% Jan 1st, 2008
% Returns texrules as part of the mdyn structure
%
% January 10, 2007
% Updated inline documentation
%
% Dec 28, 2008
% Replaced subsmaple by evalmaple
% maxdeg can be a .ndx matrix
%
% Modified May 19, 2009
% Replaced evalmaple by mysubs(,,,0)
% Moved from maple to MuPAD
%
% Modified on June 10, 2009
% Added applyRules to loop that computes A, B, and barMu
%
% Modified July 19, 2009
% Removed from Mu monomials that simplify to 0
% Optimized code for MuPAD
%
% Modified March 29, 2010

verbose=0;

MuPADloops=1;  % do more computations inside MuPAD
if MuPADloops
    mpnb=symengine;
end

if isscalar(maxdeg)
    fprintf('computing moment dynamics (maxdeg=%.1f)... ',maxdeg);
else
    fprintf('computing moment dynamics (monomials=[%d,%d])... ',size(maxdeg));
end
t0=clock;

nSpecies=length(net.species);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create symbolic species vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

symSpecies=sym(zeros(1,nSpecies));
boolSpecies=zeros(nSpecies,1);
determSpecies=zeros(nSpecies,1);
stochSpecies=zeros(nSpecies,1);
constSpecies=zeros(nSpecies,1);
for thisSpecies=1:nSpecies
  symSpecies(thisSpecies)=sym(net.species(thisSpecies).id);
    
  if strcmp(net.species(thisSpecies).type,'boolean')
    boolSpecies(thisSpecies)=1;    
  elseif strcmp(net.species(thisSpecies).type,'deterministic')
    determSpecies(thisSpecies)=1;    
  elseif strcmp(net.species(thisSpecies).type,'stochastic')
    stochSpecies(thisSpecies)=1;    
  elseif strcmp(net.species(thisSpecies).type,'constant')
    constSpecies(thisSpecies)=1;    
  else
    error('momentDynamics: unsupported species type ''%s''',net.species(thisSpecies).type)  
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select state Mu.ndx for moment dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mdyn.Mu.ndx=[];

% max 
if isscalar(maxdeg)
    if maxdeg>0
        
        % stochastic 
        thesendx=[0;find(stochSpecies)];
        vv=thesendx;
        for i=1:maxdeg-1
            vv=[kron(vv,ones(length(thesendx),1)),...
                kron(ones(size(vv,1),1),thesendx)];
        end
        % boolean & deterministic (at most once)
        for i=find(boolSpecies | determSpecies)'  
            vv=[kron(vv,ones(2,1)),...
                kron(ones(size(vv,1),1),[0;i])];
        end
        
        mm=zeros(size(vv,1),nSpecies);
        for thisSpecies=1:nSpecies
            if ~strcmp(net.species(thisSpecies).type,'constant')
                mm(:,thisSpecies)=sum(vv==thisSpecies,2);
            end
        end
        mdyn.Mu.ndx=[mdyn.Mu.ndx;mm];
    else
        mdyn.Mu.ndx=[mdyn.Mu.ndx;zeros(1,nSpecies)];
    end
    %% remove monomials with too-high degree
    if maxdeg==floor(maxdeg)
        % all monomials count  
        k=sum(mdyn.Mu.ndx,2);
    else
        % only non boolean monomials count
        k=mdyn.Mu.ndx*(~boolSpecies);
    end
    mdyn.Mu.ndx(find(k>maxdeg),:)=[];
    
    % sort as expected
    mdyn.Mu.ndx=[sum(mdyn.Mu.ndx,2),-mdyn.Mu.ndx];
    mdyn.Mu.ndx=unique(mdyn.Mu.ndx,'rows');
    mdyn.Mu.ndx=-mdyn.Mu.ndx(:,2:end);

elseif size(maxdeg,2)==nSpecies
    mdyn.Mu.ndx=maxdeg;
else
    error('momentDynamics: wrong size for ''maxdeg'' [%d,%d]',size(maxdeg,1),size(maxdeg,2));
end
    
%% remove monomials with zero degree
mdyn.Mu.ndx(find(sum(mdyn.Mu.ndx,2)<1),:)=[];

%mdyn.Mu.ndx
%boolSpecies
%determSpecies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select state Mu.mon for moment dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nMoments=size(mdyn.Mu.ndx,1);
mdyn.Mu.mon=sym(zeros(nMoments,1));
mdyn.Mu.x0=ones(nMoments,1);

if MuPADloops
    %%%% Create Mu.mon
    if prod(size(mdyn.Mu.ndx,1))>1
        ndx=char(sym(mdyn.Mu.ndx));
    else
        ndx=sprintf('matrix(1,1,%d)',mdyn.Mu.ndx);
    end
    cmd=sprintf('zip(matrix(%d,1,1)*%s,%s,_power)',size(mdyn.Mu.ndx,1),char(symSpecies),ndx);
    % it would be faster to do the prod() inside MuPAD
    mdyn.Mu.mon=1*prod(evalin(mpnb,cmd),2);  % 1* gets ride of matrix([[]])
                                             %    mdyn.Mu.mon
else
    %size(symSpecies(ones(size(mdyn.Mu.ndx,1),1),:))
    %size(symSpecies(ones(size(mdyn.Mu.ndx,1),1),:).^mdyn.Mu.ndx)
    mdyn.Mu.mon=prod(symSpecies(ones(size(mdyn.Mu.ndx,1),1),:).^mdyn.Mu.ndx,2);
end

% remove zero monomials
mdyn.Mu.mon=applyRules(net,mdyn.Mu.mon,mpnb);
k=find(mdyn.Mu.mon==0);
mdyn.Mu.mon(k)=[];
mdyn.Mu.ndx(k,:)=[];
mdyn.Mu.x0(k)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create initial condition vector  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for thisMoment=1:nMoments
  for i=1:size(mdyn.Mu.ndx,2)
      if mdyn.Mu.ndx(thisMoment,i)>0
          %          mdyn.Mu.x0(thisMoment)
          %          net.species(i).initialAmount
          %          mdyn.Mu.ndx(thisMoment,i)
          mdyn.Mu.x0(thisMoment)=double(mdyn.Mu.x0(thisMoment)*net.species(i).initialAmount^mdyn.Mu.ndx(thisMoment,i));
      end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply infinitesimal generator to Mu.mon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lpsi=applyGenerator(net,mdyn.Mu.mon,mpnb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute matrix A, vector B,vector \bar mu (indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mdyn.dotMu.sym,mdyn.dotMu.A,mdyn.dotMu.B,mdyn.Mu,mdyn.barMu,mdyn.texrules]=...
    takeExpectation(net,Lpsi,mdyn.Mu,mpnb);
    
fprintf('finished %.2fsec\n',etime(clock,t0));

if MuPADloops
    %delete(mpnb)  % should be used, but apparently create errors
end
