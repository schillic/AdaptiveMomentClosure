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
nReactions=length(net.reaction);
nDerivatives=length(net.raterule);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% select state Mu for moment dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute moment dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%

polynomial=1;

nMoments=size(mdyn.Mu.ndx,1);
mdyn.Mu.mon=sym(zeros(nMoments,1));
mdyn.Mu.x0=ones(nMoments,1);

if MuPADloops
    %%%% computation in MuPAD
    intensity=sprintf('%s,',net.reaction(:).intensity);
    intensity=['intensity:=[',intensity(1:end-1),']'];
    
    clear old new
    [old{1:nReactions,1}]=deal(net.reaction(:).old);
    old=cellfun(@(x) (['[',sprintf('%s,',x{:}),']']),old,'UniformOutput',0);
    old=regexprep(['olds:=[',sprintf('%s,',old{:}),']'],',]',']');
    [new{1:nReactions,1}]=deal(net.reaction(:).new);
    new=cellfun(@(x) (['[',sprintf('%s,',x{:}),']']),new,'UniformOutput',0);
    new=regexprep(['news:=[',sprintf('%s,',new{:}),']'],',]',']');
    evalin(mpnb,intensity);
    evalin(mpnb,old);
    evalin(mpnb,new);
end
  

if MuPADloops
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
mdyn.Mu.mon=applyRules(net,mdyn.Mu.mon);
k=find(mdyn.Mu.mon==0);
mdyn.Mu.mon(k)=[];
mdyn.Mu.ndx(k,:)=[];
mdyn.Mu.x0(k)=[];

nMoments=size(mdyn.Mu.ndx,1);

derivative=sym(zeros(nMoments,1));

for thisMoment=1:nMoments
  if verbose || ~mod(thisMoment,10)
    fprintf('\n   ...momentDynamics: computing moment %d/%d... working... ',thisMoment,nMoments);
  end

  % create initial condition vector  
  for i=1:size(mdyn.Mu.ndx,2)
      if mdyn.Mu.ndx(thisMoment,i)>0
          %          mdyn.Mu.x0(thisMoment)
          %          net.species(i).initialAmount
          %          mdyn.Mu.ndx(thisMoment,i)
          mdyn.Mu.x0(thisMoment)=double(mdyn.Mu.x0(thisMoment)*net.species(i).initialAmount^mdyn.Mu.ndx(thisMoment,i));
      end
  end
  
  %%% generator formula - loop over reactions
  if MuPADloops
      %%%% computation in MuPAD
      expression=char(mdyn.Mu.mon(thisMoment));
      cmd=sprintf('der:=0;_for(k,1,nops(intensity),1,(der:=expand(der+intensity[k]*(subs(%s,zip(olds[k],news[k],_equal))-%s))));',expression,expression);
      derivative(thisMoment)=evalin(mpnb,cmd);
  else
      %%%% computations in MATLAB
      derivative(thisMoment)=0;
      for thisReaction=1:nReactions
          derivative(thisMoment)=expand(derivative(thisMoment)+sym(net.reaction(thisReaction).intensity)*(mysubs(mdyn.Mu.mon(thisMoment),net.reaction(thisReaction).old,net.reaction(thisReaction).new,0)-mdyn.Mu.mon(thisMoment)));
          if 0
              mdyn.Mu.mon(thisMoment)
              net.reaction(thisReaction).intensity
              [net.reaction(thisReaction).old,net.reaction(thisReaction).new]
              
              [mdyn.Mu.mon(thisMoment),mysubs(mdyn.Mu.mon(thisMoment),net.reaction(thisReaction).old,net.reaction(thisReaction).new,0)]
              
              expand(derivative(thisMoment))
          end
      end  % for thisReaction
  end
  
  %%% generator formula - loop over derivatives
  for thisDerivative=1:nDerivatives
    derivative(thisMoment)=expand(derivative(thisMoment)+diff(mdyn.Mu.mon(thisMoment),net.raterule(thisDerivative).variable)*net.raterule(thisDerivative).formula);

    if 0
        mdyn.Mu.mon(thisMoment)
        net.raterule(thisDerivative).variable
        net.raterule(thisDerivative).formula
        
        diff(mdyn.Mu.mon(thisMoment),net.raterule(thisDerivative).variable)
        
        derivative(thisMoment)
    end
  end  % for thisDerivative

  % initial assignment in case next loop fails
  mdyn.dotMu.sym(thisMoment,1)=derivative(thisMoment);
end % for thisMoment


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute matrix A, vector B,vector \bar mu (indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% apply substitution rules
derivative=applyRules(net,derivative);
  
mdyn.dotMu.A=sym(zeros(nMoments));
mdyn.dotMu.B=sym(zeros(nMoments,0));

mdyn.barMu.ndx=zeros(0,nSpecies);

for thisMoment=1:nMoments
    if verbose || ~mod(thisMoment,10)
        fprintf('\n   ...momentDynamics: computing A,B %d/%d... working... ',thisMoment,nMoments);
    end
    
    if 0
        mdyn.Mu.mon(thisMoment)
        derivative(thisMoment)
        expand(derivative(thisMoment))
    end
    
    [coef,mon,idx]=monomialindexes(derivative(thisMoment),symSpecies);
    try
    catch me
        warning('momentDynamics: moment dynamics are not linear\n\t(not taking expectation of dotMu, no simplifications, not computing A, B, barMu',[]);
        mdyn.dotMu.sym=expand(mdyn.dotMu.sym);
        mdyn.dotMu.A=[];
        mdyn.dotMu.B=[];
        polynomial=0;
        break;
    end
    
    if 0
        derivative(thisMoment)
        coef
        mon
        idx  
    end
    
    if ~isempty(idx)
        %% Boolean simplification: x^n -> x  
        k=find(idx & (ones(size(idx,1),1)*boolSpecies'));
        idx(k)=1;  
        
        %% Deterministic simplification: x^n -> E[x]^n  
        k=find(idx .* (ones(size(idx,1),1)*determSpecies'));
        [id,sp]=ind2sub(size(idx),k);
        %coef  
        for i=1:length(id)
            e=zeros(1,nSpecies);
            e(sp(i))=1;    
            % add to coeficient
            coef(id(i))=coef(id(i))*expname(net,e)^idx(id(i),sp(i));
            % remove from monomials    
            idx(id(i),sp(i))=0;  
        end
        %coef  
        %  idx
        
        %% Constant simplification: E[x^n] -> x^n  
        k=find(idx .* (ones(size(idx,1),1)*constSpecies'));
        [id,sp]=ind2sub(size(idx),k);
        %coef  
        for i=1:length(id)
            e=zeros(1,nSpecies);
            e(sp(i))=1;    
            % add to coeficient
            coef(id(i))=coef(id(i))*symSpecies(sp(i))^idx(id(i),sp(i));
            % remove from monomials    
            idx(id(i),sp(i))=0;
        end
        %coef  
        %  idx
        
        %% construct matrices A and B and vector muBar
        for thisMon=1:size(idx,1)
            % monomial in Mu  
            k=find(all(mdyn.Mu.ndx==ones(nMoments,1)*idx(thisMon,:),2));
            if isempty(k)
                % monomial already in barMu
                k=find(all(mdyn.barMu.ndx==ones(size(mdyn.barMu.ndx,1),1)*idx(thisMon,:),2));
                if isempty(k)
                    if isempty(mdyn.dotMu.B)
                        mdyn.dotMu.B=sym(zeros(nMoments,1));	  
                    end  
                    mdyn.barMu.ndx(end+1,:)=idx(thisMon,:);
                    mdyn.dotMu.B(thisMoment,size(mdyn.barMu.ndx,1))=coef(thisMon);
                else
                    mdyn.dotMu.B(thisMoment,k)=mdyn.dotMu.B(thisMoment,k)+coef(thisMon);
                end      
            else
                mdyn.dotMu.A(thisMoment,k)=mdyn.dotMu.A(thisMoment,k)+coef(thisMon);
            end    
        end    % for thisMon=1:size(idx,1)
    end   % if ~isempty(idx)
end % for thisMoment

if MuPADloops
    if prod(size(mdyn.barMu.ndx))>1
        ndx=char(sym(mdyn.barMu.ndx));
    else
        ndx=sprintf('matrix(1,1,%d)',mdyn.barMu.ndx);
    end
    cmd=sprintf('zip(matrix(%d,1,1)*%s,%s,_power)',size(mdyn.barMu.ndx,1),char(symSpecies),ndx);
    % it would be faster to do the prod() inside MuPAD
    mdyn.barMu.mon=1*prod(evalin(mpnb,cmd),2);  % 1* gets ride of matrix([[]])
    %    mdyn.barMu.mon
else
    %size(symSpecies(ones(size(mdyn.barMu.ndx,1),1),:))
    %size(symSpecies(ones(size(mdyn.barMu.ndx,1),1),:).^mdyn.barMu.ndx)
    mdyn.barMu.mon=prod(symSpecies(ones(size(mdyn.barMu.ndx,1),1),:).^mdyn.barMu.ndx,2);
end

% 1 should not be in Mu, it should be in barMu
k=find(mdyn.Mu.mon==1);
if ~isempty(k)
    mdyn.Mu.ndx(k,:)=[];
    mdyn.Mu.mon(k,:)=[];
    mdyn.dotMu.A(k,:)=[];
    mdyn.dotMu.A(:,k)=[];
    mdyn.dotMu.B(k,:)=[];
end

if verbose
  fprintf('done\n');
end

if size(mdyn.barMu.mon,1)==0
  warning('\rmomentDynamics: closed moment dynamics under the assumed species'' types','momentDynamics.closure')
end

[mdyn.Mu.sym,rules1]=expname(net,mdyn.Mu.ndx);
[mdyn.barMu.sym,rules2]=expname(net,mdyn.barMu.ndx);
texrules=[rules2;rules1];
for i=1:length(net.parameter)
  id=net.parameter(i).id;
  [starts,ends]=regexp(id,'_','start','end');
  for j=length(ends):-1:1
    id=[id(1:starts(j)-1),'\\',id(ends(j):length(id))];
  end     
  texrules=[texrules;{id,net.parameter(i).idTex}];
end
for i=1:length(net.species)
  id=net.species(i).id;
  [starts,ends]=regexp(id,'_','start','end');
  for j=length(ends):-1:1
    id=[id(1:starts(j)-1),'\\',id(ends(j):length(id))];
  end     
  texrules=[texrules;{id,net.species(i).idTex}];
end
% reverse order to avoid replacement of partial names
[dummy,k]=sort(sum(char(texrules(:,1))~=' ',2),1,'descend');
mdyn.texrules=texrules(k,:);

if polynomial
    if isempty(mdyn.barMu.ndx)
        mdyn.dotMu.sym=mdyn.dotMu.A*mdyn.Mu.sym;
    else  
        mdyn.dotMu.sym=mdyn.dotMu.A*mdyn.Mu.sym+mdyn.dotMu.B*mdyn.barMu.sym;
    end
end

%mdyn.Mu.mon
%mdyn.Mu.sym
%mdyn.barMu.mon
%mdyn.barMu.sym
%mdyn.dotMu.A
%mdyn.dotMu.B 

fprintf('finished %.2fsec\n',etime(clock,t0));
