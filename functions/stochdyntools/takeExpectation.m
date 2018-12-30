function [ELpsi,A,B,Mu,barMu,texrules]=takeExpectation(net,Lpsi,Mu,mpnb)
% [ELpsi,A,B,Mu,barMu,texrules]=takeExpectation(net,Lpsi,Mu)
%
% Take expectation of the symbolic vector Lpsi and express it as a
% linear combination of the moments in Mu annd of an additional
% vector barMu:
%   E[Lpsi] = A E[Mu] + B E[bar Mu]
% This function also returns a vector with E[Mu] and E[bar Mu]
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
% Created March 29, 2010 from code previously in momentDynamics.m
%
% Modified August 25, 2010
% Prevent error when barMu is empty.

verbose=0;

MuPADloops=1;  % do more computations inside MuPAD

if MuPADloops && nargin<4
    mpnb=symengine;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create symbolic species vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSpecies=length(net.species);
nMoments=size(Mu.ndx,1);    

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
    error('takeExpectation: unsupported species type ''%s''',net.species(thisSpecies).type)  
  end  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute matrix A, vector B,vector \bar mu (indices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A=sym(zeros(nMoments));
B=sym(zeros(nMoments,0));

barMu.ndx=zeros(0,nSpecies);

polynomial=1;

for thisMoment=1:nMoments
    if verbose || ~mod(thisMoment,10)
        fprintf('\n   ...takeExpectation: computing A, B   %d/%d... working... ',thisMoment,nMoments);
    end
    
    if 0
        Mu.mon(thisMoment)
        Lpsi(thisMoment)
    end
    
    [coef,mon,idx]=monomialindexes(Lpsi(thisMoment),symSpecies,mpnb);
    try
    catch me
        warning('takeExpectation: moment dynamics are not linear\n\t(not taking expectation of dotMu, no simplifications, not computing A, B, barMu',[]);
        polynomial=0;
        break;
    end
    
    if 0
        Lpsi(thisMoment)
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
            % check if monomial is in Mu
            k=find(all(Mu.ndx==ones(nMoments,1)*idx(thisMon,:),2));
            if isempty(k) % monomial not in Mu
                % check if monomial is in barMu
                k=find(all(barMu.ndx==ones(size(barMu.ndx,1),1)*idx(thisMon,:),2));
                if isempty(k) % monomial not in barMu
                    if isempty(B)
                        B=sym(zeros(nMoments,1)); 
                   end  
                    barMu.ndx(end+1,:)=idx(thisMon,:);
                    B(thisMoment,size(barMu.ndx,1))=coef(thisMon);
                else
                    B(thisMoment,k)=B(thisMoment,k)+coef(thisMon);
                end      
            else
                A(thisMoment,k)=A(thisMoment,k)+coef(thisMon);
            end    
        end    % for thisMon=1:size(idx,1)
    end   % if ~isempty(idx)
end % for thisMoment

if MuPADloops
    %%%% Create barMu.mon
    if isempty(barMu.ndx)
        barMu.mon=sym([]);
    else
        if prod(size(barMu.ndx))>1
            ndx=char(sym(barMu.ndx));
        else
            ndx=sprintf('matrix(1,1,%d)',barMu.ndx);
        end
        cmd=sprintf('zip(matrix(%d,1,1)*%s,%s,_power)',size(barMu.ndx,1),char(symSpecies),ndx);
        % it would be faster to do the prod() inside MuPAD
        barMu.mon=1*prod(evalin(mpnb,cmd),2);  % 1* gets ride of matrix([[]])
    %    barMu.mon
    end
else
    %size(symSpecies(ones(size(barMu.ndx,1),1),:))
    %size(symSpecies(ones(size(barMu.ndx,1),1),:).^barMu.ndx)
    barMu.mon=prod(symSpecies(ones(size(barMu.ndx,1),1),:).^barMu.ndx,2);
end

% 1 should not be in Mu, it should be in barMu
k=find(Mu.mon==1);
if ~isempty(k)
    Mu.ndx(k,:)=[];
    Mu.mon(k,:)=[];
    A(k,:)=[];
    A(:,k)=[];
    B(k,:)=[];
end

if verbose
  fprintf('done\n');
end

if size(barMu.mon,1)==0
  warning('\rtakeExpectation: closed moment dynamics under the assumed species'' types','takeExpectation.closure')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take expected value of monomials in Mu, barMu, dotMu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Mu.sym,rules1]=expname(net,Mu.ndx);
[barMu.sym,rules2]=expname(net,barMu.ndx);
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
texrules=texrules(k,:);

if polynomial
    if isempty(barMu.ndx)
        ELpsi=A*Mu.sym;
    else  
        ELpsi=A*Mu.sym+B*barMu.sym;
    end
else
    ELpsi=Lpsi;
    A=[];
    B=[];
    barMu.ndx=[];
    barMu.mon=[];
    barMu.sym=[];
end

%Mu.mon
%Mu.sym
%barMu.mon
%barMu.sym
%A
%B 
