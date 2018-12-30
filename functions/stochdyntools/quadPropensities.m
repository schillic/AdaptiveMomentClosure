function [Q,b,c,s,x0]=quadPropensities(net)
% [Q,b,c,s,x0]=quadPropensities(net)
%
% Describes a network of chemical reactions in terms of their
% propensity functions (Q,b,c), stoichiometry (s), and initial
% conditions (x0).
%
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% Outputs
% -------
%
% Q  - matrix with the quadratic terms for the propensity
%      functions. Each reaction corresponds to a (# species) by (#
%      species) square matrix. These matrices are stacked on top of
%      each other so
%         size of Q = (# reactions. # species) by (# species)
%
% b  - matrix with the linear terms for the propensity functions.
%      Each reaction corresponds to a vector of length (#
%      species). These vectors are stacked on top of each other so
%         size of b = (# reactions) by (# species)
%
% c  - vector with the constant terms for the propensity functions.
%      The constants for each reaction are stacked on top of each other so 
%         size of c = (# reactions) by 1
%
% s  - matrix with the stoichiometry for the reactions. Each
%      reaction corresponds to a vector of length (#
%      species). These vector are stacked side by side so
%         size of s = (# species) by (# reactions)
%
% x0 - vector with the initial conditions in the .net  file.
%         size of x0 = (# species) by 1

% Copyright (C) 2008  Joao Hespanha

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
% Created December 19, 2008
%
% January 21, 2009
% Added documentation

mpnb=symengine;

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
x0=zeros(nSpecies,1);
for thisSpecies=1:nSpecies
  symSpecies=[symSpecies,sym(net.species(thisSpecies).id)];
  x0(thisSpecies)=round(double(net.species(thisSpecies).initialAmount));
  
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q=sym(zeros(nSpecies*nReactions,nSpecies));
b=sym(zeros(nReactions,nSpecies));
c=sym(zeros(nReactions,1));
s=sym(zeros(nSpecies,nReactions));
for thisReaction=1:nReactions
    % net.reaction(thisReaction).intensity
    % net.reaction(thisReaction).old
    % net.reaction(thisReaction).new
    % mysubs(symSpecies,net.reaction(thisReaction).old,net.reaction(thisReaction).new,0,mpnb)

    % quadratic term
    Q(1+(thisReaction-1)*nSpecies:thisReaction*nSpecies,:)=jacobian(jacobian(sym(net.reaction(thisReaction).intensity),symSpecies),symSpecies)/2;
    
    % linear term
    b(thisReaction,:)= mysubs(simplify(jacobian(sym(net.reaction(thisReaction).intensity),symSpecies)),symSpecies,zeros(nSpecies,1),0,mpnb);

    % constant term
    c(thisReaction)  = simplify(sym(net.reaction(thisReaction).intensity)-symSpecies*Q(1+(thisReaction-1)*nSpecies:thisReaction*nSpecies,:)*symSpecies.'-b(thisReaction,:)*symSpecies.');

    % stoichiometry
    s(:,thisReaction)=simplify(mysubs(symSpecies,net.reaction(thisReaction).old,net.reaction(thisReaction).new,0,mpnb)-symSpecies)';

end

%delete(mpnb)  % should be used, but apparently creates errors
