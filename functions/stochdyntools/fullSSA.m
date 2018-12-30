function [Tr,X] = fullSSA(Q,b,c,s,X0,N,T)
% [Tr,X] = ssa(Q,b,c,s,X0,N,T)
%
% Runs Gillespie's stochastic Simulation Algorithms (SSA)
%
% Inputs:
%
% Q  - Matrix with the quadratic terms for the propensity function.
%      The (# species) by (# species) matrices for each reaction are
%      stacked on top of each other so the size of Q is (# reactions
%      . # species) by # species
%
% b  - Matrix with the linear terms for the propensity function
%      The # species length vectors for each reaction are stacked on
%      top of each other so the size of b is (# reactions) by (# species)
%
% c  - Vector with the constant terms for the propensity function
%      The constants for each reaction are stacked on top of each other
%      so the size of c is (# reactions) by 1
%
% s  - Matrix with the stoichiometry for the reactions
%      The stoichiometry vectore for each reaction are stacked side
%      by side so the size of s is (# species) x (# reactions)
%
% X0 - Vector with initial conditions. When X0 is a column vector with
%      length equal to # species, all Monte Carlo simulations start
%      from the same initial condition. Otherwise, X0 should be a
%      matrix with size (# species) by (# simulations).
%
% N  - # of Monte Carlo simulations to run
%
% T  - # of reaction steps for each simulation
%
% Outputs:
%
% Tr - Matrix with the reaction times. The size of Tr is 
%      (# of reaction steps) x (# simulations)
%
% X  - Matrix with the molecule counts at the reaction times. The
%      size of X is (# species) x (# simulations) x (# of reaction steps)
%
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
% Created December 22, 2008


t0=clock;
fprintf('running fullSSA... ');

[k,n]=size(b);           % n = # species, k = # chemical reactions

Tr=zeros(T,N,'single');  % times
X=zeros(n,N,T,'single'); % populations

% all parameters must be real
Q=sparse(double(Q));
b=sparse(double(b));
c=sparse(double(c));
s=double(s);


if size(X0)==[n,1]
    Xt=X0(:,ones(1,N));
elseif size(X0)==[n,N]
    Xt=X0;
else
    error('fullSSA: wrong size for initial condition [%d,%d]\n',size(X0))
end

Expand=repmat((1:n)',[k,1]);
Sum=eye(k,'double');
Sum=sparse(Sum(:,sort(repmat((1:k)',[n,1]))));
c=c(:,ones(1,N));

ti=1;
X(:,:,ti)=Xt;
t=zeros(1,N);
while (ti<T)
    % compute propensity functions
    lambda=Sum*(Xt(Expand,:).*(Q*Xt))+b*Xt+c;
    
    % compute intensity for first reaction
    sumlambda=sum(lambda,1);
    % compute next time steps
    t=t-log(rand(1,N))./sumlambda;

    % compute probability of each reaction
    lambda=lambda./(ones(k,1)*sumlambda);
    % select reaction
    j=rand_discr(lambda);
    Xt=Xt+s(:,j);
    
    % save data
    ti=ti+1;
    Tr(ti,:)=t;
    X(:,:,ti)=Xt;
end

fprintf('finished %.2fsec\n',etime(clock,t0));

