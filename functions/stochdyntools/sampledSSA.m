function [X,Xmean,Xstd,XmeanCI,XstdCI] = sampledSSA(Q,b,c,s,x0,nMC,Ts,CI)
% [X,Xmean,Xstd,XmeanCI,XstdCI] = sampledSSA(Q,b,c,s,x0,nMC,Ts,CI)
%
% Runs Gillespie's stochastic Simulation Algorithms (SSA) for the
% network of chemical reactions with propsenties specified by Q,b,c
% and stoichiometry specified by s.
%
% This function returns all the sample paths at a set of sample
% times, as well as estimates for the means and standard deviations
% at the sample times.
%
% Inputs:
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
% x0 - vector with initial conditions. If x0 is a column vector with
%         size of x0 = (# species) by 1 
%      then all Monte Carlo simulations start from the same initial
%      condition. Otherwise, 
%         size of x0 = (# species) by (# simulations)
%      and each Monte Carlo simulations starts with the initial
%      condition given by the corresponding columns of x0
%
% nMC - # of Monte Carlo simulations to run
%
% Ts - Desired sample times for output.  The sampled paths are
%      computed exactly at all times, but their values are only
%      returned at the time in the vector Ts.
%
% CI  - (optional) desired percentage for the confidence
%       interval for the mean and standard deviation. If not
%       specified, 95% is used.
%
% Outputs:
%
% X  - Matrix with the molecule counts at the sample times. The
%      size of X is (# species) x (# simulations) x (# of sample times)
%
% Xmean - Matrix with the mean molecule counts at the sample times.
%         The size of Xmean is (# species) x (# of sample times)
%
% Xstd  - Matrix with the std. dev. molecule counts at the sample times.
%         The size of Xstd is (# species) x (# of sample times)
%
% XmeanCI - Matrix with the confidence interval for the Xmean.
%           The size of XmeanCI is (# species) x (# of sample times) x 2
%
%           The computation of XmeanCI assumes that the central
%           limit theorem is valid to determine the distribution of
%           the mean
%
% XstdCI - Matrix with the confidence interval for the Xstd.
%            The size of XstdCI is (# species) x (# of sample times) x 2
%
%            The computation of XstdCI assumes a normal distribution.
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

if nargin<8
    CI =.95;
end

t0=clock;
fprintf('running sampledSSA... ');

[nReactions,nSpecies] = size(b);
nTs   = length(Ts);

X=NaN(nSpecies,nMC,nTs,'single'); % populations

% all parameters must be real
Q=sparse(double(Q));
b=sparse(double(b));
c=sparse(double(c));
s=double(s);

if size(x0)==[nSpecies,1]
    Xt=x0(:,ones(1,nMC));
elseif size(x0)==[nSpecies,nMC]
    Xt=x0;
else
    error('sampledSSA: wrong size for initial condition [%d,%d]\n',size(x0))
end

Expand=repmat((1:nSpecies)',[nReactions,1]);
Sum=eye(nReactions,'double');
Sum=sparse(Sum(:,sort(repmat((1:nReactions)',[nSpecies,1]))));
c=c(:,ones(1,nMC));

t=zeros(1,nMC);
alive=1:nMC;
next=ones(1,nMC);
Ts=Ts(:)';
dt=zeros(1,nMC);
reaction=ones(1,nMC);
lambda=zeros(nReactions,nMC);
count=0;
while (~isempty(alive))
    count=count+1;
    if mod(count,ceil(1000000/nMC))==0
        fprintf('t=%f ',min(t));
    end
    % compute propensity functions
    lambda=Sum*(Xt(Expand,:).*(Q*Xt))+b*Xt+c;
    
    % compute intensity for first reaction
    sumlambda=sum(lambda,1);
    % compute next time steps
    t(alive)=t(alive)-log(rand(size(alive)))./sumlambda(alive);

    % save sampled data
    ksim=find((next<=nTs) & (t>Ts(min(next,nTs))));   % any simulation reached next sampling time?
    while (~isempty(ksim))
        % option 1 - loop over reactions (sub2ind seems very slow)
        if (length(ksim)>1)
            ndx=sub2ind([nMC,nTs],ksim,next(ksim));
            for j=1:nSpecies
                X(j,ndx) = Xt(j,ksim);
            end
        else
            X(:,ksim,next(ksim)) = Xt(:,ksim);
        end

        next(ksim)=next(ksim)+1; % update index of next sampling time?

        % remove simulations that terminated
        alive=find(next<=nTs);

        ksim=find((next<=nTs) & (t>Ts(min(next,nTs))));   % any simulation reached next sampling time?
    end

    % select reaction
    if ~isempty(alive)
        reaction(alive)=rand_discr(lambda(:,alive)./(ones(nReactions,1)*sumlambda(alive)));
        Xt=Xt+s(:,reaction);
    end
end

if nargout>1
    % confidence interval, along 2nd dimension
    [Xmean,Xstd,XmeanCI1,XmeanCI2,XstdCI1,XstdCI2] = meanStdCI(X,CI,2);
    Xmean = reshape(Xmean,nSpecies,nTs);
    Xstd  = reshape(Xstd,nSpecies,nTs);

    XmeanCI = zeros(nSpecies,nTs,2,'single');
    XstdCI  = zeros(nSpecies,nTs,2,'single');

    XmeanCI(:,:,1)=reshape(XmeanCI1,nSpecies,nTs);
    XmeanCI(:,:,2)=reshape(XmeanCI2,nSpecies,nTs);
    XstdCI(:,:,1)=reshape(XstdCI1,nSpecies,nTs);
    XstdCI(:,:,2)=reshape(XstdCI2,nSpecies,nTs);
    
%    Xmean = reshape(mean(X,2),nSpecies,nTs);
%    Xstd  = reshape(std(X,0,2),nSpecies,nTs);
%    XmeanCI = zeros(nSpecies,nTs,2,'single');
%    XstdCI  = zeros(nSpecies,nTs,2,'single');

%    % confidence interval for the mean  (t-distribution)
%    c=tinv(1-(1-CI)/2,nMC-1);
%    XmeanCI(:,:,1) = c*Xstd/sqrt(nMC);
%    XmeanCI(:,:,2) = Xmean+XmeanCI(:,:,1);
%    XmeanCI(:,:,1) = Xmean-XmeanCI(:,:,1);
    
%    % confidence interval for the std dev 
%    %   (chi-distribution, assumes normal population)
%    XstdCI(:,:,1)=sqrt(nMC-1)*Xstd/sqrt(chi2inv(1-(1-CI)/2,nMC-1));
%    XstdCI(:,:,2)=sqrt(nMC-1)*Xstd/sqrt(chi2inv((1-CI)/2,nMC-1));
end

fprintf('finished %.2fsec\n',etime(clock,t0));

