function [Xmean,Xstd,XmeanCI1,XmeanCI2,XstdCI1,XstdCI2] = meanStdCI(X,CI,dim)
% [Xmean,Xstd,XmeanCI1,XmeanCI2,XstdCI1,XstdCI2] = meanStdCI(X,CI,dim)
% 
% Estimates the mean and standard devition, providing confidence
% intervals 
%%
% Inputs:
%
% X  - Vector or n-dimensional Matrix with data. 
%      1) When vector returns a scalar mean and standard devition
%      2) When a matrix returns a vector of means and standard
%         deviations, taken aong the given dimension
%
% CI  - desidered confidence interval
%
% dim - (optional) dimension along which the mean is take for a matrix data
%       input. Default = 1;

%
% Outputs:
%
% Xmean - scalar/matrix with means
%
% Xstd  - scalar/matrix with standard deviations
%
% XmeanCI1,XmeanCI2 - scalars/matrices with the confidence interval
%                     for the means 
%
% XstdCI1,XstdCI2 - scalars/matrices with the confidence intervals
%                   for the standard deviations (assumes a normal distribution)
%
% Copyright (C) 2008  Joao Hespanha

% somewhat replicates stats/normfit

    
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
% Created January 18, 2009

if nargin<3
    dim=1;
end

% mean and std. dev.

N=size(X,dim);
Xmean = mean(X,dim);
Xstd  = std(X,0,dim);

% confidence interval for the mean  (t-distribution)
c=tinv(1-(1-CI)/2,N-1);
XmeanCI1 = c*Xstd/sqrt(N);
XmeanCI2 = Xmean+XmeanCI1;
XmeanCI1 = Xmean-XmeanCI1;

% confidence interval for the std dev 
%   (chi-distribution, assumes normal population)
XstdCI1=sqrt(N-1)*Xstd/sqrt(chi2inv(1-(1-CI)/2,N-1));
XstdCI2=sqrt(N-1)*Xstd/sqrt(chi2inv((1-CI)/2,N-1));
