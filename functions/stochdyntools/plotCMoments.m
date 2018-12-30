function [handles,averages,stddevs]=...
      plotCMoments(net,mdyn,t,mu,varargin)
% [handles,average,stddev]=plotCMoments(net,mdyn,t,mu,...
%             symExpression_1,style_1,symExpression_2,style_2,...,assumezero)
% 
% Given a time series mu of uncentered moments, typically computed
% using ode23s(), this function computes and plots the mean and
% standard deviation of a given expression involving the species populations.
% 
% Internally, this function uses getCMoments() to compute means and standard
% deviations.
% 
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% mdyn is a structure characterizing the moment dynamics, typically
%      computed using closureDynamics(), momentDynamics(), or momentClosure(). 
%      Details can be obtained with 'help mdyn'
%
% mu is a matrix containing a time series of Mu or [Mu;Phi], typically
%    computed using ode23s(). Each row of mu corresponds to a
%    different time instant and each column to the different moment. 
%
% symExpression_1, symExpression_2,... are symbolic expressions whose
%      means and standard deviations will be plotted. The means of
%      each symExpression_i will be plotted with the style defined by
%      style_i and dotted lines will be used to depict the mean
%      plus/minus one standard deviation. 
%  
%      This function uses getCMoments() so the reader is referred to
%      the documentation of that function for details on the
%      computation of means and standard deviations. 
% 
% style_1,style_2,... are character strings defining the style of the
%      line to be used, as in Matlab's plot() command.
%
% assumezero is an optional boolean variable that, when nonzero, 
%    instructs the function to assume zero all moments that are not
%    available to compute the mean and/or standard deviations.
%
% Output
% ------
%
% handles is a vector with the handles of all the lines plotted.
%         These handles are useful, for example to produce legends.
%
% average is a time series with the expected values of all
%         expressions, with one time instant per row and one
%         expression per column. 
%
% stddev is a time series with the standard deviations of all
%        expressions, with one time instant per row and one expression
%        per column. 
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
% Created May 8, 2007 
%
% June 27, 2007
% mdyn passed as a single structure
%
% December 26, 2007
% made compatible with 'lna' moment closure
%
% January 10, 2007
% updated inline documentation
%
% July 18, 2009
% added assumezero input

if mod(nargin,2)
    assumezero=varargin{end};
    varargin(end)=[];
else
    assumezero=0;
end

handles=[];
oldhold=ishold;
averages=[];
stddevs=[];
for i=1:2:length(varargin)
  [average,stddev]=getCMoments(net,mdyn,mu,varargin{i},assumezero);
  averages=[averages,average];  
  stddevs=[stddevs,stddev];
  
  % replace '--' or '-. or '-' by ':' (dotted) to mark +- one std. dev.
  stddev_style=regexprep(varargin{i+1},'--|-.|:|-',':');

  handles=[handles;plot(t,averages(:,end),varargin{i+1},t,...
			[averages(:,end)-stddevs(:,end),averages(:,end)+stddevs(:,end)], ...
			stddev_style)];
  hold on
end

% return to previous hold state
if ~oldhold
  hold off
end
