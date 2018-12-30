function mdyn=momentClosure(net,mdyn,method)
% mdyn=momentClosure(net,mdyn,method)
%
% This function takes exact (open) moment dynamics of the form
%           \dot Mu = A Mu + B barMu       ($)
% where 
%       Mu    - column vectors containing the moments of interest
%       A,B   - appropriate matrices
%       barMu - column vector with higher-order moments
% computed by momentdynamics() and uses the moment closure technique
% specified by the input parameter 'method' to compute approximate
% (closed) moment dynamics of the form 
%           \dot Mu = A Mu + B F(\Mu)      ($$)
%     where 
%       Mu    - column vector that approximates Mu
%       F()   - moment closure function
% or 
%           \dot Phi = G(Phi)
%           \dot Mu  = A Mu + B F(Phi,Mu)  ($$$)
%     where 
%       Mu    - column vector that approximates Mu
%       F()   - moment closure function
%       Phi   - state of an additional dynamical system
% 
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% mdyn is a structure characterizing the exact moment dynamics in ($),
%      typically computed by momentDynamics().  
%      Details can be obtained with 'help mdyn'
%
% method (optional) is a string that specifies the method that should
%        be used for moment closure.  This parameter is optional, and
%        in its absence the first method below is used. The following
%        methods are currently recognized: 
%    'dm' or 'derivativematching' --- moment closure obtained by
%         matching the derivatives of the exact and approximate moment
%         dynamics. 
%    'zc' or 'zerocumulants' --- moment closure obtained by assuming
%         that the high-order cumulants corresponding to all unknown
%         moments are equal to zero. 
%    'ld' or 'lowdispersion' --- moment closure obtained by assuming
%         that the high-order normalized centered moments
%         corresponding to all unknown moments are equal to zero.
%         For second order closure this is the same as 'zc'
%    'qd' or 'quasideterministic' --- moment closure obtained by
%         assuming that the high-order quasi-deterministic normalized
%         centered moments corresponding to all unknown moments are
%         equal to zero. For this approximation, the closed dynamics
%         are of the form ($$$) where Phi is the solution to the
%         determinstic dynamics (in molecule counts).  
%    {'lna','Volume'} or {'linearnoiseapproximation','Volume'} or
%    {'vankampen','Volume'} --- Van Kampen's linear noise
%         approximation. The second element of the cell specifies the
%         variable to be used as the volume. This variable should have
%         be declared in the 'parameters' section of the .net file.
%         To achieve moment closure, one considers the limit as this
%         variable converges to infinity. For this limit to result in
%         a closed set of moment equations, the all reactions should
%         be elementary and their rates should depend on the volume as
%         follows:  
%             0   -> *       rate = c Volume
%             X   -> *       rate = c X
%            2 X  -> *       rate = c X (X-1)/ Volume
%           X + Y -> *       rate = c X Y / Volume
%         For this approximation, the closed dynamics are of the form
%         ($$$), where Phi is the solution to the determinstic
%         dynamics (in concentrations). 
%    'zv' or 'zerovariance' --- moment closure assuming a negligible
%         variance for the populations, i.e., assuming that for any
%         species X1
%             E[ X1 f(X2,...Xk) ] = E[ X1 ] E[ f(X2,...Xk) ]
%         for any function $f(.)$ of the remaining species.
%
%         More details are provided in the PDF documentation.
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
% October 27, 2006
% Added support for multiple truncation methods 
%
% May 12, 2007
% Added 'zerocumulants' and 'lowdispersion' closures 
%
% June 27, 2007
% Returns everything in a single structure
%
% December 25, 2007
% Added 'lna' closure
% December 31, 2007
% Added 'quasidet' closure
% Allows shortened names
%
% January 10, 2007
% Updated inline documentation

if nargin<3
  method='derMatch';
end

if ~iscell(method)
  method={method};  
end

fprintf('computing moment closure functions (''%s'')... ',method{1});
t0=clock;

computeDotMu=1;
switch lower(method{1})
 case {'dm','dermatch','derivativematching'}
    mdyn.approxDotMu.barMu=derMatchClosure(net,mdyn);
 case {'zv','zerovar','zerovariance'}
    mdyn.approxDotMu.barMu=zeroVarianceClosure(net,mdyn);
 case {'zc','zerocum','zerocumulants'}
    mdyn.approxDotMu.barMu=zeroCumulantsClosure(net,mdyn);
 case {'ld','lowdisp','lowdispersion'}
    mdyn.approxDotMu.barMu=lowDispersionClosure(net,mdyn);
 case {'qd','quasidet','quasideterministic'}
    mdyn=quasiDeterministicClosure(net,mdyn);
    computeDotMu=0;   
 case {'lna' ,'linearnoiseapproximation','vankampen'}
    if (length(method)<2)
      error('momentClosure: volume variable not specified for ''%s'' method',method{1});
    end      
    mdyn=lnaClosure(net,mdyn,sym(method{2}));
    computeDotMu=0;   
 otherwise
    error('momentClosure: method ''%s'' not implemented',method{1})
end

% compute approxDotMu
if computeDotMu
  if isempty(mdyn.dotMu.B)
    mdyn.approxDotMu.sym=simplify(mdyn.dotMu.A*mdyn.Mu.sym);
  else
    mdyn.approxDotMu.sym=simplify(mdyn.dotMu.A*mdyn.Mu.sym...
	+mdyn.dotMu.B*mdyn.approxDotMu.barMu);
  end
end

fprintf('finished %.2fsec\n',etime(clock,t0));

% compute approxDotMu.A, approxDotMu.C
mdyn.approxDotMu.A=jacobian(mdyn.approxDotMu.sym,mdyn.Mu.sym);
mdyn.approxDotMu.c=expand(mdyn.approxDotMu.sym-mdyn.approxDotMu.A*mdyn.Mu.sym);
