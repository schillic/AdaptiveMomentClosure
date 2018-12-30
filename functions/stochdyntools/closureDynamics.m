function mdyn=closureDynamics(net,maxdeg,method,funname,symParameters, VOLUME)
% mdyn=closureDynamics(net,maxdeg,method,funname,symParameters)
%
% This function computes the exact (open) moment dynamics
%           \dot Mu = A Mu + B barMu       ($)
%     where 
%       Mu    - column vectors containing the moments of interest
%       A,B   - appropriate matrices
%       barMu - column vector with higher-order moments
% and then uses the moment closure technique specified by the input
% parameter 'method' to compute approximate (closed) moment dynamics
% of the form 
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
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% maxdeg is an integer that specifies the largest degree for the
%        uncentered moments in Mu. If 'maxdeg' is not an integer, then
%        boolean variables are not taken into account for the degree
%        of a moment. In this case, more moments are included in Mu
%        for the same value of 'maxdeg'.
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
%         variance for the populations, i.e., assuming that the
%         equality (2) above holds for all species.
%
%         More details are provided in the PDF documentation.
%
% funname (optional) is a string containing the filename used to
%         create an m-file that computes the left-hand side of the
%         differential equations in ($$) or ($$$). For approximate
%         dynamics of the form ($$$), the Mu components of the state
%         appear above the Phi components. When the function created
%         by closureDynamics() is called with no arguments, it outputs
%         the names of the state variables (as a string array).
%  
%         The function created can be used as input to ode23s() to
%         simulate the moment dynamics (see PDF documentation).
%  
% symParameters (optional) is a vector of symbolic variables
%         corresponding to parameters that appear in the rate
%         expressions (e.g., rate constants). This vector is used only
%         in the creation of the function 'funname' and it allows
%         inclusion of additional inputs to this function that may be
%         needed to specify numerical values for parameters that
%         appear in the moment dynamics. These values will override
%         any default values specified in net.parameter 
%
% Output
% ------
%
% mdyn is a structure characterizing the exact moment dynamics in ($)
%      and the approximate moment dynamics in ($$) or ($$$). 
%      Details can be obtained with 'help mdyn'
%     
% Under the hood
% --------------
%
% The function closureDynamics() essentially combines the
% functionality of the three lower-level functions: momentDynamics(),
% momentClosure, and sym2mfile(). Documentation for these functions is
% provided in the corresponding m-files.
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
% November 4, 2006
% Added support for net.parameter
%
% May 7, 2007
% Added initial condition as output
%
% May 12, 2007
% Added 'zerocumulants' and 'lowdispersion' closures 
%
% June 27, 2007
% Returns everything in a single structure
% Returns texrules
%
% December 25, 2007
% Added 'lna' closure
%
% January 10, 2007
% Updated inline documentation

if nargin<3
  method='derMatch';
end

mdyn=momentDynamics(net,maxdeg);

% Christian: Fixed bug:
% The variable 'method' had to be a cell here for lna closure, but a string for
% the comment printing below.
% As this is complicated to pass anyway, I now pass the second entry as yet
% another argument to the function.
if (exist('VOLUME', 'var'))
    methodCell = {method, VOLUME};
else
    methodCell = {method};
end
mdyn=momentClosure(net,mdyn,methodCell);

if nargin>3
  if nargin<5
    symParameters=[];
  end

  preamble=sprintf('%%%% created by closureDynamics(net,maxdeg=%d,method=''%s'',funname=''%s'',symParameters=''%s'')\n\n',maxdeg,method,funname,char(symParameters));

  dynamics2mfile(net,mdyn,funname,symParameters,preamble);
end
