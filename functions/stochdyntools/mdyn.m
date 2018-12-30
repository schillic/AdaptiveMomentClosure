% The mdyn structure characterizes 
% (1)  the exact moment dynamics
%           \dot Mu = A Mu + B barMu       ($)
%      where 
%        Mu    - column vectors containing the moments of interest
%        A,B   - appropriate matrices
%        barMu - column vector with higher-order moments
% (2a) approximate moment dynamics of the form
%           \dot Mu = A Mu + B F(\Mu)      ($$)
%      where 
%        Mu    - column vector that approximates Mu
%        F()   - moment closure function
% or 
% (2b) approximate moment dynamics of the form
%           \dot Phi = G(Phi)              ($$$ a)
%           \dot Mu  = A Mu + B F(Phi,Mu)  ($$$ b)
%      where 
%        Mu    - column vector that approximates Mu
%        F()   - moment closure function
%        Phi   - state of an additional dynamical system
% 
% The mdyn structure includes:
% . the entries of the Mu vector in the exact and approximate
%   moment dynamics ('Mu' field),
% . the exact (open) moment dynamics in ($) ('dotMu' field),
% . the approximate (closed) moment dynamics in ($$) ('approxDotMu'
%   field) or in ($$$) ('approxDotMu' and 'dotPhi' fields), and
% . set of formatting rules to produce the moment dynamics in LaTeX
%   format ('texrules' field).
%
% A detailed explanation of the key fields of this structure
% follows. This structure may contain additional fields not described
% here, which are used internally by the toolbox. 
%
% Exact moment dynamics
% ---------------------
%
% The following fields describe the exact moment dynamics in ($):
%
% 'Mu' is a structure describing the entries of the vector Mu in ($), ($$), and
% ($$$ b). The following table describes the key fields of this
% structure for a network of chemical reactions with species X1,X2,X3,...:
% 
%            |                                    | value of row corresponding
%  field     | type                               | to moment 
%            |                                    | E[X1^{m1}X2^{m2} X3^{m3}...]
%  ----------+------------------------------------+-----------------------------
%  Mu.sym    | col vector of symbolic variables   | mu_X1m1_X2m2_X3m3...
%  Mu.mon    | col vector of symbolic expressions | X1^m1*X2^m2*X3^m3...
%  Mu.ndx    | matrix of integers with one column | [m1,m2,m3,...]
%            | per species (in the order they were|
%            | declared in the .net file)         | 
%  Mu.x0     | col vector of doubles              | initial value from .net file
%
% 'barMu' is a structure describing the entries of the vector barMu in ($).
% The following table describes the key fields of this structure for a
% network of chemical reactions with species X1,X2,X3,...:
%            |                                    | value of row corresponding
%  field     | type                               | to moment 
%            |                                    | E[X1^{m1}X2^{m2} X3^{m3}...]
%  ----------+------------------------------------+-----------------------------
%  barMu.sym | col vector of symbolic variables   | mu_X1m1_X2m2_X3m3...
%  barMu.mon | col vector of symbolic expressions | X1^m1*X2^m2*X3^m3...
%  barMu.ndx | matrix of integers with one column | [m1,m2,m3,...]
%            | per species (in the order they were|
%            | declared in the .net file)         | 
%
% 'dotMu' is a structure describing the exact moment dynamics ($). It
% contains the following fields: 
%  field     | type                               | value
%  ----------+------------------------------------+-----------------------------
%  dotMu.A   | matrix of symbolic expressions     | matrix A in ($)
%  dotMu.B   | matrix of symbolic expressions     | matrix B in ($)
%  dotMu.sym | col vector of symbolic expressions | whole right-hand side of ($)
%
%  The right-hand side of ($) can thus be obtained using any one of
%  the following two symbolic expressions:
%    . dotMu.A * Mu.sym + dotMu.B * barMu.sym
%    . dotMu.sym
%
% Approximate moment dynamics
% ---------------------------
%
% The following fields describe the approximate moment dynamics in
% ($$) or ($$$ b), depending on the moment closure technique used. It
% contains the following fields: 
%  field             | type                    | value
%  ------------------+-------------------------+---------------------
%  approxDotMu.barMu | col vec. of symb.expres.| approx.value for barMu in ($)
%  approxDotMu.sym   | col vec. of symb.expres.| right-hand side of ($$)/($$$ b)
%  approxDotMu.A     | matrix of symb. expres. | Jacobian of approxDotMu.sym 
%                    |                         | with respect to Mu.sym
%  approxDotMu.c     | vector of symb. expres. | approxDotMu.sym
%                    |                         |    -approxDotMu.A*Mu.sym
%
% The right-hand side of ($$) or ($$$ b) can thus be obtained using
% any one of the following three symbolic expressions:
%   . dotMu.A * Mu.sym + dotMu.B * approxDotMu.barMu
%   . approxDotMu.sym
%   . approxDotMu.A * Mu.sym + approxDotMu.c
% The last representation is especially useful for the 'qd' and 'lna'
% moment closure methods because in this case the matrices approxDotMu.A and
% approxDotMu.c do not depend on the entries of Mu.sym
%
% Currently, the moment closure method 'lna' does not return the field
% approxDotMu.barMu
%
% The following fields are only needed to describe the approximate
% moment dynamics in ($$$). 
%
% 'Phi' is a structure describing the entries of the vector Phi in
% ($$$ a). The following table describes the key fields of this
% structure for a network of chemical reactions with species X1,X2,X3,...:
%            |                                    | value of row corresponding
%  field     | type                               | to moment 
%            |                                    | E[X1^{m1}X2^{m2} X3^{m3}...]
%  ----------+------------------------------------+-----------------------------
%  Phi.sym   | col vector of symbolic variables   | phi_X1m1_X2m2_X3m3...
%  Phi.mon   | col vector of symbolic expressions | X1^m1*X2^m2*X3^m3...
%            |                                    |    ('qd' method)    
%            |                                    | (X1/V)^m1(X2/V)^m2...
%            |                                    |    ('lna' method)    
%  Phi.x0    | col vector of doubles              | initial value from .net file
%
% 'dotPhi' is a vector of symbolic expressions with the derivative of
% Phi as in ($$$ a) 
%
% The following additional fields are only returned by the 'lna'
% moment closure method: 
%
% 'Xi' is a structure describing the entries of the first-order
% "perturbation" vector Xi in for the linear noise approximation
% (see PDF documentation).
%
% 'dotXi' is a symbolic vector with the exact derivative of Xi.
%
% 'approxDotXi' is a structure describing the approximate derivative
% of Xi, obtained by making the volume converge to infinity. It
% contains the following fields: 
%  field             | type                    | value
%  ------------------+-------------------------+---------------------
%  approxDotXi.sym   | col vec. of symb.expres.| approximate derivative of Xi
%  approxDotXi.A     | matrix of symb. expres. | Jacobian of approxDotXi.sym 
%                    |                         | with respect to Xi.sym
%  approxDotXi.c     | vector of symb. expres. | approxDotXi.sym
%                    |                         |    -approxDotXi.A*Xi.sym
%
% LaTeX formatting
% ----------------
%
% The following field is used to produced LaTeX-formatted versions of
% any symbolic expression involving reaction parameters, moments, or
% exact and approximate moment dynamics.  These rules should be used
% by the function mylatex() described.
%
% 'texrules' is a cell array of strings with one row per
% transformation rule and two columns: 
%   . the first column contains a regular expression, following the
%     syntax of Matlab's regexp() 
%   . the second column contains the string to replace the regular
%     expression. 

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
% Created January 9, 2008

