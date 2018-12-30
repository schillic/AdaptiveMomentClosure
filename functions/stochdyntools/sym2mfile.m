function sym2mfile(filename,symFunction,symVariables,symParameters,preamble)
% sym2mfile(filename,symFunction,symVariables,symParameters,preamble)
%
% This function creates an m-file that computes a given
% vector-valued symbolic function. When the created function is
% called with no arguments, it outputs the names of the symbolic
% input variables (as a string array).
%
% Input
% -----
%
% filename : filename for the m-file to be created
%
% symFunction : vector of symbolic expressions with the function to
%               be computed by the m-file.
%
% symVariables : vector of symbolic variables to be passed as the
%                first input to the m-file
%
% symParameters : vector of symbolic variables to be passed as
%                 subsequent inputs to the m-file
%
% preamble: (optional) string to be included in the m-file before
%           any computations. Can be used, e.g., to set parameter
%           values.
%
% Output
% ------
%
% No output returned, but a m-file is created. 
%
% ATTENTION: any existing m-file with the given name will be overwritten.
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
% Added support preamble
%
% February 2, 2009
% Added 'rehash'
%
% February 27, 2010
% only clear function just created
%
% June 12, 2009
% corrected code when nargin<5

if nargin<4
  symParameters=sym([]);
end

if nargin<5
  preamble='';
end

fid=fopen(sprintf('%s.m',filename),'w');

%% check types
try
  for i=1:length(symVariables)
    eval(sprintf('%s=1;',char(symVariables(i))));
  end
catch
  error('symVariables(%d) ''%s'' is not a symbolic variable',i,char(symVariables(i)))
end
try
  for i=1:length(symParameters)
    eval(sprintf('%s=1;',char(symParameters(i))));
  end
catch
  error('symParameters(%d) ''%s'' is not a symbolic variable',i,char(symParameters(i)))
end

%% produce function header
fprintf(fid,'function [dx,magic]=%s(x',filename);
for i=1:length(symParameters)
  fprintf(fid,',%s',char(symParameters(i)));  
end  

%% include preamble
fprintf(fid,')\n\n%s\n',preamble);

magic=round(100000*rand(1));
%% return variable names when nargin==0
fprintf(fid,'if nargin==0\n  %% return variable names when nargin==0\n  syms ',filename);
for i=1:length(symVariables)-1
  fprintf(fid,'%s ',char(symVariables(i)));  
end  
fprintf(fid,'%s\n  dx={',char(symVariables(end)));
for i=1:length(symVariables)-1
  fprintf(fid,'''%s'';',char(symVariables(i)));  
end  
fprintf(fid,'''%s''};\n  magic=%d;  %% random number used to check if correct function is being called\nelse\n',char(symVariables(end)),magic);  

%% compute function
fprintf(fid,'  %% compute function\n');
new=[];
old=[];
for i=1:length(symVariables)
  old=[old;cellstr(char(symVariables(i)))];  
  new=[new;cellstr(sprintf('x(%d)',i))];
end
%old
%new
dx=mysubs(symFunction,old,new,1);

for i=1:length(symFunction)
  fprintf(fid,'  dx(%d,1) = %s;\n',i,char(dx(i)));
end  
fprintf(fid, 'end\n');

fclose(fid);

% update list of known files and check timestamps, 
% does NOT replace clear functions
rehash

% to make sure that the function just created is indeed used
%clear functions %% appears to create trouble with 2009b
clear(filename)

[dummy,magicrc]=eval(filename);

if magic~=magicrc
    error('sym2mfile: function is not returning the expected value');
end