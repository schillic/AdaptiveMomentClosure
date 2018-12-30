function [exact,approx,MC]=plotCompare(net,yPlots,xPlots,Tmax,varargin)
% [exact,approx,MC]=plotCompare(net,yPlots,xPlots,Tmax,plot1,plot2,...)
% 
% Plots comparisons between alternative moment closure methods and
% moment estimates obtained by Monte Carlos simulations
% 
% Inputs
% ------
%
% net is a structure describing the network of chemical reactions.
%     It is typically obtained from a .net file using net=readNet(filename) 
%
% yPlots,xPlots provide the size of a grid where the plots will be displayed
%
% Tmax is the final simulation time for all plots
%
% plot1, plot2, ... are lists that describe each plot. The number
%      of plots should not exceed yPlots*xPlots (and generally
%      should be exactly equal to this number).
%
%      Each list has one of the two following forms
%      1) { {'plot title',maximum value for y-axis}, ...
%           'MCdistribution', # of MC simulations , 'species' }
%      
%         In this case an estimate of the distribution for the given
%         species at the final time will be ploted. This estimate is
%         obtained by running the given number of Monte Carlo
%         simulations.
%
%      2) { {'plot title',maximum value for y-axis}, ...
%           {'closure method 1',maxdeg,'symExpression_1','plot style 1'},...
%           {'closure method 2',maxdeg,'symExpression_2','plot style 2'},...
%           ... }
%
%         In this case, the mean and standard deviations of the
%         expressions given will be plotted. The means will be plotted
%         with the given style and dotted lines will be used to depict
%         the mean plus/minus one standard deviation. The plot styles
%         are are character strings defining the style of the line to
%         be used, as in Matlab's plot() command.
%
%         The function reuses (and saved) all previous computations so that
%         each new plot may not require recomputation of moment
%         closure or the solution to the corresponding ODE.
%  
%         This function uses getCMoments() so the reader is referred
%         to the documentation of that function for details on the
%         computation of means and standard deviations.
% 

% Output
% ------
%
% exact is an array contained the exact moment dynamics for
%       different degrees
%
% approx is an array containing the closed moment dynamcis for
%       different degrees and methods
%
% MC is and array containing the results of Monte Carlos runs 

clf
exact=cell(0);
MC=cell(0);
approx=cell(0);

options=odeset('RelTol',1e-6);%,'AbsTol',1e-6);%,'Stats','on');
for thisSubPlot=1:nargin-4
  subplot(yPlots,xPlots,thisSubPlot)
  
  PlotTitle=varargin{thisSubPlot}{1}{1};
  Xminmax=varargin{thisSubPlot}{1}{2};
  
  if length(varargin{thisSubPlot})>1
  if iscell(varargin{thisSubPlot}{2})
    % list of line plots 
    hndl=[];  
    leg={};    
    for thisLine=2:length(varargin{thisSubPlot})
      if length(varargin{thisSubPlot}{thisLine}) ~= 4
          error('plotCompare: cell array of length 4 expected to describe plot\n\t{method,maxdeg/#simulations,species,style}');
      end
      
      %%%%%%%%%%%%%%%%%%%%%%	
      %% Moment closure plot
      %%%%%%%%%%%%%%%%%%%%%%	
      method=varargin{thisSubPlot}{thisLine}{1};     
      maxdeg=varargin{thisSubPlot}{thisLine}{2};
      species=varargin{thisSubPlot}{thisLine}{3};
      style=varargin{thisSubPlot}{thisLine}{4};
      
      if iscell(method)
	  methodField=method{1};
      else
	  methodField=method;
      end
      
      if ~strcmp(methodField,'MCsample') && ~strcmp(methodField,'MCmean')
        % compute moment dynamics if needed	
        if length(exact)<maxdeg || ~isstruct(exact{maxdeg})
            exact{maxdeg}.mdyn=momentDynamics(net,maxdeg);
        end
        % compute moment closure if needed	
        if length(approx)<maxdeg || ~isfield(approx{maxdeg},methodField)
            closure.mdyn=momentClosure(net,exact{maxdeg}.mdyn,method);
            preamble=sprintf('%%%% created by plotCompare(maxdeg=%d,method=''%s'',funname=''%s'',symParameters=''%s'')\n\n',maxdeg,method,'fun','');
            dynamics2mfile(net,closure.mdyn,'fun',[],preamble)
            if isfield(closure.mdyn,'dotPhi')
                x0=[closure.mdyn.Mu.x0;closure.mdyn.Phi.x0];
            else  
                x0=closure.mdyn.Mu.x0;
            end
            t0=clock;
            fprintf('Solving ODE (%d order)... ',size(fun,1));
            [closure.t,closure.mu]=ode23s(@(t,x)fun(x),[0,Tmax],x0,options);
            fprintf('finished %.2fsec\n',etime(clock,t0));
            if length(approx)<maxdeg
                approx{maxdeg}=struct(methodField,closure);
            else
                approx{maxdeg}=setfield(approx{maxdeg},methodField,closure);
            end	    
        else
            closure=getfield(approx{maxdeg},methodField);	  
        end
        [hndl(1+3*(thisLine-2):3*(thisLine-1)),av,st]=plotCMoments(...
            net,closure.mdyn,closure.t,closure.mu,species,style);
        axis([0,Tmax,Xminmax])
        %        grid on	
        drawnow    
        hold on
        if isnan(st(end))	
            leg{end+1}=sprintf('%s\\approx%g (%s,%d)',...
                               mylatex(sym(species),closure.mdyn.texrules),...
                               smartround(av(end)),methodField,maxdeg);
        else
            leg{end+1}=sprintf('%s\\approx%g\\pm%g (%s,%d)',...
                               mylatex(sym(species),closure.mdyn.texrules),...
                               smartround(av(end)),smartround(st(end)),methodField,maxdeg);
        end	  
      else
        %%%%%%%%%%%%%%%
	%% Monte Carlo 
	%%%%%%%%%%%%%%%
        if length(net.raterule)>0
            error('cannot simulate continuous rates\n');
        end
        if length(exact)<1 || ~isstruct(exact{1})  % to get texrules
	    exact{1}.mdyn=momentDynamics(net,1);
        end
        
        nMC=maxdeg;
        
        thisMC=length(MC)+1;
        [MC{thisMC}.Q,MC{thisMC}.b,MC{thisMC}.c,MC{thisMC}.s,MC{thisMC}.X0]=...
            quadPropensities(net);
        MC{thisMC}.Q=double(subsParameters(net,MC{thisMC}.Q));
        MC{thisMC}.b=double(subsParameters(net,MC{thisMC}.b));
        MC{thisMC}.c=double(subsParameters(net,MC{thisMC}.c));
        MC{thisMC}.s=double(subsParameters(net,MC{thisMC}.s));
        
        MC{thisMC}.Ts=(0:Tmax/100:Tmax)';
        MC{thisMC}.species=sampledSSA(MC{thisMC}.Q,MC{thisMC}.b,MC{thisMC}.c,MC{thisMC}.s,MC{thisMC}.X0,nMC,MC{thisMC}.Ts);
	
	for i=1:length(net.species)
          if strcmp(species,net.species(i).id)
            if strcmp(method,'MCsample')
              hndl(1+3*(thisLine-2):3*(thisLine-1))=plot(MC{thisMC}.Ts,...
                                                         reshape(MC{thisMC}.species(i,1:min(end,3),:),[],length(MC{thisMC}.Ts))',style);
            else % MCmean
                data=reshape(MC{thisMC}.species(i,:,:),nMC,[])';
                av=mean(data,2);
                sd=std(data,0,2);
                
                % replace '--' or '-. or '-' by ':' (dotted) to mark +- one std. dev.
                stddev_style=style;
                k=regexp(stddev_style,'--|-.');
                if ~isempty(k)
                    stddev_style(k:k+1)=':';
                else
                    k=regexp(stddev_style,'-');
                    if ~isempty(k)
                        stddev_style(k)=':';
                    end
                end
                
                hndl(1+3*(thisLine-2):3*(thisLine-1))=plot(MC{thisMC}.Ts,...
                                                           av,style,MC{thisMC}.Ts,[av-sd,av+sd],stddev_style);
            end
            break
          end
        end	    
        axis([0,Tmax,Xminmax])
        %	grid on	
	drawnow
	hold on	  
        if strcmp(method,'MCsample')
          leg{end+1}=sprintf('%s (MC sample)',...
                             mylatex(sym(species),exact{1}.mdyn.texrules));
        else % MCmean
          leg{end+1}=sprintf('%s\\approx%g\\pm%g (%d MC)',...
                             mylatex(sym(species),exact{1}.mdyn.texrules),...
                             smartround(av(end)),smartround(sd(end)),nMC);
        end
      end		
    end      
    legend(hndl(1:3:end),leg,'location','best');
  else % if iscell(varargin{thisSubPlot}{2})
      if ~strcmp(varargin{thisSubPlot}{2},'MCdistribution')
          error('plotCompare: ''MCdistribution'' expected instead of ''%s''',...
                varargin{thisSubPlot}{2}) 
      end	  
      
    %%%%%%%%%%%%%%%%%%%%%%	
    %% distribution plot    
    %%%%%%%%%%%%%%%%%%%%%%	
    nMCs=varargin{thisSubPlot}{3};
    species=varargin{thisSubPlot}{4};

    if length(net.raterule)>0
        error('cannot simulate continuous rates\n');
    end
    
    if length(exact)<1 || ~isstruct(exact{1})  % to get texrules
        exact{1}.mdyn=momentDynamics(net,1);
    end
    
    thisMC=length(MC)+1;
    [MC{thisMC}.Q,MC{thisMC}.b,MC{thisMC}.c,MC{thisMC}.s,MC{thisMC}.X0]=quadPropensities(net);
    MC{thisMC}.Q=double(subsParameters(net,MC{thisMC}.Q));
    MC{thisMC}.b=double(subsParameters(net,MC{thisMC}.b));
    MC{thisMC}.c=double(subsParameters(net,MC{thisMC}.c));
    MC{thisMC}.s=double(subsParameters(net,MC{thisMC}.s));
    
    MC{thisMC}.distribution=sampledSSA(MC{thisMC}.Q,MC{thisMC}.b,MC{thisMC}.c,MC{thisMC}.s,MC{thisMC}.X0,nMCs,Tmax);
    MC{thisMC}.distribution=reshape(MC{thisMC}.distribution,length(net.species),nMCs);

    for i=1:length(net.species)
      if strcmp(species,net.species(i).id)
	[n,x]=hist(MC{thisMC}.distribution(i,:),min(20,floor(Xminmax(2)-Xminmax(1)+1)));
	break
      end
    end	    
    barh(x,n);
    axis([0,1.1*max(n),Xminmax])
    legend(sprintf('%s\\approx%g\\pm%g (%d MC)',...
	mylatex(sym(species),exact{1}.mdyn.texrules),...
	smartround(mean(MC{thisMC}.distribution(i,:))),...
	smartround(std(MC{thisMC}.distribution(i,:))),nMCs),'location','best');
    %    grid on	
    hold on
%    xx=min(x):(max(x)-min(x))/1000:max(x);
%    pdf=getDistribution(net,mdyn3,mu3(end,:),'X4',xx,'lognormal');
%    npdf=size(MC{thisMC}.distribution,2)*pdf*(x(2)-x(1));
%    plot(npdf,xx)
  end % if iscell(varargin{thisSubPlot}{2})
  else % if length(varargin{thisSubPlot})>1
      axis([0,Tmax,Xminmax])
  end
  title(PlotTitle)
  drawnow    
end  % for

function y=smartround(x)

ndigits=4;

scale=10^min(0,floor(log10(x))-ndigits+1);
y=scale*round(x/scale);