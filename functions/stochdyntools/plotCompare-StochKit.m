function [exact,approx,MC]=plotCompare-StochKit(net,yPlots,xPlots,Tmax,varargin)

clf
exact=cell(0);
MC=struct;
approx=cell(0);
sampleMC=0;
distributionMC=0;

options=odeset('RelTol',1e-6);%,'AbsTol',1e-6);%,'Stats','on');
for thisSubPlot=1:nargin-4
  subplot(yPlots,xPlots,thisSubPlot)
  
  PlotTitle=varargin{thisSubPlot}{1}{1};
  Xminmax=varargin{thisSubPlot}{1}{2};
  
  if iscell(varargin{thisSubPlot}{2})
    % list of line plots 
    hndl=[];  
    leg={};    
    for thisLine=2:length(varargin{thisSubPlot})
      if length(varargin{thisSubPlot}{thisLine}) == 4

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
	
	% compute moment dynamics if needed	
	if length(exact)<maxdeg || ~isstruct(exact{maxdeg})
	  exact{maxdeg}.mdyn=momentDynamics(net,maxdeg);
	end
	% compute moment closure if needed	
	if length(approx)<maxdeg || ~isfield(approx{maxdeg},methodField)
	  closure.mdyn=momentClosure(net,exact{maxdeg}.mdyn,method);
	  dynamics2mfile(net,closure.mdyn,'fun',[])
	  if isfield(closure.mdyn,'dotPhi')
	    x0=[closure.mdyn.Mu.x0;closure.mdyn.Phi.x0];
	  else  
	    x0=closure.mdyn.Mu.x0;
	  end
	  [closure.t,closure.mu]=ode23s(@(t,x)fun(x),[0,Tmax],x0,options);
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
	grid on	
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
	if ~strcmp(varargin{thisSubPlot}{thisLine}{1},'sampleMC')
	  error('plotCompare: ''%s'' method without 3 arguments',...
	      varargin{thisSubPlot}{thisLine}{1}) 
	end	  

	%%%%%%%%%%%%%%%%%%%%%%	
	%% Monte Carlo plot
	%%%%%%%%%%%%%%%%%%%%%%	
	species=varargin{thisSubPlot}{thisLine}{2};
	style=varargin{thisSubPlot}{thisLine}{3};

	if ~sampleMC	
	  if length(exact)<1 || ~isstruct(exact{1})
	    exact{1}.mdyn=momentDynamics(net,1);
	  end
	  
	  if ~distributionMC
	    net2stochKit(net,'ProblemDefinition',exact{1}.mdyn.Mu.x0);
	  end	    
	  cmd='make single';
	  fprintf('executing: ''%s''...\n',cmd);
	  [rc,res]=system(cmd); 
	  if rc
	    disp(res)	    
	    error('\nFailed to make stochKit executable.\n');
	  end
	  sampleMC=1;
	  cmd=sprintf('rm -f single.txt; nice single %0g single.txt',Tmax);
	  fprintf('executing: ''%s''...\n',cmd);
	  [rc,res]=system(cmd);	  
	  if rc
	    disp(res)	    
	    error('\nFailed to run stochKit.\n');
	  end
	  MC.sample=load('single.txt');
	end  
	
	for i=1:length(net.species)
	  if strcmp(species,net.species(i).id)
	    hndl(1+3*(thisLine-2):3*(thisLine-1))=plot(...
		MC.sample(:,1),...
		MC.sample(:,1+i),style);
	    break
	  end
	end	    
	axis([0,Tmax,Xminmax])
	grid on	
	drawnow
	hold on	  
	leg{end+1}=sprintf('%s (MC)',...
	    mylatex(sym(species),exact{1}.mdyn.texrules));
      end		
    end      
    legend(hndl(1:3:end),leg,'location','best');
    title(PlotTitle)
    drawnow    
  else
    %%%%%%%%%%%%%%%%%%%%%%	
    %% distribution plot    
    %%%%%%%%%%%%%%%%%%%%%%	
    nMCs=varargin{thisSubPlot}{3};
    species=varargin{thisSubPlot}{4};

    if ~distributionMC
      if length(exact)<1 || ~isstruct(exact{1})
	exact{1}.mdyn=momentDynamics(net,1);
      end
      
      if ~sampleMC	
	net2stochKit(net,'ProblemDefinition',exact{1}.mdyn.Mu.x0);
      end
      cmd='make stats';
      fprintf('executing: ''%s''...\n',cmd);
      [rc,res]=system(cmd);	  
      if rc
	disp(res)
	error('\nFailed to make stochKit executable.\n');
      end
      distributionMC=1;
      cmd=sprintf('rm -f stats.txt; nice stats %0g %d stats.txt',Tmax,nMCs);
      fprintf('executing: ''%s''...\n',cmd);
      [rc,res]=system(cmd);
      if rc
	disp(res)
	error('\nFailed to run stochKit.\n');
      end
      MC.distribution=load('stats.txt');
    end      

    for i=1:length(net.species)
      if strcmp(species,net.species(i).id)
	[n,x]=hist(MC.distribution(:,i),min(20,floor(Xminmax(2)-Xminmax(1)+1)));
	break
      end
    end	    
    barh(x,n);
    axis([0,1.1*max(n),Xminmax])
    legend(sprintf('%s\\approx%g\\pm%g (MC,%d)',...
	mylatex(sym(species),exact{1}.mdyn.texrules),...
	smartround(mean(MC.distribution(:,i))),...
	smartround(std(MC.distribution(:,i))),nMCs),'location','best');
    grid on	
    drawnow    
    hold on
%    xx=min(x):(max(x)-min(x))/1000:max(x);
%    pdf=getDistribution(net,mdyn3,mu3(end,:),'X4',xx,'lognormal');
%    npdf=size(MC.distribution,1)*pdf*(x(2)-x(1));
%    plot(npdf,xx)

  end    
  
end  

function y=smartround(x)

ndigits=4;

scale=10^min(0,floor(log10(x))-ndigits+1);
y=scale*round(x/scale);