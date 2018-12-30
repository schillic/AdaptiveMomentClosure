%time horizon
Tmax = 500;
Ts = 1;
tIn = 0:Ts:Tmax;

%reaction parameters
a = 2*10^(-2);
b = 2*10^(-3);
c = 2*10^(-4);
d = 2*10^(-1);
e = 2*10^(-4);
f = 2*10^(-2);

%stoichiometry matrix
S = [1,0,0;
     -1,1,0;
     1,-1,0;
     0,0,1;
     -1,0,0;
     0,0,-1];

%number of simulated trajectories
n=1000;

%X (Transcription Factor)
Xa=zeros(n,Tmax/Ts+1);
%Y (Gene)
Ya=zeros(n,Tmax/Ts+1);
%Z (Protein)
Za=zeros(n,Tmax/Ts+1);

X0 = 0;
Y0 = 0;
Z0 = 0;

z0 = [0,0,0];
    
for k=1:n

    if mod(k,100) == 0
        disp(k)
    end

    T=Ts;
    X=zeros(1,Tmax/Ts+1);
    Y=zeros(1,Tmax/Ts+1);
    Z=zeros(1,Tmax/Ts+1);
    t=0;

    z=z0;

    while t<Tmax
      
        
      %reaction propensities
      h1=a;
      h2=b*z(1)*(1-z(2));
      h3=c*z(2);
      h4=d*z(2);
      h5=e*z(1)*z(3);
      h6=f*z(3);
      h0=h1+h2+h3+h4+h5+h6;	
      
      %simulate the waiting time to the next reaction
      if h0 ~= 0
        t1 = -log(rand)/h0;
        t = t+t1;
      else t = Tmax;
      end

      %save the state at which the trajectory is before the next time grid point is crossed
      while t >= T
          if t >= T && T <= Tmax
              ind = round(T/Ts + 1);
              X(ind)=z(1);
              Y(ind)=z(2);
              Z(ind)=z(3);
              T = T + Ts;
          elseif T > Tmax
              break
          end
      end


    %simulate which reaction fires next
    p=[h1/h0,h2/h0,h3/h0,h4/h0,h5/h0,h6/h0];
    uni=rand(1);
    cumprob=[0 cumsum(p)];
    ind=max(find((uni>cumprob)));
    %update the state using the stoichiometry matrix S
    z = z + S(ind,1:end);
    
    end
    
    X(1) = z0(1);
    Y(1) = z0(2);
    Z(1) = z0(3);

    %save the trajectory
    Xa(k,:) = X;
    Ya(k,:) = Y;
    Za(k,:) = Z;
    
    %plot individual trajectories
%     figure(1)
%     hold on
%     plot(Z)
    
end

%compute means of the simulated trajectories  
MX1 = mean(Xa,1);
MY1 = mean(Ya,1);
MZ1 = mean(Za,1);

%compute variances of the simulated trajectories
VarX = var(Xa,1);
VarY = var(Ya,1);
VarZ = var(Za,1);

%compute covariances of the simulated trajectories
CovXY = mean(Xa.*Ya,1) - mean(Xa).*mean(Ya);
CovXZ = mean(Xa.*Za,1) - mean(Xa).*mean(Za);
CovYZ = mean(Ya.*Za,1) - mean(Ya).*mean(Za);

color = 'b';

figure(1)
hold on
plot(tIn,MX1,color)

figure(2)
hold on
plot(tIn,MY1,color)

figure(3)
hold on
plot(tIn,MZ1,color)

figure(4)
hold on
plot(tIn,VarX,color)

figure(5)
hold on
plot(tIn,VarY,color)

figure(6)
hold on
plot(tIn,VarZ,color)



