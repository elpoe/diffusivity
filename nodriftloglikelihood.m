function [out]=hmmloggennn(obs,k,dt,Dmin,Dmax)


%obs=load('randomdata.txt');

states=10; %number of states
N=length(obs); %number of steps
P0=zeros(states,1);
P0(:)=deal(1/states); %uniform initial state

%generate matrix with rates by shifting the row vector k 
kvec=zeros(1,states);
kvec(1)=-2*k;
kvec(2)=k;
kvec(end)=k;
kalt=zeros(1,states);
kalt(1)=-k;
kalt(2)=k;

K=zeros(states,states);
for i=1:states
    if i==1
       K(i,:)=kalt;
    elseif i==states
       kalt2=-1*circshift(kalt,states-2);
       K(i,:)=kalt2;
    else 
    K(i,:)=circshift(kvec,i-1);
    end
end


%%

TRANS=abs(expm(K*dt)); %transitionmatrix

Dvec=linspace(Dmin,Dmax,states)'; %diffusionskonstanterne er lige fordelt mellem minimum og maximumsværdien

z=4*dt;
z1=pi*z;
gen=@(x,D) -dot(x,x)./(z.*D)-log(z1.*D); %general 2D brownian motion without drift
genD=@(x) gen(x,Dvec(:)); %general 2D brownian motion w/o drift, with the different diffusion values as input

TRNSLOG=log(TRANS'); %log to transposed transition matrix

frwrd=log(zeros(states,N+1)); %frwrd will contain all calculations for each time step

for i=1:states
    frwrd(i,1)=log(P0(i));
end

for n=2:N+1
    for i=1:states
      
        %frwrd(i,n)=util_logsumexp(frwrd(:,n-1)+log(TRANS(i,:)')+gen(obs(n-1,:),Dvec(:))); 
        frwrd(i,n)=util_logsumexp(frwrd(:,n-1)+TRNSLOG(:,i)+genD(obs(n-1,:)));
 
    end
end
out=util_logsumexp(frwrd(:,N+1));
end
