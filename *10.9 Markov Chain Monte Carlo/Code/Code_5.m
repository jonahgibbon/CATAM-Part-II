digits(8)
Data={"Team","2024","2023","2022","2021","2020"
"Arsnal",83,90,78,87,81
"Asten Villa",47,56,45,50,60
"Blackborn Rovers",42,44,60,46,56
"Boltin Wandrers",58,53,44,40,60
"Charlston Athletic",46,53,49,44,45
"Chelsea Buns",95,79,67,64,71
"Evraton",61,39,59,43,46
"Fullem",44,52,48,44,53
"Livurpule",58,60,64,80,54
"Manchester Ununited",77,75,83,77,70
"Middlesbro",55,48,49,45,51
"Newcassel Divided",44,56,69,71,64
"Slothampton",32,47,52,45,55
"Tottenham Coldspur",52,45,50,50,59};
K=size(Data,1)-1;
T=size(Data,2)-1;

sigma_0=10;
alpha_0=10^-5;
beta_0=10^-3;
mu_0=60;
tau_0=20;


List=[20,30,50,100,150,250];
NCounter=1;
MeanProb=zeros(K,size(List,2));
MeanMu=zeros(K,size(List,2));
VarProb=zeros(K,size(List,2));
VarMu=zeros(K,size(List,2));

while NCounter<=size(List,2)
N=List(1,NCounter);
Number_of_Runs=500;

Prob=zeros(K,Number_of_Runs);
Mu=zeros(K,Number_of_Runs);

SuperCounter=1;
while SuperCounter<=Number_of_Runs

LineLength=fprintf('List: %d, Run: %.1f%% complete',N,SuperCounter*100/Number_of_Runs);
Sim=Gibbs(Data,sigma_0,alpha_0,beta_0,mu_0,tau_0,N);

Counter=2;
Carry=zeros(K,1);
while Counter<=N+1
    SubCounter=1;
    while SubCounter<=K
        if Sim(SubCounter,Counter)>Sim(2*K+1,Counter)
            Carry(SubCounter,1)=Carry(SubCounter,1)+1;
        end
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end
Prob(:,SuperCounter)=Carry/N;
Carry=sum(Sim,2);
Mu(:,SuperCounter)=Carry(1:K,1)/N;
fprintf(repmat('\b',1,LineLength))
SuperCounter=SuperCounter+1;
end

MeanProb(:,NCounter)=sum(Prob,2)/Number_of_Runs;
VarProb(:,NCounter)=sum((Prob-MeanProb(:,NCounter)).^2,2)/(Number_of_Runs-1);
MeanMu(:,NCounter)=sum(Mu,2)/Number_of_Runs;
VarMu(:,NCounter)=sum((Mu-MeanMu(:,NCounter)).^2,2)/(Number_of_Runs-1);
NCounter=NCounter+1;
end
disp(MeanMu)
disp(VarMu)
disp(MeanProb)
disp(VarProb)


figure
Counter=1;
while Counter<=K
    scatter(List,VarMu(Counter,:)/VarMu(Counter,1),'k')
    hold on
    Counter=Counter+1;
end
x=linspace(List(1,1)-5,List(1,end)+10);
y=List(1,1)./x;
plot(x,y,'LineWidth',2)
xlim([List(1,1)-5,List(1,end)+10])
xlabel('N')
ylabel('\mu Sample Varience (normalised)')
print('Image_2','-depsc')

figure
Counter=1;
while Counter<=K
    scatter(List,VarProb(Counter,:)/VarProb(Counter,1),'k')
    hold on
    Counter=Counter+1;
end
x=linspace(List(1,1)-5,List(1,end)+10);
y=List(1,1)./x;
plot(x,y,'LineWidth',2)
xlim([List(1,1)-5,List(1,end)+10])
xlabel('N')
ylabel('Prob Sample Varience (normalised)')
print('Image_3','-depsc')


function Simulations = Gibbs(Data,sigma_0,alpha_0,beta_0,mu_0,tau_0,N)
K=size(Data,1)-1;
T=size(Data,2)-1;

%Standardizes parameters 
sigma_0=sigma_0^2;
tau_0=tau_0^2;

%Creates inital state and memory for simulations
InitialState=zeros(2*K+1,1);
InitialState(1:K)=mu_0;
InitialState(K+1:2*K)=beta_0/alpha_0;
InitialState(2*K+1)=60;

Simulations=zeros(2*K+1,N+1);
Simulations(:,1)=InitialState;
Counter=2;

while Counter<=N+1
    SubCounter=1;
    while SubCounter<=2*K+1
        if SubCounter<=K
            %   SIMULATING MU PARAMETERS
            sigma_k=Simulations(SubCounter+K,Counter-1);
            sumData=sum([Data{SubCounter+1,2:T+1}]);
            theta=Simulations(2*K+1,Counter-1);
            Mean=(sigma_k^(-1)*sumData+theta*sigma_0^(-1))/(T*sigma_k^(-1)+sigma_0^(-1));
            Var=1/(T*sigma_k^(-1)+sigma_0^(-1));
            Simulations(SubCounter,Counter)=normrnd(Mean,sqrt(Var));
        elseif SubCounter<=2*K
            %   SIMULATING SIGMA PARAMETERS
            sumData=sum(([Data{SubCounter-K+1,2:T+1}]-Simulations(SubCounter-K,Counter)).^2);
            Mean=alpha_0+T/2;
            Var=beta_0+sumData/2;
            Simulations(SubCounter,Counter)=(gamrnd(Mean,1/Var))^-1;
        else
            %   SIMULATING THETA PARAMETER
            sumData=sum(Simulations(1:K,Counter));
            Mean=(sigma_0^-1*sumData+mu_0*tau_0^-1)/(K*sigma_0^-1+tau_0^-1);
            Var=1/(K*sigma_0^-1+tau_0^-1);
            Simulations(SubCounter,Counter)=normrnd(Mean,sqrt(Var));
        end
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end
end