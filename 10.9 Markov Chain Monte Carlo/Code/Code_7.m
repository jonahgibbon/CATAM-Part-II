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

CalculatedMean=[81.4968825139968;52.1355802069355;50.672260968429;52.5080817745987;...
    48.1389564329937;69.784992226184;51.6374396731938;48.7565228578938;...
    62.2933104774652;74.386133868267;49.8382175215858;59.9520363668573;...
    48.2862648696441;51.5355709368392];
CalculatedVar=[1.9394709;0.43097878;0.67722952;0.66446875;0.20170042;2.7968475;0.89288769;...
    0.25573182;0.96594047;0.74396463;0.18181906;0.93427989;0.89550653;0.2827326];

StartingValue=zeros(2*K+1,1);
StartingValue(1:2*K+1)=60;
StartingValue(K+1:2*K)=100;
StartingValue(2*K+1)=60;

List=[0,70];
N=30;

NCounter=1;
MeanProb=zeros(K,size(List,2));
MeanMu=zeros(K,size(List,2));
VarProb=zeros(K,size(List,2));
VarMu=zeros(K,size(List,2));

while NCounter<=size(List,2)
M=List(1,NCounter);
Number_of_Runs=500;

Prob=zeros(K,Number_of_Runs);
Mu=zeros(K,Number_of_Runs);

SuperCounter=1;
while SuperCounter<=Number_of_Runs

LineLength=fprintf('List: %d, Run: %.1f%% complete',M+N,SuperCounter*100/Number_of_Runs);
Sim=Gibbs(Data,sigma_0,alpha_0,beta_0,mu_0,tau_0,M+N,StartingValue);

Counter=2+M;
Carry=zeros(K,1);
while Counter<=M+N+1
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
Carry=sum(Sim(:,M+2:M+N+1),2);
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
latex(sym(vpa(VarMu)))

Norms=[List;vecnorm(MeanMu-CalculatedMean);vecnorm(VarMu-CalculatedVar)];
disp(Norms)
latex(sym(vpa(Norms)))

disp(MeanProb)
disp(VarProb)
latex(sym(vpa(VarProb)))



function Simulations = Gibbs(Data,sigma_0,alpha_0,beta_0,mu_0,tau_0,N,StartVal)
K=size(Data,1)-1;
T=size(Data,2)-1;

%Standardizes parameters 
sigma_0=sigma_0^2;
tau_0=tau_0^2;

Simulations=zeros(2*K+1,N+1);
Simulations(:,1)=StartVal;
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