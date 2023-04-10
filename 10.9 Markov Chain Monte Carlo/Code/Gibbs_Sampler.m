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

StartVal=zeros(2*K+1,1);
StartVal(1:2*K+1)=60;
StartVal(K+1:2*K)=100;

N=5;
tic
Sim=Gibbs(Data,sigma_0,alpha_0,beta_0,mu_0,tau_0,N,StartVal)
latex(sym(vpa(Sim)))
toc

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
    lineLength = fprintf('%.1f%% complete',Counter*100/N);
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
    fprintf(repmat('\b',1,lineLength))
    Counter=Counter+1;
end
end