digits(16)
NewPoint=[-1;-1.3];
Min=[2^(-2/3)-1,-2^(-1/3)];
N=2;
LambdaSigFig=5;

Iterates=10;

tic
Iterations=zeros(Iterates+1,N+4);
Iterations(1,1)=0;
Counter=1;
while Counter<=N
    Iterations(1,Counter+1)=NewPoint(Counter);
    Counter=Counter+1;
end

Iterations(1,N+3)=vpa(func(NewPoint));
Iterations(1,N+4)=norm(grad(NewPoint));

Differences=zeros(floor(Iterates/N),7);
Index=1;
LowerPoint=NewPoint;

n=0;
while n<Iterates
    if mod(n,N)==0
        s=-grad(NewPoint);
    else
        s=-grad(NewPoint)+norm(grad(NewPoint))^2/norm(grad(OldPoint))^2*s;
    end
    OldPoint=NewPoint;
    

%     %Plot of $phi(lambda)$ with prompt
    Inc=10000000;
    Dif=vpa(func(NewPoint+Inc*s)-func(NewPoint));
    while Dif>0
            Inc=Inc/10;
            Dif=vpa(func(NewPoint+Inc*s)-func(NewPoint));
    end
    XPhi=linspace(0,5*Inc);
    YPhi=zeros(1,100);
    SubCounter=1;
    while SubCounter<=100
        YPhi(SubCounter)=func(NewPoint+XPhi(SubCounter)*s);
        SubCounter=SubCounter+1;
    end
    figure
    plot(XPhi,YPhi)
    prompt=input(sprintf('Enter %d-th lambda value: ',n+1));
    close

%     prompt=[];
    if isempty(prompt)==1
        [NewPoint,Iterations(n+1,N+2)]=lambda(NewPoint,s,LambdaSigFig);
    else
        NewPoint=vpa(NewPoint+prompt*s);
        Iterations(n+1,N+2)=prompt;
    end

    Iterations(n+2,1)=n+1;
    Counter=1;
    while Counter<=N
        Iterations(n+2,Counter+1)=NewPoint(Counter);
        Counter=Counter+1;
    end
    Iterations(n+2,N+3)=vpa(func(NewPoint));
    Iterations(n+2,N+4)=vpa(norm(grad(NewPoint)));

    if and(mod(n+1,N)==0,n~=0)
        Differences(Index,1)=Index-1;
        Differences(Index,2)=vpa(func(NewPoint))-vpa(func(LowerPoint));
        Differences(Index,5)=norm(vpa(NewPoint)-vpa(LowerPoint));
        if Index~=1
            Differences(Index,3)=Differences(Index,2)/Differences(Index-1,2);
            Differences(Index,4)=Differences(Index,2)/Differences(Index-1,2)^2;
            Differences(Index,6)=Differences(Index,5)/Differences(Index-1,5);
            Differences(Index,7)=Differences(Index,5)/Differences(Index-1,5)^2;
        end
        LowerPoint=vpa(NewPoint);
        Index=Index+1;
    end

% TERMINATION PROCESS
%     if Iterations(n+2,N+4)<=10^-7
%         disp('Termination - Gradient small')
%         break
%     end
%     if norm(NewPoint-OldPoint)<=10^-10
%         disp('Termination - Iterations close together')
%         break
%     end
    n=n+1;
end
toc

disp(Iterations)
latex(sym(vpa(Iterations)))

% disp(Differences)
% latex(sym(vpa(Differences)))

% PLOT ITERATES
figure
x=linspace(-1.5,0);
y=linspace(-1.5,0);
SupCounter=1;
Z=zeros(100,100);
while SupCounter<=100
    SubCounter=1;
    while SubCounter<=100
        Z(SubCounter,SupCounter)=func([x(SupCounter);y(SubCounter)]);
        SubCounter=SubCounter+1;
    end
    SupCounter=SupCounter+1;
end
contour(x,y,Z)
colorbar
hold on
scatter(Min(1),Min(2),'+','Red','LineWidth',1)
plot(Iterations(:,2),Iterations(:,3),'k')
scatter(Iterations(:,2),Iterations(:,3),'k')
print('Image','-depsc')


%ERROR ANAYLISIS FOR FUNCTION (5)
% figure
% x=zeros(floor(Iterates/2),1);
% y=zeros(floor(Iterates/2),1);
% Counter=1;
% while Counter<=floor(Iterates/2)
%     x(Counter)=log(Iterations(2*Counter-1,N+3));
%     y(Counter)=log(Iterations(2*Counter+1,N+3));
%     Counter=Counter+1;
% end
% X=x(1:end-2);
% Y=y(1:end-2);
% m=(sum(X.*Y)-sum(X)*sum(Y)/(size(X,1)))/(sum(X.^2)-sum(X)^2/(size(X,1)));
% b=(sum(Y)-m*sum(X))/size(X,1);
% z=m*x+b;
% Z=m*X+b;
% plot(x,z)
% hold on
% scatter(x,y)
% xlabel('log E_{kN}')
% ylabel('log E_{(k+1)N}')
% legend('Least Squares Regression Line','','Location','northwest')
% print('Regression Line','-depsc')
% R=1-sum((Y-Z).^2)/(sum((Y-(sum(Y)/size(Y,1))).^2));

function [Approx,L] = lambda(x,s,d)
    Inc=10000000;
    Dif=vpa(func(x+Inc*s)-func(x));
    while Dif>0
            Inc=Inc/10;
            Dif=vpa(func(x+Inc*s)-func(x));
    end
    Counter=1;
    Approx=x;
    L=0;
    while Counter<=d
        Inc=Inc/10;
        SubCounter=1;
        Dif=vpa(func(Approx+SubCounter*Inc*s)-func(Approx+(SubCounter-1)*Inc*s));
        while Dif<0
            SubCounter=SubCounter+1;
            Dif=vpa(func(Approx+SubCounter*Inc*s)-func(Approx+(SubCounter-1)*Inc*s));
        end
        L=L+(SubCounter-2)*Inc;
        Approx=vpa(Approx+(SubCounter-2)*Inc*s);
        Counter=Counter+1;
    end
end

%FUNCTIONS TO INTERCHANGE

function answer = func(x)
answer=vpa(x(1)+x(2)+x(1)^2/4-x(2)^2+(x(2)^2-x(1)/2)^2);
end

function answer = grad(x)
answer=vpa([1+x(1)-x(2)^2;1-4.*x(2)*(-x(2)^2+x(1)/2)-2.*x(2)]);
end

% function answer = func(x)
% answer=(1-x(1))^2+80*(x(2)-x(1)^2)^2;
% end
% 
% function answer = grad(x)
% answer=[320*x(1)^3+(2-320*x(2))*x(1)-2;160*(x(2)-x(1)^2)];
% end