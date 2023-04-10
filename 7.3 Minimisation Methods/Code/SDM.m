digits(16)
Point=[-1;-1.3];
Min=[2^(-2/3)-1;-2^(-1/3)];
LambdaSigFig=5;
Iterates=10;

tic
Iterations=zeros(1,10);
Iterations(1,1)=0;
Iterations(1,2)=Point(1);
Iterations(1,3)=Point(2);
Iterations(1,5)=vpa(func(Point));
Iterations(1,6)=norm(gradient(Point));

Counter=1;
n=1;
while n<Iterates+1
s=-gradient(Point);
Counter=Counter+1;
Iterations=[Iterations;zeros(1,10)];
Iterations(n+1,1)=n;

% %Plot of $phi(lambda)$ with prompt
Inc=10000000;
Dif=vpa(func(Point+Inc*s)-func(Point));
while Dif>0
        Inc=Inc/10;
        Dif=vpa(func(Point+Inc*s)-func(Point));
end
XPhi=linspace(0,5*Inc);
YPhi=zeros(1,100);
SubCounter=1;
while SubCounter<=100
    YPhi(SubCounter)=func(Point+XPhi(SubCounter)*s);
    SubCounter=SubCounter+1;
end
figure
plot(XPhi,YPhi)
prompt=input(sprintf('Enter %d-th lambda value: ',n));
close

% prompt=[];
if isempty(prompt)==1
    [Point,Iterations(Counter-1,4)]=lambda(Point,s,LambdaSigFig);
else
    Point=vpa(Point+prompt*s);
    Iterations(Counter-1,4)=prompt;
end


Iterations(Counter,2)=Point(1);
Iterations(Counter,3)=Point(2);
Iterations(Counter,5)=vpa(func(Point));
Iterations(Counter,6)=norm(gradient(Point));
Iterations(Counter,7)=Iterations(Counter,5)-Iterations(Counter-1,5);
Iterations(Counter,8)=Iterations(Counter,7)/Iterations(Counter-1,7);
Iterations(Counter,9)=norm(Point-[Iterations(Counter-1,2);Iterations(Counter-1,3)]);
if Counter>=3
    Iterations(Counter,10)=Iterations(Counter,9)/Iterations(Counter-2,9);
end

%Termination Process
% if norm(s)<=10^-7
%     disp('Termination - Gradient small')
%     break
% end
% if (Iterations(n+1,2)-Iterations(n,2))^2+...
%         (Iterations(n+1,3)-Iterations(n,3))^2<=10^-20
%     disp('Termination - Iterations close together')
%     break
% end
n=n+1;
end
disp(Iterations)
latex(sym(vpa(Iterations)))

% disp(norm(Point-[1;1]))

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

toc

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
answer=x(1)+x(2)+x(1)^2/4-x(2)^2+(x(2)^2-x(1)/2)^2;
end

function answer = gradient(x)
answer=[1+x(1)-x(2)^2;1-4*x(2)*(-x(2)^2+x(1)/2)-2*x(2)];
end

% function answer = func(x)
% answer=(1-x(1))^2+80*(x(2)-x(1)^2)^2;
% end
% 
% function answer = gradient(x)
% answer=[320*x(1)^3+(2-320*x(2))*x(1)-2;160*(x(2)-x(1)^2)];
% end