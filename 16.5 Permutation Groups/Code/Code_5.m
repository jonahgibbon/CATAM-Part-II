n=5;
Generators={[1 2 3 4 5;1 3 4 2 5],[1 2 3 4 5;5 2 3 4 1]};
Alpha=2;

Orb=Orbit(n,Generators,Alpha);
disp(Orb)

function [answer,siz] = Orbit(n,Generators,Alpha)
answer=cell(n,2);
answer{Alpha,1}=Alpha;
answer{Alpha,2}=[1:n;1:n];
Queue=[Alpha];

Counter=1;
UpperLimit=1;
while Counter<=UpperLimit
    GenCounter=1;
    while GenCounter<=size(Generators,2)
        Image=Generators{1,GenCounter}(2,Queue(Counter));
        if isempty(answer{Image,1})==1
            UpperLimit=UpperLimit+1;
            Queue=[Queue,Image];
            answer{Image,1}=Image;
            answer{Image,2}=Multiply(Generators{1,GenCounter},answer{Queue(Counter),2});
        end
        GenCounter=GenCounter+1;
    end
    Counter=Counter+1;
end
siz=size(Queue,2);
end

function answer = Multiply(A,B)
    Counter=1;
    answer=zeros(2,size(A,2));
    while Counter<=size(A,2)
        answer(1,Counter)=Counter;
        answer(2,Counter)=A(2,B(2,Counter));
        Counter=Counter+1;
    end
end