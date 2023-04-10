tic
StartN=1;
FinishN=10;
SampleSize=1000;

ExactPN=[1 1;2 3/4;3 1/2;4 3/8;5 19/40;6 53/120];
Tally=zeros(FinishN-StartN+1,2);
n=StartN;
while n<=FinishN
Counter=1;
SubTally=0;
while Counter<=SampleSize
    lineLength = fprintf('n = %d, %.1f%% complete',n,Counter*100/SampleSize);
    Perms=random(2,n);
    Data=Order(n,Perms);
    if Data(1,1)==factorial(n)
        SubTally=SubTally+1;
    end
    Counter=Counter+1;
    fprintf(repmat('\b',1,lineLength))
end
Tally(n-StartN+1,1)=n;
Tally(n-StartN+1,2)=SubTally/SampleSize;
n=n+1;
end
disp(Tally)

figure
plot(linspace(StartN,FinishN,FinishN-StartN+1),Tally(:,2),'-o','LineWidth',1)
hold on
scatter(ExactPN(:,1),ExactPN(:,2),'+','LineWidth',1)
plot(linspace(StartN,FinishN,FinishN-StartN+1),linspace(0.75,0.75,FinishN-StartN+1),'LineWidth',1)
legend('Estimated P_{n}','Exact P_{n}','k=0.75')
xlabel('n')
ylabel('Probability')
xlim([StartN,FinishN])
ylim([0,1])
print('Image_1','-depsc')
toc

function answer = random(quantity,n)
answer=cell(1,0);
Counter=1;
while Counter<=quantity
    Perm=[1:n;1:n];
    SubCounter=1;
    while SubCounter<=n-1
        Carry=Perm(2,SubCounter);
        RandIndex=randi(n-SubCounter+1)+SubCounter-1;
        Perm(2,SubCounter)=Perm(2,RandIndex);
        Perm(2,RandIndex)=Carry;
        SubCounter=SubCounter+1;
    end
    answer={answer{:,:} Perm};
    Counter=Counter+1;
end
end

function Data=Order(n,Generators)
Data=zeros(1,5);
Data(1,2)=size(Generators,2);
GenArray=Generators;
[~,GenArray]=SAOS(n,GenArray);
Data(1,3)=size(GenArray,2);

while size(GenArray,2)~=0
    Counter=1;
    while Counter<=n
        [~,Size]=Orbit(n,GenArray,Counter);
        if Size>1
            break
        end
        Counter=Counter+1;
    end
    if Counter==n+1
        break
    end
    Data(end,4)=Counter;
    Data(end,5)=Size;
    Data=[Data;0,0,0,0,0];
    GenArray=Stabiliser(n,GenArray,Counter);
    Data(end,2)=size(GenArray,2);
    [~,GenArray]=SAOS(n,GenArray);
    Data(end,3)=size(GenArray,2);
end
Data(end,1)=1;
Counter=size(Data,1)-1;
while Counter>=1
    Data(Counter,1)=Data(Counter+1,1)*Data(Counter,5);
    Counter=Counter-1;
end
end

function S = Stabiliser(n,Generators,Alpha)
OrbitArray=Orbit(n,Generators,Alpha);
S=cell(1,0);

Counter=1;
while Counter<=n
    SubCounter=1;
    while SubCounter<=size(Generators,2)
        if size(OrbitArray{Counter,1},1)==0
            break
        end
        %Calculates yt
        Perm=Multiply(Generators{1,SubCounter},OrbitArray{Counter,2});
        %Equals Varphi(yt)
        Varphi=OrbitArray{Perm(2,Alpha),2};
        S={S{1,:} Multiply(Inverse(Varphi),Perm)};
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end
end

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

function [Array,NewGenerators] = SAOS(n,OldGenerators)
Identity=[1:n;1:n];
NewGenerators=cell(1,0);

Array=cell(n,n);
Counter=1;
while Counter<=size(OldGenerators,2)
    RowCounter=1;
    %Checks if current permutation is the identity
    while RowCounter<=n     
        if OldGenerators{Counter}==Identity
            break
        end
        %Checks if the permutation fixes Counter
        if OldGenerators{Counter}(2,RowCounter)~=RowCounter
            %If cell is empty, fills cell in array
            if isempty(Array{RowCounter,OldGenerators{Counter}...
                (2,RowCounter)})==1
                Array{RowCounter,OldGenerators{Counter}(2,RowCounter)}=...
                        OldGenerators{Counter};
                NewGenerators={NewGenerators{1,:} OldGenerators{Counter}};
                    break
            %If cell isn't empty, modifys permutation to fix Counter
            else
                OldGenerators{Counter}=Multiply(Inverse(Array{...
                    RowCounter,OldGenerators{Counter}(2,RowCounter)}),...
                    OldGenerators{Counter});
            end
        end
        RowCounter=RowCounter+1;
    end
    Counter=Counter+1;
end
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

function answer = Inverse(A)
    answer=zeros(2,size(A,2));
    answer(1,:)=A(2,:);
    answer(2,:)=A(1,:);
    Counter=1;
    while Counter<=size(A,2)
        if answer(1,Counter)==Counter
            Counter=Counter+1;
        else
            CarryVector=answer(:,answer(1,Counter));
            answer(:,answer(1,Counter))=answer(:,Counter);
            answer(:,Counter)=CarryVector;
        end
    end
end