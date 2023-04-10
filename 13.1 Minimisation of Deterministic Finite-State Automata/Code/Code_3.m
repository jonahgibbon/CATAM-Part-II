n=4;
k=2;
Everything=GenerateEverything(n,k);
disp(size(Everything,1))

function TransitionTables = GenerateEverything(n,k)
States=(1:n)';
Fill=cell(1,k);
Counter=1;
while Counter<=k
    Fill{1,Counter}=States;
    Counter=Counter+1;
end
Filler=Fill;
[Filler{:}] = ndgrid(Fill{:});
Bank=cell2mat(cellfun(@(m)m(:),Filler,'uni',0));
BankSize=size(Bank,1);

States=(1:BankSize)';
Fill=cell(1,n);
Counter=1;
while Counter<=n
    Fill{1,Counter}=States;
    Counter=Counter+1;
end
Filler=Fill;
[Filler{:}] = ndgrid(Fill{:});
TableIndex=cell2mat(cellfun(@(m)m(:),Filler,'uni',0));
Counter=1;
TransitionTables=cell(0,1);
Tally=0;
while Counter<=size(TableIndex,1)
    Table=zeros(n,k);
    SubCounter=1;
    while SubCounter<=n
        Table(SubCounter,:)=Bank(TableIndex(Counter,SubCounter),:);
        SubCounter=SubCounter+1;
    end
    M=Accessible(Table,n,k);
    if sum(M(:,2))==n
        Tally=Tally+1;
        TransitionTables{Tally,1}=Table;
    end
    Counter=Counter+1;
end
end

function AccessibleStates = Accessible(TransTab,n,k)
Queue=zeros(1,n);
Queue(1,1)=1;
LengthQueue=1;
Counter=1;
AccessibleStates=zeros(n,2);
AccessibleStates(:,1)=1:n;
AccessibleStates(1,2)=1;

while Counter<=LengthQueue
    SubCounter=1;
    while SubCounter<=k
        Element=TransTab(Queue(Counter),SubCounter);
        if AccessibleStates(Element,2)==0
            AccessibleStates(Element,2)=1;
            Queue(LengthQueue+1)=Element;
            LengthQueue=LengthQueue+1;
        end
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end
end