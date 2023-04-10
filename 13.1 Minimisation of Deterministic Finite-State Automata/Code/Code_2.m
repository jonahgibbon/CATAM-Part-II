TransTab=[
2 5 1
10 10 1
1 2 0
8 3 0
9 2 0
9 2 0
1 9 1
4 6 0
3 6 0
9 2 1];
n=size(TransTab,1);
k=size(TransTab,2)-1;

[m,n]=Clean(TransTab,n,k)

function [TransTab,n] = Clean(TransTab,n,k)
AccessibleStates=Accessible(TransTab,n,k);
ChopCount=0;
Counter=1;
while Counter<=n
    if AccessibleStates(Counter,2)==0
        TransTab(Counter-ChopCount,:)=[];
        SubCounter=1;
        while SubCounter<=n-ChopCount-1
            Sub2Counter=1;
            while Sub2Counter<=k
                if TransTab(SubCounter,Sub2Counter)>Counter-ChopCount
                    TransTab(SubCounter,Sub2Counter)=TransTab(SubCounter,Sub2Counter)-1;
                end
                Sub2Counter=Sub2Counter+1;
            end
            SubCounter=SubCounter+1;
        end
        ChopCount=ChopCount+1;
    end
    Counter=Counter+1;
end
n=n-ChopCount;
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