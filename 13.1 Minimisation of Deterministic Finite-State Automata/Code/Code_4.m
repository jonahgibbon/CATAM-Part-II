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
9 2 1
    ];
n=size(TransTab,1);
k=size(TransTab,2)-1;

tic
[TransTab,n]=Clean(TransTab,n,k);
[TransTab,n]=Minimise(TransTab,n,k)
toc

function [TransTab,n] = Minimise(TransTab,n,k)
HopTable=Hopcroft(TransTab,n,k);
SubCounter=n;
while SubCounter>1
    Counter=1;
    while Counter<=SubCounter-1
        if HopTable(SubCounter,Counter)==0
            Sub2Counter=1;
            while Sub2Counter<=k
                 Sub3Counter=1;
                 while Sub3Counter<=n
                     if TransTab(Sub3Counter,Sub2Counter)==SubCounter
                         TransTab(Sub3Counter,Sub2Counter)=Counter;
                     end
                     Sub3Counter=Sub3Counter+1;
                 end
                 Sub2Counter=Sub2Counter+1;
            end
        end
        Counter=Counter+1;
    end
    SubCounter=SubCounter-1;
end
TransTab=Clean(TransTab,n,k);
n=size(TransTab,1);
end

function HopTable = Hopcroft(TransTab,n,k)
HopTable=zeros(n,n);
Queue=zeros(0,2);
Counter=1;
while Counter<=n-1
    SubCounter=Counter+1;
    while SubCounter<=n
        if TransTab(SubCounter,k+1)~=TransTab(Counter,k+1)
            HopTable(SubCounter,Counter)=1;
            Queue=[Queue; SubCounter,Counter];
        end
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end
Lists=cell(n,n);
Counter=1;
while Counter<=n-1
    SubCounter=1+Counter;
    while SubCounter<=n
        if HopTable(SubCounter,Counter)==0
            Sub2Counter=1;
            while Sub2Counter<=k
                XState=TransTab(SubCounter,Sub2Counter);
                YState=TransTab(Counter,Sub2Counter);
                if XState==YState
                else
                    if XState<YState
                    else
                        Z=YState;
                        YState=XState;
                        XState=Z;
                    end
                    Lists{YState,XState}=[Lists{YState,XState}; SubCounter,Counter];
                end
                Sub2Counter=Sub2Counter+1;
            end
        end
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end
QueueLength=size(Queue,1);
Counter=1;
while Counter<=QueueLength
    CurrentList=Lists{Queue(Counter,1),Queue(Counter,2)};
    SubCounter=1;
    while SubCounter<=size(CurrentList,1)
        if HopTable(CurrentList(SubCounter,1),CurrentList(SubCounter,2))==0
            HopTable(CurrentList(SubCounter,1),CurrentList(SubCounter,2))=1;
            Queue=[Queue;CurrentList(SubCounter,:)];
            QueueLength=QueueLength+1;
        end
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end
end

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