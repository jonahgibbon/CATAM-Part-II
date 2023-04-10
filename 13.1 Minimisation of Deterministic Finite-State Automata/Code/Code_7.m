n=4;
k=2;

tic
Everything=GenerateEverything(n,k);
States=(0:1)';
Fill=cell(1,n);
Counter=1;
while Counter<=n
    Fill{1,Counter}=States;
    Counter=Counter+1;
end
Filler=Fill;
[Filler{:}] = ndgrid(Fill{:});
Bank=cell2mat(cellfun(@(m)m(:),Filler,'uni',0));

Fill=cell(1,2);
Fill{1,1}=(1:size(Everything,1))';
Fill{1,2}=(1:size(Bank,1))';
Filler=Fill;
[Filler{:}] = ndgrid(Fill{:});
TableIndex=cell2mat(cellfun(@(m)m(:),Filler,'uni',0));


MinimalDFAs=cell(0,1);
Tally=0;
Counter=1;
while Counter<=size(TableIndex,1)
    Table=[Everything{TableIndex(Counter,1),1} Bank(TableIndex(Counter,2),:)'];
    [Table,TabSize]=Minimise(Table,n,k);
    if TabSize==n
        Tally=Tally+1;
        MinimalDFAs{Tally,1}=Table;
    end
    Counter=Counter+1;
end


fData=zeros(2^n,2);
fData(:,1)=0:2^n-1;
SuperCounter=1;
while SuperCounter<=size(MinimalDFAs,1)
    lineLength=fprintf('Completion: %.1f%%',SuperCounter*100/size(MinimalDFAs,1));
    TransTab=MinimalDFAs{SuperCounter,1};
    Length=n;
    Flag=0;
    while Length<=2*n-1
        Words=GenerateWords(Length,k);
        Counter=1;
        while Counter<=k^Length
            if runWord(TransTab,Words(Counter,:),k)==1
                Flag=1;
                break
            end
            Counter=Counter+1;
        end
        Length=Length+1;
        if Flag==1
            break
        end
    end
    fprintf(repmat('\b',1,lineLength))
    if Flag==0
        if TransTab(1,k+1)==1
            LanguageSize=1;
        else
            LanguageSize=0;
        end
        Length=1;
        while Length<=n-1
            Words=GenerateWords(Length,k);
            Counter=1;
            while Counter<=k^Length
                if runWord(TransTab,Words(Counter,:),k)==1
                    LanguageSize=LanguageSize+1;
                end
                Counter=Counter+1;
            end
            Length=Length+1;
        end
        if LanguageSize==1
            disp('LanguageSize=1')
            disp(MinimalDFAs{SuperCounter,1})
            disp(latex(sym(MinimalDFAs{SuperCounter,1})))
            disp('*****')
        end
        if LanguageSize==2^(n-1)-1
            disp('LanguageSize=2^(n-1)-1')
            disp(MinimalDFAs{SuperCounter,1})
            disp(latex(sym(MinimalDFAs{SuperCounter,1})))
            disp('*****')
        end
        fData(LanguageSize+1,2)=fData(LanguageSize+1,2)+1;
    end
    SuperCounter=SuperCounter+1;
end
fData(:,2)=fData(:,2)./factorial(n-1);
disp(fData)
latex(sym(fData'))
toc




function answer = GenerateWords(Length,k)
States=(1:k)';
Fill=cell(1,Length);
Counter=1;
while Counter<=Length
    Fill{1,Counter}=States;
    Counter=Counter+1;
end
Filler=Fill;
[Filler{:}] = ndgrid(Fill{:});
answer=cell2mat(cellfun(@(m)m(:),Filler,'uni',0));
Counter=1;
end


function verdict = runWord(TransTab,Word,k)
    Length=size(Word,2);
    Counter=1;
    State=1;
    while Counter<=Length
        State=TransTab(State,Word(1,Counter));
        Counter=Counter+1;
    end
    verdict=TransTab(State,k+1);
end

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
Counter=1;
while Counter<=n-1
    SubCounter=Counter+1;
    while SubCounter<=n
        if TransTab(SubCounter,3)~=TransTab(Counter,3)
            HopTable(SubCounter,Counter)=1;
        end
        SubCounter=SubCounter+1;
    end
    Counter=Counter+1;
end

flag=1;
while flag==1
    flag=0;
    Counter=1;
    while Counter<=n-1
        SubCounter=Counter+1;
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
                        if HopTable(YState,XState)==1
                            HopTable(SubCounter,Counter)=1;
                            flag=1;
                        else
                        end
                    end
                    Sub2Counter=Sub2Counter+1;
                end
            end
            SubCounter=SubCounter+1;
        end
        Counter=Counter+1;
    end
end
end

function TransTab = Clean(TransTab,n,k)
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