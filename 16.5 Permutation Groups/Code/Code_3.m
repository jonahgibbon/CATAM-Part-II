[Array,NewGenerators]=SAOS(5,{[1 2 3 4 5;2 1 3 4 5],[1 2 3 4 5;2 3 1 4 5],...
    [1 2 3 4 5;1 3 5 4 2],[1 2 3 4 5;5 2 1 4 3],[1 2 3 4 5;1 2 3 4 5]})


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