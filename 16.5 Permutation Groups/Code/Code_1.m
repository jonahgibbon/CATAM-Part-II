A=[1 2 3 4;3 4 2 1]
B=[1 2 3 4;2 1 3 4]

Multiply(A,B)
Inverse(A)

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
    answer(1,:)=A(1,:);
    Counter=1;
    while Counter<=size(A,2)
        answer(2,A(2,Counter))=Counter;
        Counter=Counter+1;
    end
end