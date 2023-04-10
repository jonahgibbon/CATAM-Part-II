Count=1;
digits(10)
Data=zeros(4,6);
while Count<=6
    SubCount=1;
    while SubCount<=4
        Data(SubCount,Count)=2^Count * Count^(Count*SubCount);
        SubCount=SubCount+1;
    end
    Count=Count+1;
end
disp(Data')
latex(sym(vpa(Data')))
