figure
x=linspace(-1.5,1.5,1000);
y=linspace(-1.5,1.5,1000);
[X,Y]=meshgrid(x,y);
Z=func(X,Y);
contour(X,Y,Z,-1:0.4:7)
colorbar
hold on
scatter(2^(-2/3)-1,-2^(-1/3),'+','Red','LineWidth',1)
legend('','Minimum')
print('Image','-depsc')

function answer = func(x,y)
answer=x+y+x.^2./4-y.^2+(y.^2-x./2).^2;
end

% function answer = func(x,y)
% answer=(1-x).^2+80*(y-x.^2).^2;
% end
