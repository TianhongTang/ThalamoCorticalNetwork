xdata=[0,1,2,3,4,5,6,7];
ydata=[100,7,5,4,3.5,3.3,3.5,1.7];
plot(xdata,ydata)
[x,resnorm]=lsqcurvefit('myfun',[1,2,10],xdata,ydata)
