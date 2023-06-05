xlim = mean(pic.xi)+0.2*[-1 1];
twpe = 0;
A = mean(pic.twpelim(twpe).xlim(xlim).A,1);
By = mean(pic.twpelim(twpe).xlim(xlim).By,1);
Bx = mean(pic.twpelim(twpe).xlim(xlim).Bx,1);



maxBy = max(By);

[I,J] = find(By>maxBy*0.5);


%plot(A,By,A(J),By(J),'.',[min(A) max(A)],0.5*maxBy*[1 1],'-')

plotyy(pic.zi,[Bx;By],pic.zi,A)