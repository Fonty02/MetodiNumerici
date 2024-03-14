function [x1,n_it] = newton(atol,rtol,max_it,x0,J,f)
n_it=0;
errore=true;
while errore
    dk=feval(J,x0)\feval(f,x0); %invece di -f
    x1=x0-dk; %invece di +f
    n_it=n_it+1;
    errore=norm(abs(x1-x0)./(atol+rtol*abs(x1)),'inf')>1 && n_it<max_it;
    x0=x1;
end

