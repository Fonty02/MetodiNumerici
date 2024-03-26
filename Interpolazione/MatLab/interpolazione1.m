function a=interpolazione1(x,f)
    diff_div=f;
    a(1)=diff_div(1);
    for i=2:4
        diff_div=(diff_div(2:end)-diff_div(1:end-1))./(x(i:end)-x(1:end-(i-1)));  %in C/Java avremmo fattto un altro for interno per calcolare elemento per elemento
        a(i)=diff_div(1);
    end
end