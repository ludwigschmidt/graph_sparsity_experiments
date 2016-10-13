function y = Afun_spg(x,mode,A,At)

if mode==1
    y = A(x);
elseif mode==2
    y = At(x);
end