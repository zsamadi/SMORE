function nkl=nchoosekLog(n, k)

nkl=zeros(length(n),1);

for in=1:length(n)
    ki=k(in);
    ni=n(in);

    if ki<=ni

        nfact=sum(log((1:ni)));
        k0=max(ki,1);
        kfact=sum(log((1:k0)));
        nmk0=max(ni-ki, 1);
        nmk=sum(log((1:nmk0)));

        nkli=nfact-kfact-nmk;
    else
        nkli=-30;
    end
    nkl(in)=nkli;
end


