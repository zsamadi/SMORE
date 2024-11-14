function [Rkl, Rbh, Rcld]=getKLCorr(X)
X=X./sum(X);

X(X==0)=1e-15;

numCols=size(X, 2);
Rkl=zeros(numCols);
Rbh=zeros(numCols);


for ic =1:size(X, 2)
    for jc=1:size(X, 2)
        Rkl(ic, jc)=sum(X(:, ic).*(log2(X(:, ic))-log2(X(:, jc))));
        Rbh(ic, jc)=-log(sum(sqrt(X(:, ic).*X(:, jc))));

    end
end

Rcld=eye(numCols);

for ic =1:size(X, 2)
    icd=1:size(X, 2);
    icd(ic)=[];
    XD=X(:, ic)-X(:, icd);
    [~,XDI]=min(abs(XD),[], 2);
    for icdi =1:length(icd)
        XDI1=sum(XDI==icdi)/length(XDI);
        Rcld(ic, icd(icdi))=XDI1;
    end
end

Rcld=(Rcld+Rcld.')/2;





