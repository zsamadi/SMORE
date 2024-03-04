function [u, iu, ju]=uniqueSorted(sortedArray, maxElem)

w=size(sortedArray, 2);
if w>1
    sortedVector=sum(sortedArray.*maxElem.^(w-1:-1:0), 2);
else
    sortedVector=sortedArray;
end

sortedVector2=sortedVector(2:end);
sortedVectorDiff=sortedVector2-sortedVector(1:end-1);
ind2=(2:length(sortedVector)).';
iu=[1;ind2(sortedVectorDiff~=0)];
u=sortedArray(iu, :);

if nargout>2
    iU2R=[iu(2:end)-iu(1:end-1); length(sortedVector)-iu(end)+1];
    ju=repelem((1:length(iU2R)).', iU2R,1);
end