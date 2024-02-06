
function isRadial=checkRadial(nodeList, gStruct)

stDist=gStruct.coordinates(nodeList(2:end), :)-gStruct.coordinates(nodeList(1), :);
stDist=sum(stDist.^2, 2);
isRadialSt=all(stDist(2:end)>stDist(1:end-1));
if isRadialSt
    rvDist=gStruct.coordinates(nodeList(1:end-1), :)-gStruct.coordinates(nodeList(end), :);
    rvDist=sum(rvDist.^2, 2);
    isRadialRv=all(rvDist(1:end-1)>rvDist(2:end));
    isRadial=isRadialRv;
else
    isRadial=false;
end
end
