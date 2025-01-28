function txt = displayCoordinates(~,info)
    x = info.Position(1);
    y = info.Position(2);
    myDatatipText = "(%s, %s)";
    txt = sprintf(myDatatipText, num2str(x), num2str(y));
end