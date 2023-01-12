function addPvalue(locX,locY,pV,pCritical, color)

if pV<=pCritical(1)

    if pV>=0.001
        text(locX, locY, append('p = ', sprintf('%.3f', pV)), HorizontalAlignment="center", Color=color, FontWeight="bold");
    else
        text(locX, locY, append('p < ', sprintf('%.3f', 0.001)), HorizontalAlignment="center", Color=color, FontWeight="bold");
    end

elseif pV<=pCritical(2) & pV>pCritical(1) 

    text(locX, locY, append('p = ', sprintf('%.3f', pV)), HorizontalAlignment="center", Color=color, FontWeight="normal");

else
end

end