function [whichPH, order4ES]=findPH(modalities)


nModalities=numel(modalities);
[~,alphabeticalOrder]=sort(modalities);

loop=1;
for i=1:nModalities
    for j=1:nModalities
        if i~=j
            orderPH(loop,1)=alphabeticalOrder(i);
            orderPH(loop,2)=alphabeticalOrder(j);
            loop=loop+1;
        end
    end
end

loop=1;
for i=1:nModalities
    for j=1:nModalities
        if i<j
            order4ES(loop,1)=i;
            order4ES(loop,2)=j;
            loop=loop+1;
        end
    end
end

[order4ES(:,1),b]=sort(order4ES(:,1),'descend');
order4ES(:,2)=order4ES(b,2);

for i=1:size(order4ES,1)
    for j=1:size(orderPH)
        if isequal(order4ES(i,:),orderPH(j,:))
            whichPH(i,1)=j;
        end
    end
end




end

