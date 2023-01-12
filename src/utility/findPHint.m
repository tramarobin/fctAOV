function [whichPHint]=findPHint(modalities1,modalities2)

nModalities(1)=numel(modalities1);
nModalities(2)=numel(modalities2);
[~,alphabeticalOrder]=sort(modalities1);
totalPH=nModalities(1)*nModalities(2);

% cut the first effect by the second effect
loop=1;
for i=1:nModalities(1)
    for j=1:nModalities(1)
        if i~=j
            orderPH(loop,1)=alphabeticalOrder(i);
            orderPH(loop,2)=alphabeticalOrder(j);
            loop=loop+1;
        end
    end
end

loop=1;
for i=1:nModalities(1)
    for j=1:nModalities(1)
        if i<j
            order(loop,1)=i;
            order(loop,2)=j;
            loop=loop+1;
        end
    end
end

[order(:,1),b]=sort(order(:,1),'descend');
order(:,2)=order(b,2);

for i=1:size(order,1)
    for j=1:size(orderPH)
        if isequal(order(i,:),orderPH(j,:))
            whichPH(i,1)=j;
        end
    end
end

Nph=2*nModalities(2)*size(order,1);
subGroupInInt=[1:size(orderPH,1):Nph]-1;
[~,alphabeticalOrder]=sort(modalities2);

for i=1:numel(alphabeticalOrder)
    for j=1:numel(alphabeticalOrder)
        if alphabeticalOrder(i)==j
            positionInOrder(j)=i;
        end
    end
end

whichPHint=[];
for i=positionInOrder
    whichPHint=[whichPHint; subGroupInInt(i)+whichPH];
end

end

