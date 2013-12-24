% Teager function to remove noise
function y = deco_teager(x)
l = 1;
p = 2;
q = 0;
s = 3;
for k = 4:size(x,1)
    %y(k,:) = x(k,:).*x(k,:) - x(k-1,:).*x(k-2,:);
    y(k,:) = x(k-l,:).*x(k-p,:) - x(k-q,:).*x(k-s,:);
end

end