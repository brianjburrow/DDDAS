disp("Checking HDMR first order expansion")
n = 100;
counter = 1;
for iVar1 = 1:n
    counter = counter + 1;
end
disp(counter);
disp(countEvals(100, 2, 1))

disp("Checking HDMR second order expansion")
n = 100;
counter = 1;
for iVar1 = 1:n
    for iVar2 = iVar1:n
        counter = counter + 1;
    end
end
disp(counter);
disp(countEvals(100, 2, 2))

disp("Checking HDMR third order expansion")
n = 100;
counter = 1;
for iVar1 = 1:n
    for iVar2 = iVar1:n
        for iVar3 = iVar2:n
            counter = counter + 1;
        end
    end
end


function count = countEvals(n, s, l)
    count = 1;
    for iOrder = 1:l
        count = count + nchoosek(n, iOrder)*(s - 1)^iOrder;
    end
end


