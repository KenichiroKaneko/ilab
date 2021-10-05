A = [1, 2, 3, 4, 6];

for i = 1:6

    if (find(A == i))
        disp(i)
    else
        disp('hello');
    end

end
