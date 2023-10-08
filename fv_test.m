prev = 1;

for i=1:4
    [error1, x1, y1] = hw1_fv(10*i, 10*i);
    max_int_error = max(error1(:));
    fprintf("curr_err=%d | prev_err=%d\n", max_int_error, prev);

    prev = max_int_error;
end