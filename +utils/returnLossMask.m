function returnLossMask(f_low, f_high)
    x = [f_low f_low f_high f_high ];
    y = [0     10    10     0      ];
    fill(x, y, 'Green');
    alpha(0.5);
end