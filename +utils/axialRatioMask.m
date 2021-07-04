function axialRatioMask(f_low, f_high)
    x = [f_low, f_low, f_high, f_high   ];
    y = [0,     3,     3,      0        ];
    
    fill(x, y, 'g-');
    alpha(0.2);
end