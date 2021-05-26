function cubicPoly(x)
    y = x.^3 + x.^2 + x + 1; 
    return y
end

q = (cubicPoly,0,1)