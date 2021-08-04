function TF = isrow(x)
TF = false;
if ndims(x) == 2
    [r, c] = size(x);
    if c == 1
        TF = true;
    end
end

