function y = bitxor_toni( array )

if size(array,2) == 1
    y = array;
else
    y = bitxor(array(1),bitxor_toni(array(2:end)));
end