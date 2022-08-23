function output= align_vectors_to_learning(input, learning_ind, length_to_pad)

[dim1, dim2]=size(input);

if dim2>dim1
    input=input';
end

            % create a vector +/- sessions around learning (length_to_pad)
            input_to_pad= [NaN(1,length_to_pad)'; input; NaN(1,length_to_pad)'];
            learning_ind_to_pad=learning_ind+length_to_pad;
            output=input_to_pad(learning_ind_to_pad-length_to_pad:learning_ind_to_pad+length_to_pad);

end