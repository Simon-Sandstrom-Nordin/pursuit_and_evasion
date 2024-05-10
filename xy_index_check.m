function index = xy_index_check(index_to_be_checked, N_matrix)
    if index_to_be_checked <= 0
        index = 1;
    elseif index_to_be_checked > N_matrix
        index = N_matrix;
    else
        index = index_to_be_checked;
    end
end
