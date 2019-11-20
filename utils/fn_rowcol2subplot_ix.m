function subplot_ix = fn_rowcol2subplot_ix(n_row,n_col,row_ix,col_ix)
% Return the subplot index of a certain (row,col) position in a subplot
% array of size (n_rows,n_cols)

% fill by column, so do that then transpose
subplot_idx = reshape(1:n_row*n_col,[n_col n_row])';
subplot_ix = subplot_idx(row_ix,col_ix);

end
