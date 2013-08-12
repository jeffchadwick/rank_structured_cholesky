function write_sparse_matrix( A, filename )
  [i,j,s] = find(A);

  % Construct column pointers
  p = int32(zeros( size(A,2) + 1, 1 ));
  for col_idx = 1:size(A,2)
    p(col_idx + 1) = p(col_idx) + nnz(A(:,col_idx));
  end

  % Make rows zero-indexed
  i = i - 1;

  % Write to the file
  fid = fopen(filename, 'w');

  % This is a hack (using uint64 in place of size_t)
  fwrite(fid, nnz(A), 'uint64');
  fwrite(fid, size(A,1), 'uint64');
  fwrite(fid, size(A,2), 'uint64');

  % Write column offsets
  fwrite(fid, p(:,1), 'int32');

  % Row indices
  fwrite(fid, i(:,1), 'int32');

  % Numerical data
  fwrite(fid, s(:,1), 'double');
  
  fclose(fid);
end
