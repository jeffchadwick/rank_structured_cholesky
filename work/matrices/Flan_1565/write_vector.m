function v = write_vector( x, fname )
  x = double( x );
  fid = fopen (fname, 'w');
  n_col = 1;
  n_row = size(x, 1);
  fwrite (fid, n_row, 'int32');
  fwrite (fid, x(:,1), 'double');
  fclose (fid);
end
