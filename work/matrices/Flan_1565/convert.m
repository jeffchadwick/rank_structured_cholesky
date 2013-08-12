addpath('../../../matlab');

load('Flan_1565.mat');

write_sparse_matrix( Problem.A, 'Flan_1565_system.bcsm' );

z = ones( size(Problem.A, 1), 1 );

write_vector( z, 'Flan_1565_rhs.vector' );

quit;
