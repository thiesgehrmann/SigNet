
% Some taken from Sepideh

beta = beta %Parameter should be given!!!

N_file = 'analysis/diffusion_input_network_sparse.tsv';
S_file = 'analysis/diffusion_input_scores.tsv';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Reading the data');
load analysis/diffusion_input_network_sparse.tsv;
WGGs = spconvert(diffusion_input_network_sparse);
WGGs = full(WGGs);

nGenes = max(size(WGGs))

fd = fopen(S_file);
S= fscanf(fd, '%f');
fclose(fd);
S=S';

disp('Done reading data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = WGGs;

D = sum(A)';
Laplacian = diag(D) - A;
disp('Computing eigenvalues');
tic
[T, DIAG] = eig(-Laplacian);
toc
DIAG = diag(DIAG);


disp('Computing inverse');
tic
Tinv = inv(T);
toc

disp(sprintf('Performing diffusion for beta=%f', beta));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
  temp1                = beta*DIAG;
  temp1(temp1<(-2^10)) = -2^10;
  K                    = T*diag(exp(temp1))*Tinv;
  K(K< 0)              = 0;
  %% Diffusion Score
  SD                   = S*K;
  SD                   = SD'
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Save
outf_network = sprintf('analysis/diffusion_output_network.%0.5f.tsv', beta);
outf_scores  = sprintf('analysis/diffusion_output_scores.%0.5f.tsv', beta);
outf_mat     = sprintf('analysis/diffusion_output_network_and_scores.%0.5f.mat', beta);

fd = fopen(outf_scores, 'w');
fprintf(fd, '%0.7f\n', SD);
fclose(fd);
disp('done writing scores');

[i,j,val]           = find(K > 10^-10);
network_sparse_dump = [i,j,val];

fd = fopen(outf_network,'w')
fprintf( fd,'%d\t%d\t%0.7f\n', network_sparse_dump' );
fprintf( fd,'%d\t%d\t%0.7f', nGenes, nGenes, K(nGenes, nGenes));
fclose(fd);
disp('done writing network')

