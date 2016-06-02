
% Some taken from Sepideh

N_file = 'analysis/diffusion_input_network_sparse.tsv';
S_file = 'analysis/diffusion_input_scores.tsv';

BETA = linspace(0, 0.03, 10)

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

A=diag(sum(WGGs));
L=WGGs-A;

disp('Computing eigenvalues');
tic
[T,D] = eig(L, 'nobalance');
toc

save('analysis/eig.mat', 'T', 'D')

disp('Computing inverse');
tic
invT  = inv(T);
save('analysis/inv.mat', 'invT')
toc

for beta = BETA

  beta = 0;

  disp(sprintf('Performing diffusion for beta=%f', beta));

  tic
    temp1                = beta*diag(D);
    temp1(temp1<(-2^10)) = -2^10;
    K                    = T*diag(exp(temp1))*invT;
    K(K< 0)              = 0;
    %% Diffusion Score
    SD                   = S*K;
  toc;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save
  outf_network = sprintf('analysis/diffusion_output_network.%0.5f.tsv', beta);
  outf_scores  = sprintf('analysis/diffusion_output_scores.%0.5f.tsv', beta);
  outf_mat = sprintf('analysis/diffusion_output_network_and_scores.%0.5f.mat', beta);

  save(outf_mat, 'K', 'SD')

  %save(outf_scores, 'SD', '-ascii');

  %[i,j,val]           = find(K);
  %[i,j,val] = find(WGGs);
  %network_sparse_dump = [i,j,val];

  %fd = fopen(outf_network,'w')
  %fprintf( fd,'%d\t%d\t%0.5f\n', network_sparse_dump' );
  %fprintf( fd,'%d\t%d\t%0.5f', nGenes, nGenes, 0);
  %fclose(fd);

end





