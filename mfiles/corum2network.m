function [corum_network, unqprots] = corum2network(corum_filename)

% Read raw CORUM file
fid = fopen(corum_filename,'r');
fgetl(fid); % header

% read body of file
cc = 0;
organism = cell(10000,1);
complexes = cell(10000,1);
while ~feof(fid)
    t1 = strsplit(fgetl(fid),'\t');
    cc = cc+1;
    organism{cc} = t1{3};
    complexes{cc} = t1{6};
end
organism = organism(1:cc);
complexes = complexes(1:cc);
fclose(fid);


% Just human complixes with >=2 members
I = ismember(lower(organism), 'human') & contains(complexes, ';');
complexes = complexes(I);


% All unique proteins
cc = 0;
unqprots = cell(10^5,1);
for ii = 1:length(complexes)
    tmp = strsplit(complexes{ii},';');
    for jj = 1:length(tmp)
        cc = cc+1;
        unqprots{cc} = tmp{jj};
    end
end
unqprots = unique(unqprots(1:cc));


% Split complexes into pairwise list
corum_network = zeros(length(unqprots));
for ii = 1:length(complexes)
  cmplx = complexes{ii};
  cmplx = strrep(cmplx,',',' ');
  cmplx = strrep(cmplx,';',' ');
  cmplx = strsplit(cmplx,' ');
  cmplx = cmplx(not(cellfun('isempty',cmplx)));
  cmplx = unique(cmplx);
  Nprot = length(cmplx);

  if Nprot<2
    continue
  end
  
  I = find(ismember(unqprots, cmplx));
  for jj = 1:length(I)
      for kk = 1:length(I)
          if jj==kk; continue; end
          corum_network(I(jj),I(kk)) = 1;
      end
  end
end


