function [outMatch, outIdx]=tomtom(Q, T)
range=100;
QSize=size(Q);
numCells=QSize(2);
query_len=QSize(1);
targets_len=numel(T)/numCells;

target_count=targets_len/query_len;

QE=repmat(Q, 1,targets_len);

TV=T.';
TV=(TV(:)).';
TVE=repmat(TV, query_len, 1);


EDT=(QE-TVE).^2;
EDT=EDT.';
scoreMatrix=reshape(EDT, numCells, targets_len, query_len);
scoreMatrix=squeeze(sum(scoreMatrix, 1));
scoreMatrix=-sqrt(scoreMatrix.');


index_1=floor(targets_len/2);

index_2=ceil(targets_len/2);

scoreMatrix_sorted=sort(scoreMatrix, 2);

quantile=-(scoreMatrix_sorted(:, index_1)+scoreMatrix_sorted(:, index_2))/2;

scoreMatrixShifted=scoreMatrix+quantile;
large=max(max(scoreMatrixShifted));
small=min(min(scoreMatrixShifted));

scores_offset=floor(small);

small=floor(small);

if large-small>0
    scores_scale=floor(range/(large-small));
else
    scores_scale=floor(range/2);
end




scoreMatrixShifted=round(scores_scale*(scoreMatrixShifted-scores_offset));



% TQDScQuz=TQDScQuz.';
% 
% seedCellTi=scoreMatrixShifted(:)==(1:range);
% PWM1T=seedCellTi;
% PWM1T=reshape(PWM1T,query_len,targets_len, range);
% freqs=squeeze(sum(PWM1T, 2))/targets_len;
% freqs=[zeros(query_len, 1),freqs];
% freqs=freqs/targets_len;



address_index=1;
scores=scoreMatrixShifted;
uniform_bg=1/targets_len;
pv_lookup_matrix=zeros((query_len*(query_len + 1))/2, query_len*range+1);
pmf_matrix=zeros((query_len*(query_len + 1))/2, query_len*range+1);
reference_matrix=zeros(query_len, query_len);

  for start = query_len:-1:1
      pdf_new=zeros(1, query_len*range+1);
pdf_new(1)=1;

    for i = 1: query_len-start+1
      maxv = (i-1) * range;
      pdf_old=pdf_new;
      pdf_new=zeros(1, query_len*range+1);
      kv=(1: maxv+1);
      for j = 1: targets_len
        s = scores(start + i-1, j); 
        pdf_new(kv+s) = pdf_new(kv+s) + (pdf_old(kv) * uniform_bg);

      end


      pv=cumsum(pdf_new(end:-1:1));
      pv=pv(end:-1:1);

      pv_lookup_matrix(address_index, :)= pv;

      pmf_matrix(address_index,:)= pdf_new;


      reference_matrix(start, start + i-1)= address_index;
      address_index=address_index+1;
    end
  end
  idx=1;


% check=1;

  overlap_num=(2*query_len-1);

  temp_score=zeros(target_count*overlap_num,1);

  pv_index=zeros(target_count*overlap_num,1);

  i_start=(query_len:-1:query_len-overlap_num+1).';
  i_start=max(i_start, 1);
  i_start=repmat(i_start, target_count,1);

  i_stop=(overlap_num:-1:1).';
  i_stop=min(i_stop, query_len);
  i_stop=repmat(i_stop, target_count,1);

  overlap=i_stop-i_start+1;

for in=1:target_count

    scoret=scores(:, (in-1)*query_len+1:in*query_len);
    for idiag=-query_len+1:query_len-1
        temp_score(idx)=sum(diag(scoret,idiag));
        pv_index(idx) = reference_matrix(i_start(idx), i_stop(idx));
        idx=idx+1;

    end
    
end

pvalue=-(temp_score  + overlap * scores_offset * scores_scale);

% pvalue_notComplete = pv_lookup_matrix(pv_index, temp_score);

target_index=repmat((1:target_count),overlap_num, 1) ;
target_index=target_index(:);


[pvalueOverlapSorted, idxo]=sortrows([target_index, pvalue,-overlap]);

optimalPvalue=pvalueOverlapSorted(1:overlap_num:end,2);
optimalOverlap=-pvalueOverlapSorted(1:overlap_num:end,3);
optimalOffset=mod(idxo(1:overlap_num:end)-1, overlap_num)+1-query_len;

match.pvalue=-optimalPvalue;
match.evalue=zeros(size(optimalPvalue));

match.lengths=query_len*ones(target_count,1);
match.overlap=optimalOverlap;
match.offset=optimalOffset;
  check=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  norient = 1;

  pmf_len = size(pmf_matrix,2);

  ext_pmf_len = ceil(pmf_len - (query_len-1) * scores_offset * scores_scale);


  num_pv_lookup_array = size(pmf_matrix,1);



%   extended_pmf_matrix = zeros(num_pv_lookup_array, ext_pmf_len);
pmf_indexs=reference_matrix.';
pmf_indexs=pmf_indexs(:);
pmf_indexs=pmf_indexs(pmf_indexs>0);
start_colr=repelem((1: query_len).', (query_len:-1:1));
end_colr1=repmat((1: query_len).',query_len, 1);
end_colr0=repmat((1: query_len),query_len, 1);
end_colr0=end_colr0(:);


end_colr=end_colr1(end_colr1>=end_colr0);

offset_correctionN=scores_offset * (query_len - 1 - (end_colr-start_colr)) * scores_scale;
  extended_pmf_matrix = zeros(num_pv_lookup_array, ext_pmf_len);
  jV=(1:pmf_len);

  [~, idxSt]=sort(pmf_indexs);

offset_correctionN=offset_correctionN(idxSt);

for ij=1:length(pmf_indexs)

  extended_pmf_matrix(ij, jV-offset_correctionN(ij))=pmf_matrix(ij, jV);
end

% check=1;





  pmf_new = zeros(ext_pmf_len,1);

tlen.nuinique=1;
tlen.lengths=query_len*ones(target_count,1);
[~, ~, tlen.lookup]=unique(tlen.lengths);
pv_of_max_algn_matrix=zeros(tlen.nuinique, ext_pmf_len);
for j=1:tlen.nuinique
    pmf_new(1)= 1.0;
    len = tlen.lengths(j);
    nalgns = query_len + len - 1;

    for i =0: nalgns-1
      if (i < query_len) 
        start_col = query_len - i ;
        end_col = min(query_len, start_col+len) ;
      else 
        start_col = 1;
        end_col = min(query_len, len + query_len - i - 1) ;
      end





      pmf_index = reference_matrix(start_col, end_col);
      pmf_current = (extended_pmf_matrix(pmf_index,:)).';
      cdf_current=cumsum(pmf_current);

      pmf_old=pmf_new;
      cdf_old=cumsum(pmf_old);

      pmf_new(1)= pmf_current(1) * pmf_old(1);
      pmf_new(2:end)=pmf_current(2:end).*cdf_old(1:end-1)+pmf_old(2:end).*cdf_current(1:end-1)+pmf_current(2:end).*pmf_old(2:end);


      
    end

    pmf_old=cumsum(pmf_new(end:-1:1));  
    pmf_old=pmf_old(end:-1:1);  

    pv_of_max_algn_matrix(j,:)= pmf_old.';

end

  for i = 1:target_count 
    score = match.pvalue(i) - query_len * scores_offset * scores_scale;
    truncated_pvalue = pv_of_max_algn_matrix(tlen.lookup(i), score+1);
    match.pvalue(i) =  truncated_pvalue;
    match.evalue(i) = truncated_pvalue * target_count / norient;
  end

  [outMatch, outIdx]=sortrows([match.pvalue, -match.overlap]);
outMatch(:, 2)=-outMatch(:, 2);

outMatch=[outMatch, match.offset(outIdx)];



% 
% function An=compRow(Ap, freqsi)
% tQuantz=length(freqsi);
% An=zeros(1,tQuantz);
% for x=1:tQuantz
%     A1x=Ap(x-1:-1:1).*freqsi(1:x-1);
%     An(x)=sum(A1x);
% 
% end 

