#!/usr/bin/env python

import argparse
import logging
import time
import uuid

from pathlib import Path


import pandas as pd
import numpy as np
import Levenshtein as ed


ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
ch.setFormatter(formatter)
logger = logging.getLogger('maund.'+str(uuid.uuid1())[2:8])
logger.setLevel(logging.DEBUG)
logger.addHandler(ch)


write_to_tsv = lambda df, output : df.to_csv(output,header=False,index=False,sep='\t')
toCharArray = lambda seq : np.array([seq]).astype('S').view('S1')
revertedSeq = lambda seq : seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
def mismatch(seq1,seq2):
    n = min(len(seq1),len(seq2))
    m = max(len(seq1),len(seq2))
    mismatch = np.sum(toCharArray(seq1[:n])!=toCharArray(seq2[:n]))
    return mismatch + m-n

def matchUpto1(target, seq):
    half_len = int(len(target)/2)
    fst = target[:half_len]
    snd = target[half_len:]
    
    i_fst = seq.find(fst)    
    while(i_fst!=-1):
        beg=i_fst
        end=beg+len(target)
        alignedSeq = seq[beg:end]
        if mismatch(alignedSeq, target)<2:
            return beg
        i_fst = seq.find(fst,end)
    
    i_snd = seq.find(snd)
    while(i_snd!=-1):
        beg=i_snd-half_len
        end=beg+len(target)
        alignedSeq = seq[beg:end]
        if mismatch(alignedSeq, target)<2:
            return beg
        i_snd = seq.find(snd,end)
    return -1


parser = argparse.ArgumentParser(description='MAUND')
#parser.add_argument("-v", "--verbosity", action="count", default=0,
#                    help="increase output verbosity")
parser.add_argument('aseq')
parser.add_argument('rgen')
parser.add_argument('files', nargs='*')
parser.add_argument('-c','--comparison_range', type=int, default=60)
parser.add_argument('-b','--window_beg', type=int, default=4)
parser.add_argument('-e','--window_end', type=int, default=7)
parser.add_argument('-ib','--idxseq_beg', type=int, default=13)
parser.add_argument('-ie','--idxseq_end', type=int, default=22)
parser.add_argument('-t','--target_nt',  default='A', choices=['A','C','G','T'])
parser.add_argument('-mcut','--mismatch_cutoff', type=int, default=4)
args = parser.parse_args()

t1=time.time()
seq_wt,RGEN_seq, f_input = args.aseq, args.rgen, args.files
seq_wt=seq_wt.upper()     # target seq of 250 bp
RGEN_seq=RGEN_seq.upper() # 23bp for RGEN

i_window_beg = args.window_beg - 1
i_window_end = args.window_end
i_idxseq_beg = args.idxseq_beg - 1
i_idxseq_end = args.idxseq_end
nt_target    = args.target_nt

mismatch_cutoff = args.mismatch_cutoff

len_rgen = len(RGEN_seq)
idx_cleavage = len_rgen-6
len_crange = args.comparison_range
filt_n = 1                                       # cutoff for count number
len_indicator_seq = 15

for file_path in f_input:
    file_path = Path(file_path).resolve()
    file_name = file_path.name
    logger.info("Begin: {} with {}".format(file_name, RGEN_seq))
    out_name= '{}.{}.maund.out.'.format(file_name,RGEN_seq)
        
    i_for =seq_wt.find(RGEN_seq)
    i_rev = seq_wt.find(revertedSeq(RGEN_seq))
    if i_for != -1 : 
        i = i_for
        start_pos=i+idx_cleavage-len_crange
        end_pos  =i+idx_cleavage+len_crange        
        is_rev_match = False
    elif i_rev != -1 : 
        i = i_rev
        start_pos=i+(len_rgen-idx_cleavage)-len_crange
        end_pos  =i+(len_rgen-idx_cleavage)+len_crange        
        is_rev_match = True
    else:
        logger.error('Cannot find target seqence in amplicon sequence.')
        raise Exception('Cannot find target seqence in amplicon sequence.')
    if start_pos < 0:
        start_pos=0
    if end_pos > len(seq_wt):
        end_pos = len(seq_wt)
    s_seq = seq_wt[i:i+len_rgen]
    seq_range = seq_wt[start_pos:end_pos]
    
    pri_for = seq_range.upper()[:len_indicator_seq]                   # 15bp primer < 1bp mismatch
    pri_back = seq_range.upper()[-len_indicator_seq:]                 # 15bp primer < 1bp mismatch
    length_range=len(seq_range)
    
    logger.info("Is reverted match? : {}".format(is_rev_match))
    logger.info("Use amplicon[{}:{}] as comparison range".format(start_pos,end_pos))
    ###################################################  0
    with file_path.open("r") as f1 :
        s1 = f1.read().splitlines()
    # Get seqence lines
    seq_lines = s1[1::4]
    s_N = 'N'
    
    df_z1=pd.DataFrame(seq_lines, columns=['frag'])
    df_z2=df_z1[-df_z1.frag.str.contains(s_N)].copy()
    df_z2['i_beg']=df_z2.frag.apply(lambda x : matchUpto1(pri_for,x))
    tmp = df_z2[df_z2.i_beg!=-1].copy()
    tmp['i_end']=tmp.frag.apply(lambda x : matchUpto1(pri_back,x))
    df_z3 = tmp[tmp.i_beg<tmp.i_end].apply(lambda x : x.frag[x.i_beg:x.i_end+len(pri_back)], axis=1).copy()
    df_z4 = df_z3.value_counts()
    
    #write_to_tsv(df_z1, out_name+"z1_M1.fastq.txt")
    #write_to_tsv(df_z2, out_name+"z2_M1.fastq.txt")
    #write_to_tsv(df_z3, out_name+"z3_M1.fastq.txt")
    #write_to_tsv(df_z4.reset_index(), out_name+"z4_re_M1.fastq.txt")    #z4 is Series
    
    # Write output_all
    output_window  = out_name+"_window.txt"
    output_aligned = out_name+"_aligned.txt"    

    #Filter    
    df_seq_all = df_z4.reset_index()
    df_seq_all.columns = ['seq','n_seq']
    df_seq_all['seq_len']=df_seq_all.seq.str.len()
    df_seq_all=df_seq_all[df_seq_all.n_seq>filt_n]    
    df_seq_mut      = df_seq_all[(-df_seq_all.seq.str.contains(s_seq))&(df_seq_all.seq_len!=length_range)]
    df_seq_wt       = df_seq_all[(-df_seq_all.seq.str.contains(s_seq))&(df_seq_all.seq_len==length_range)]
    df_seq_same_len = df_seq_all[df_seq_all.seq_len==length_range]    
    
    write_to_tsv(df_seq_all, out_name+"_all.txt")   
    write_to_tsv(df_seq_mut, out_name+"_mut.txt")
    write_to_tsv(df_seq_wt,  out_name+"_WT_subst.txt")
    write_to_tsv(df_seq_same_len, out_name+"_same_length.txt")
    
    def col_sum(col):
        if col.empty:
            return 0
        return col.sum()
    tot_count      = col_sum(df_seq_all.n_seq)
    mut_count      = col_sum(df_seq_mut.n_seq)
    wt_subst_count = col_sum(df_seq_wt.n_seq)
    same_len_count = col_sum(df_seq_same_len.n_seq)       

    f2_cols     = "#target_seq No_A No_C No_G No_T No_Non_Target ratio_A ratio_C ratio_G ratio_T ratio_non_target t_ratio_A t_ratio_C t_ratio_G t_ratio_T ratio_non_target_to_allWT".split()    
    f2_cols_rev = "#target_seq No_T No_G No_C No_A No_Non_Target ratio_T ratio_G ratio_C ratio_A ratio_non_target t_ratio_T t_ratio_G t_ratio_C t_ratio_A ratio_non_target_to_allWT".split()    
    
    wt_len_total=tot_count-mut_count
    
    if is_rev_match : 
        rgen_cridx_beg = len_crange + idx_cleavage - len_rgen
        rgen_cridx_end = len_crange + idx_cleavage
    else :
        rgen_cridx_beg = len_crange - idx_cleavage
        rgen_cridx_end = len_crange - idx_cleavage + len_rgen
        
    def countPerLocation2(df_target_regions, s_seq):
        zero = lambda nt : {'targetRegion':nt*len(s_seq), 'n_seq':0}
        df_reads = df_target_regions[['targetRegion','n_seq']].append([zero('A'),zero('C'),zero('G'),zero('T')], ignore_index=True)

        df_sseq = pd.DataFrame(np.array([s_seq]).astype('S%d'%len_rgen).view('int8'))
        df = pd.DataFrame(df_reads.targetRegion.astype('S%d'%len_rgen).values.view('int8').reshape(-1,len_rgen))
        df.index = df_reads.index
        #df2 : idx nc n_seq
        df2 = df.stack().reset_index(level=1).merge(df_reads[['n_seq']],left_index=True,right_index=True)
        df2.columns = ['idx', 'nc', 'n_seq']
        df3 = pd.pivot_table(df2, values='n_seq', index=['idx'], columns=['nc'], aggfunc=np.sum, fill_value=0)
        df_subst = df3.unstack().drop(list(zip(*[df_sseq[0],df_sseq.index])), errors='ignore').groupby(level=1).sum()
        df_count = pd.concat([df3,df_subst],axis=1)
        return df_count

    def countPerLocation(df_reads, is_rev_match, s_seq):
        a = rgen_cridx_beg
        b = rgen_cridx_end
        df = df_reads.copy()
        df['targetRegion'] = df_reads.seq.str.slice(a,b)    
        return countPerLocation2(df, s_seq)

    seq_wt_count  = countPerLocation(df_seq_wt, is_rev_match, s_seq)
    seq_wt_ratio = seq_wt_count/wt_subst_count
    samelen_ratio = countPerLocation(df_seq_same_len, is_rev_match, s_seq)/same_len_count

    df=pd.concat([pd.DataFrame(list(s_seq)), seq_wt_count.astype('int'), seq_wt_ratio.round(4), samelen_ratio.round(4)],axis=1)
    df.columns = f2_cols
    if is_rev_match:
        df = df[::-1]
        df["#target_seq"] = df["#target_seq"].map({"A":"T", "T":"A", "C":"G","G":"C"})
        df.columns = f2_cols_rev
        df = df[f2_cols]
    df.to_csv(out_name+"_subst_result.txt",index=False,sep='\t')

    
    def hasReplacement(nt, ref, seq):
        for r,s in zip(ref,seq):
            if r==nt and r!=s:
                return True
        return False
    
    #fineMatch : edit distance to target sequence is less than cutoff
    #well-aligned : no indel and find match to target sequence.
    noIndelInTargetSeqRegion = lambda ops : 0 == len([op for op in ops if op[0]!='replace' and rgen_cridx_beg<=op[1] and op[1]<rgen_cridx_end])
    isFineMatch = lambda ops : len([op for op in ops if rgen_cridx_beg<=op[1] and op[1]<rgen_cridx_end]) < mismatch_cutoff
    
    def findIdx(editops, idx):
        ins  = [op for op in editops if op[0]=='insert' and op[1]<idx]
        dels = [op for op in editops if op[0]=='delete' and op[1]<idx]
        return idx + len(ins)-len(dels)
    
    df_all = df_seq_all.copy()
    df_all['editops']=df_all.seq.apply(lambda x : ed.editops(seq_range,x))
    df_fineMatch   = df_all[df_all.editops.apply(isFineMatch)].copy()
    df_wellAligned = df_fineMatch[df_fineMatch.editops.apply(noIndelInTargetSeqRegion)].copy()    
    
    i_rgen = seq_range.find(s_seq)
    logger.info("Target seq range is located in comparison_range[{}:{}]".format(i_rgen,i_rgen+len_rgen))
    if is_rev_match :
        w_beg = i_rgen + len(s_seq) - i_window_end
        w_end = i_rgen + len(s_seq) - i_window_beg
    else:
        w_beg = i_rgen + i_window_beg
        w_end = i_rgen + i_window_end
    logger.info("{} base editing in window : comparison_range[{}:{}] = {}".format(nt_target,w_beg,w_end,seq_range[w_beg:w_end]))
    if df_wellAligned.empty :
        logger.info("No well-aligned case. Skip to generate _window and _aligned")
        n_total = 0
        n_mutated = 0
        mutation_ratio = np.nan
    else :
        df_wellAligned['window']=df_wellAligned.apply(lambda x : x.seq[findIdx(x.editops,w_beg):findIdx(x.editops,w_end)], axis=1)
        df_Win = df_wellAligned[['window','n_seq']].groupby('window').sum().sort_values(by='n_seq',ascending=False).reset_index()
        if is_rev_match:
            df_Win.window=df_Win.window.apply(revertedSeq)
        df_Win.to_csv(output_window,sep='\t',index=False)
    
        rgen_window = RGEN_seq[i_window_beg:i_window_end]
        n_mutated = col_sum(df_Win[df_Win.window.apply(lambda x : hasReplacement(nt_target,rgen_window,x))].n_seq)
        n_total   = col_sum(df_Win.n_seq)
        mutation_ratio = n_mutated/n_total
    
    #targetRegion : region that corresponds to the target sequence region
        tr_beg = i_rgen
        tr_end = i_rgen + len(s_seq)
        df_wellAligned['targetRegion']=df_wellAligned.apply(lambda x : x.seq[findIdx(x.editops,tr_beg):findIdx(x.editops,tr_end)], axis=1)
        counts = countPerLocation2(df_wellAligned[df_wellAligned.targetRegion.str.len()==len(RGEN_seq)], s_seq)
        n_total = df_wellAligned.n_seq.sum()
        counts=pd.concat([pd.DataFrame(list(s_seq),columns=['target_seq']),counts,counts/n_total],axis=1).round(4)
        cols_out     = "target_seq No_A No_C No_G No_T No_Non_Target ratio_A ratio_C ratio_G ratio_T ratio_non_target".split()
        cols_out_rev = "target_seq No_T No_G No_C No_A No_Non_Target ratio_T ratio_G ratio_C ratio_A ratio_non_target".split()
        if is_rev_match:
            counts = counts[::-1]
            counts.target_seq = counts.target_seq.map({"A":"T", "T":"A", "C":"G","G":"C"})
            counts.columns = cols_out_rev
            counts = counts[cols_out]
        else:
            counts.columns = cols_out
        counts.to_csv(output_aligned,index=False,sep='\t')
    
    if is_rev_match:
        index_beg = i_rgen + len(s_seq) - i_idxseq_end
        index_end = i_rgen + len(s_seq) - i_idxseq_beg
    else:
        index_beg = i_rgen + i_idxseq_beg
        index_end = i_rgen + i_idxseq_end
    index_seq = seq_range[index_beg:index_end]
    logger.info("Index sequence = comparison_range[{}:{}] = {}".format(index_beg,index_end,index_seq))
    df_indels = df_seq_all[(-df_seq_all.seq.str.contains(index_seq))&(df_seq_all.seq_len!=len(seq_range))]
    n_indels = col_sum(df_indels.n_seq)
    n_all    = df_seq_all.n_seq.sum()
    indel_ratio = n_indels/n_all    
    
    with open(out_name+"Miseq_summary.txt",'w') as fsumm:
        #fsumm.write('{}\t{}\t{}\t{}\t{.4f}\n'.format("input_file", "target_seq", "window_mutated", "window_total", "window_ratio"))
        form = '{}\t{}\t{}\t{}\t{:.4f}\t{}\t{}\t{:.4f}\n'
        fsumm.write(form.format(file_name, RGEN_seq, n_mutated,n_total, mutation_ratio, n_indels, n_all, indel_ratio))
    
t2=time.time()
logger.info('Finished: {:.4f} sec'.format(t2-t1))

