## #
## #  =====================================================================================
## # 
## #        Filename:  fastcon.c
## # 
## #     Description:
## # 
## #         Version:  0.1
## #         Created:  07/20/2013 17:00:00
## #        Revision:  none
## #        Compiler:  gcc
## # 
## #          Author:  Jason Chin,
## #         Company:
## # 
## #  =====================================================================================
## #
## # #################################################################################$$
## # # Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
## # #
## # # All rights reserved.
## # #
## # # Redistribution and use in source and binary forms, with or without
## # # modification, are permitted (subject to the limitations in the
## # # disclaimer below) provided that the following conditions are met:
## # #
## # #  * Redistributions of source code must retain the above copyright
## # #  notice, this list of conditions and the following disclaimer.
## # #
## # #  * Redistributions in binary form must reproduce the above
## # #  copyright notice, this list of conditions and the following
## # #  disclaimer in the documentation and/or other materials provided
## # #  with the distribution.
## # #
## # #  * Neither the name of Pacific Biosciences nor the names of its
## # #  contributors may be used to endorse or promote products derived
## # #  from this software without specific prior written permission.
## # #
## # # NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
## # # GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
## # # BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
## # # WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
## # # OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
## # # DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
## # # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## # # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## # # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
## # # USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
## # # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
## # # OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
## # # OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
## # # SUCH DAMAGE.
## # #################################################################################$$
## # 

import
  common

const
  UINT8_MAX = 255
  UINT_MAX = uint32.high


type
  align_tag_t* = object
    t_pos*: seq_coor_t
    delta*: uint8
    q_base*: char
    p_t_pos*: seq_coor_t       ## # the tag position of the previous base
    p_delta*: uint8            ## # the tag delta of the previous base
    p_q_base*: char            ## # the previous base
    q_id*: uint32

  align_tags_t* = object
    length*: seq_coor_t
    align_tags*: ptr align_tag_t

  align_tag_col_t* = object
    size*: uint16
    n_link*: uint16
    p_t_pos*: ptr seq_coor_t    ## # the tag position of the previous base
    p_delta*: ptr uint8         ## # the tag delta of the previous base
    p_q_base*: cstring         ## # the previous base
    link_count*: ptr uint16
    count*: uint16
    best_p_t_pos*: seq_coor_t
    best_p_delta*: uint8
    best_p_q_base*: uint8      ## # encoded base
    score*: cdouble

  msa_base_group_t* = object
    base*: ptr align_tag_col_t

  msa_delta_group_t* = object
    size*: uint8
    max_delta*: uint8
    delta*: ptr msa_base_group_t

  msa_pos_t* = ptr msa_delta_group_t

proc get_align_tags*(aln_q_seq: cstring; aln_t_seq: cstring; aln_seq_len: seq_coor_t;
                    srange: ptr aln_range; q_id: uint32; t_offset: seq_coor_t): ptr align_tags_t =
  var p_q_base: char
  var tags: ptr align_tags_t
  var
    i: seq_coor_t
    j: seq_coor_t
    jj: seq_coor_t
    k: seq_coor_t
    p_j: seq_coor_t
    p_jj: seq_coor_t
  tags = calloc(1, sizeof((align_tags_t)))
  tags.length = aln_seq_len
  tags.align_tags = calloc(aln_seq_len + 1, sizeof((align_tag_t)))
  i = srange.s1 - 1
  j = srange.s2 - 1
  jj = 0
  p_j = - 1
  p_jj = 0
  p_q_base = '.'
  k = 0
  while k < aln_seq_len:
    if aln_q_seq[k] != '-':
      inc(i)
      inc(jj)
    if aln_t_seq[k] != '-':
      inc(j)
      jj = 0
    if j + t_offset >= 0 and jj < UINT8_MAX and p_jj < UINT8_MAX:
      (tags.align_tags[k]).t_pos = j + t_offset
      (tags.align_tags[k]).delta = uint8(jj)
      (tags.align_tags[k]).p_t_pos = p_j + t_offset
      (tags.align_tags[k]).p_delta = uint8(p_jj)
      (tags.align_tags[k]).p_q_base = p_q_base
      (tags.align_tags[k]).q_base = aln_q_seq[k]
      (tags.align_tags[k]).q_id = q_id
      p_j = j
      p_jj = jj
      p_q_base = aln_q_seq[k]
    inc(k)
  ## # sentinal at the end
  ## #k = aln_seq_len;
  tags.length = k
  (tags.align_tags[k]).t_pos = -1 #UINT_MAX
  (tags.align_tags[k]).delta = UINT8_MAX
  (tags.align_tags[k]).q_base = '.'
  (tags.align_tags[k]).q_id = UINT_MAX
  return tags

proc free_align_tags*(tags: ptr align_tags_t) =
  free(tags.align_tags)
  free(tags)

proc allocate_aln_col*(col: ptr align_tag_col_t) =
  col.p_t_pos = cast[ptr seq_coor_t](calloc(col.size, sizeof((seq_coor_t))))
  col.p_delta = cast[ptr uint8](calloc(col.size, sizeof((uint8))))
  col.p_q_base = cast[cstring](calloc(col.size, sizeof((char))))
  col.link_count = cast[ptr uint16](calloc(col.size, sizeof((uint16))))

proc realloc_aln_col*(col: ptr align_tag_col_t) =
  col.p_t_pos = cast[ptr seq_coor_t](realloc(col.p_t_pos,
      (col.size) * sizeof((seq_coor_t))))
  col.p_delta = cast[ptr uint8](realloc(col.p_delta, (col.size) * sizeof((uint8))))
  col.p_q_base = cast[cstring](realloc(col.p_q_base, (col.size) * sizeof((char))))
  col.link_count = cast[ptr uint16](realloc(col.link_count,
      (col.size) * sizeof((uint16))))

proc free_aln_col*(col: ptr align_tag_col_t) =
  free(col.p_t_pos)
  free(col.p_delta)
  free(col.p_q_base)
  free(col.link_count)

proc allocate_delta_group*(g: ptr msa_delta_group_t) =
  var
    i: cint
    j: cint
  g.max_delta = 0
  g.delta = cast[ptr msa_base_group_t](calloc(g.size, sizeof((msa_base_group_t))))
  i = 0
  while i < g.size:
    g.delta[i].base = cast[ptr align_tag_col_t](calloc(5, sizeof((align_tag_col_t))))
    j = 0
    while j < 5:
      g.delta[i].base[j].size = 8
      allocate_aln_col(addr((g.delta[i].base[j])))
      inc(j)
    inc(i)

proc realloc_delta_group*(g: ptr msa_delta_group_t; new_size: uint16) =
  var
    i: cint
    j: cint
    bs: cint
    es: cint
  bs = g.size
  es = new_size
  g.delta = cast[ptr msa_base_group_t](realloc(g.delta,
      new_size * sizeof((msa_base_group_t))))
  i = bs
  while i < es:
    g.delta[i].base = cast[ptr align_tag_col_t](calloc(5, sizeof((align_tag_col_t))))
    j = 0
    while j < 5:
      g.delta[i].base[j].size = 8
      allocate_aln_col(addr((g.delta[i].base[j])))
      inc(j)
    inc(i)
  g.size = new_size

proc free_delta_group*(g: ptr msa_delta_group_t) =
  ## #manything to do here
  var
    i: cint
    j: cint
  i = 0
  while i < g.size:
    j = 0
    while j < 5:
      free_aln_col(addr((g.delta[i].base[j])))
      inc(j)
    free(g.delta[i].base)
    inc(i)
  free(g.delta)

proc update_col*(col: ptr align_tag_col_t; p_t_pos: seq_coor_t; p_delta: uint8;
                p_q_base: char) =
  var updated: cint = 0
  var kk: cint
  inc(col.count, 1)
  kk = 0
  while kk < col.n_link:
    if p_t_pos == col.p_t_pos[kk] and p_delta == col.p_delta[kk] and
        p_q_base == col.p_q_base[kk]:
      inc(col.link_count[kk])
      updated = 1
      break
    inc(kk)
  if updated == 0:
    if col.n_link + 1 > col.size:
      if col.size < (UINT16_MAX > 1) - 1:
        col.size = col.size * 2
      else:
        inc(col.size, 256)
      assert(col.size < UINT16_MAX - 1)
      realloc_aln_col(col)
    kk = col.n_link
    col.p_t_pos[kk] = p_t_pos
    col.p_delta[kk] = p_delta
    col.p_q_base[kk] = p_q_base
    col.link_count[kk] = 1
    inc(col.n_link)

proc get_msa_working_sapce*(max_t_len: cuint): ptr msa_pos_t =
  var msa_array: ptr msa_pos_t
  var i: cuint
  msa_array = calloc(max_t_len, sizeof(ptr msa_pos_t))
  i = 0
  while i < max_t_len:
    msa_array[i] = calloc(1, sizeof((msa_delta_group_t)))
    msa_array[i].size = 8
    allocate_delta_group(msa_array[i])
    inc(i)
  return msa_array

proc clean_msa_working_space*(msa_array: ptr msa_pos_t; max_t_len: cuint) =
  var
    i: cuint
    j: cuint
    k: cuint
  var col: ptr align_tag_col_t
  i = 0
  while i < max_t_len:
    j = 0
    while j < msa_array[i].max_delta + 1:
      k = 0
      while k < 5:
        col = msa_array[i].delta[j].base + k
        ## #
        ## #                for (c =0; c < col->size; c++) {
        ## #                    col->p_t_pos[c] = 0;
        ## #                    col->p_delta[c] = 0;
        ## #                    col->p_q_base[c] = 0;
        ## #                    col->link_count[c] =0;
        ## #                }
        ## #                
        col.n_link = 0
        col.count = 0
        col.best_p_t_pos = 0
        col.best_p_delta = 0
        col.best_p_q_base = 0
        col.score = 0
        inc(k)
      inc(j)
    msa_array[i].max_delta = 0
    inc(i)

const
  STATIC_ALLOCATE* = true

## ##undef STATIC_ALLOCATE

proc get_cns_from_align_tags*(tag_seqs: ptr ptr align_tags_t; n_tag_seqs: cuint;
                             t_len: cuint; min_cov: cuint): ptr consensus_data =
  var
    i: seq_coor_t
    j: seq_coor_t
  var t_pos: seq_coor_t = 0
  var coverage: ptr cuint
  var local_nbase: ptr cuint
  var consensus: ptr consensus_data
  ## #char * consensus;
  var c_tag: ptr align_tag_t
  coverage = calloc(t_len, sizeof(cuint))
  local_nbase = calloc(t_len, sizeof(cuint))
  when not defined(STATIC_ALLOCATE):
    var msa_array: ptr msa_pos_t = nil
    msa_array = calloc(t_len, sizeof(ptr msa_pos_t))
    i = 0
    while i < t_len:
      msa_array[i] = calloc(1, sizeof((msa_delta_group_t)))
      msa_array[i].size = 8
      allocate_delta_group(msa_array[i])
      inc(i)
  else:
    var msa_array: ptr msa_pos_t = nil
    if msa_array == nil:
      var xxx: pid_t = syscall(__NR_gettid)
      fprintf(stderr, "STATIC allocating: %d\x0A", xxx)
      msa_array = get_msa_working_sapce(100000)
    assert(t_len < 100000)
  ## # loop through every alignment
  fprintf(stderr, "XX %d\x0A", n_tag_seqs)
  i = 0
  while i < n_tag_seqs:
    ## # for each alignment position, insert the alignment tag to msa_array
    j = 0
    while j < tag_seqs[i].length:
      c_tag = tag_seqs[i].align_tags + j
      var delta: cuint
      delta = c_tag.delta
      if delta == 0:
        t_pos = c_tag.t_pos
        inc(coverage[t_pos])
      if delta > msa_array[t_pos].max_delta:
        msa_array[t_pos].max_delta = delta
        if msa_array[t_pos].max_delta + 4 > msa_array[t_pos].size:
          realloc_delta_group(msa_array[t_pos], msa_array[t_pos].max_delta + 8)
      var base: cuint = - 1
      case c_tag.q_base
      of 'A':
        base = 0
      of 'C':
        base = 1
      of 'G':
        base = 2
      of 'T':
        base = 3
      of '-':
        base = 4
      ## # Note: On bad input, base may be -1.
      update_col(addr((msa_array[t_pos].delta[delta].base[base])), c_tag.p_t_pos,
                 c_tag.p_delta, c_tag.p_q_base)
      inc(local_nbase[t_pos])
      inc(j)
    inc(i)
  ## # propogate score throught the alignment links, setup backtracking information
  var g_best_aln_col: ptr align_tag_col_t = 0
  var g_best_ck: cuint = 0
  var g_best_t_pos: seq_coor_t = 0
  var kk: cint
  var ck: cint
  ## # char base;
  var best_i: cint
  var best_j: cint
  var best_b: cint
  var best_ck: cint = - 1
  var score: cdouble
  var best_score: cdouble
  var g_best_score: cdouble
  ## # char best_mark;
  var aln_col: ptr align_tag_col_t
  g_best_score = - 1
  i = 0
  while i < t_len:
    ## #loop through every template base
    ## #printf("max delta: %d %d\n", i, msa_array[i]->max_delta);
    j = 0
    while j <= msa_array[i].max_delta:
      ## # loop through every delta position
      kk = 0
      while kk < 5:
        ## # loop through diff bases of the same delta posiiton
        ## #
        ## #                    switch (kk) {
        ## #                        case 0: base = 'A'; break;
        ## #                        case 1: base = 'C'; break;
        ## #                        case 2: base = 'G'; break;
        ## #                        case 3: base = 'T'; break;
        ## #                        case 4: base = '-'; break;
        ## #                    }
        ## #                    
        aln_col = msa_array[i].delta[j].base + kk
        if aln_col.count >= 0:
          best_score = - 1
          best_i = - 1
          best_j = - 1
          best_b = - 1
          ck = 0
          while ck < aln_col.n_link:
            ## # loop through differnt link to previous column
            var pi: cint
            var pj: cint
            var pkk: cint
            pi = aln_col.p_t_pos[ck]
            pj = aln_col.p_delta[ck]
            case aln_col.p_q_base[ck]
            of 'A':
              pkk = 0
            of 'C':
              pkk = 1
            of 'G':
              pkk = 2
            of 'T':
              pkk = 3
            of '-':
              pkk = 4
            else:
              pkk = 4
            if aln_col.p_t_pos[ck] == - 1:
              score = cast[cdouble](aln_col.link_count[ck]) -
                  cast[cdouble](coverage[i] * 0.5)
            else:
              score = msa_array[pi].delta[pj].base[pkk].score +
                  cast[cdouble](aln_col.link_count[ck]) -
                  cast[cdouble](coverage[i] * 0.5)
            ## # best_mark = ' ';
            if score > best_score:
              best_score = score
              aln_col.best_p_t_pos = best_i = pi
              aln_col.best_p_delta = best_j = pj
              aln_col.best_p_q_base = best_b = pkk
              best_ck = ck
              ## # best_mark = '*';
            inc(ck)
          aln_col.score = best_score
          if best_score > g_best_score:
            g_best_score = best_score
            g_best_aln_col = aln_col
            g_best_ck = best_ck
            g_best_t_pos = i
            ## #printf("GB %d %d %d %d\n", i, j, ck, g_best_aln_col);
        inc(kk)
      inc(j)
    inc(i)
  assert(g_best_score != - 1)
  ## # reconstruct the sequences
  var index: cuint
  var bb: char = '$'
  var ck: cint
  var cns_str: cstring
  var eqv: ptr cint
  var score0: cdouble
  consensus = calloc(1, sizeof((consensus_data)))
  consensus.sequence = calloc(t_len * 2 + 1, sizeof((char)))
  consensus.eqv = calloc(t_len * 2 + 1, sizeof(cuint))
  cns_str = consensus.sequence
  eqv = consensus.eqv
  index = 0
  ck = g_best_ck
  i = g_best_t_pos
  while 1:
    if coverage[i] > min_cov:
      case ck
      of 0:
        bb = 'A'
      of 1:
        bb = 'C'
      of 2:
        bb = 'G'
      of 3:
        bb = 'T'
      of 4:
        bb = '-'
    else:
      case ck
      of 0:
        bb = 'a'
      of 1:
        bb = 'c'
      of 2:
        bb = 'g'
      of 3:
        bb = 't'
      of 4:
        bb = '-'
    ## # Note: On bad input, bb will keep previous value, possibly '$'.
    score0 = g_best_aln_col.score
    i = g_best_aln_col.best_p_t_pos
    if i == - 1 or index >= t_len * 2: break
    j = g_best_aln_col.best_p_delta
    ck = g_best_aln_col.best_p_q_base
    g_best_aln_col = msa_array[i].delta[j].base + ck
    if bb != '-':
      cns_str[index] = bb
      eqv[index] = cast[cint](score0) - cast[cint](g_best_aln_col.score)
      ## #printf("C %d %d %c %lf %d %d\n", i, index, bb, g_best_aln_col->score, coverage[i], eqv[index] );
      inc(index)
  ## # reverse the sequence
  i = 0
  while i < index div 2:
    cns_str[i] = cns_str[i] xor cns_str[index - i - 1]
    cns_str[index - i - 1] = cns_str[i] xor cns_str[index - i - 1]
    cns_str[i] = cns_str[i] xor cns_str[index - i - 1]
    eqv[i] = eqv[i] xor eqv[index - i - 1]
    eqv[index - i - 1] = eqv[i] xor eqv[index - i - 1]
    eqv[i] = eqv[i] xor eqv[index - i - 1]
    inc(i)
  cns_str[index] = 0
  ## #printf("%s\n", cns_str);
  when not defined(STATIC_ALLOCATE):
    i = 0
    while i < t_len:
      free_delta_group(msa_array[i])
      free(msa_array[i])
      inc(i)
    free(msa_array)
  else:
    clean_msa_working_space(msa_array, t_len + 1)
  free(coverage)
  free(local_nbase)
  return consensus

## #const unsigned int K = 8;

proc generate_consensus*(input_seq: cstringArray; n_seq: cuint; min_cov: cuint;
                        K: cuint; min_idt: cdouble): ptr consensus_data =
  var j: cuint
  var seq_count: cuint
  var aligned_seq_count: cuint
  var lk_ptr: ptr kmer_lookup
  var sa_ptr: seq_array
  var sda_ptr: seq_addr_array
  var kmer_match_ptr: ptr kmer_match
  var arange: ptr aln_range
  var aln: ptr alignment
  var tags_list: ptr ptr align_tags_t
  ## #char * consensus;
  var consensus: ptr consensus_data
  var max_diff: cdouble
  max_diff = 1.0 - min_idt
  seq_count = n_seq
  var xxx: pid_t = syscall(__NR_gettid)
  fprintf(stderr, "XX n_seq:%d pid:\x0A", n_seq, xxx)
  var
    start: timespec
    `end`: timespec
  clock_gettime(CLOCK_MONOTONIC_RAW, addr(start))
  var delta_us: uint64_t
  ## #for (j=0; j < seq_count; j++) {
  ## #    printf("seq_len: %u %u\n", j, strlen(input_seq[j]));
  ## #};
  ## #fflush(stdout);
  tags_list = calloc(seq_count, sizeof(ptr align_tags_t))
  lk_ptr = allocate_kmer_lookup(1 shl (K * 2))
  sa_ptr = allocate_seq(cast[seq_coor_t](strlen(input_seq[0])))
  sda_ptr = allocate_seq_addr(cast[seq_coor_t](strlen(input_seq[0])))
  add_sequence(0, K, input_seq[0], strlen(input_seq[0]), sda_ptr, sa_ptr, lk_ptr)
  ## #mask_k_mer(1 << (K * 2), lk_ptr, 16);
  mydelta_us = 0
  aligned_seq_count = 0
  j = 1
  while j < seq_count:
    ## #printf("seq_len: %ld %u\n", j, strlen(input_seq[j]));
    kmer_match_ptr = find_kmer_pos_for_seq(input_seq[j], strlen(input_seq[j]), K,
        sda_ptr, lk_ptr)
    const
      INDEL_ALLOWENCE_0 = 6
    arange = find_best_aln_range(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5)
    ## # narrow band to avoid aligning through big indels
    ## #printf("1:%ld %ld %ld %ld\n", arange_->s1, arange_->e1, arange_->s2, arange_->e2);
    ## #arange = find_best_aln_range2(kmer_match_ptr, K, K * INDEL_ALLOWENCE_0, 5);  // narrow band to avoid aligning through big indels
    ## #printf("2:%ld %ld %ld %ld\n\n", arange->s1, arange->e1, arange->s2, arange->e2);
    const
      INDEL_ALLOWENCE_1 = 0.1
    if arange.e1 - arange.s1 < 100 or arange.e2 - arange.s2 < 100 or
        abs((arange.e1 - arange.s1) - (arange.e2 - arange.s2)) >
        (int)(0.5 * INDEL_ALLOWENCE_1 *
        (arange.e1 - arange.s1 + arange.e2 - arange.s2)):
      free_kmer_match(kmer_match_ptr)
      free_aln_range(arange)
      continue
    const
      INDEL_ALLOWENCE_2 = 150
    aln = align(input_seq[j] + arange.s1, arange.e1 - arange.s1,
              input_seq[0] + arange.s2, arange.e2 - arange.s2, INDEL_ALLOWENCE_2, 1)
    if aln.aln_str_size > 500 and
        (cast[cdouble](aln.dist div cast[cdouble](aln.aln_str_size))) < max_diff:
      tags_list[aligned_seq_count] = get_align_tags(aln.q_aln_str, aln.t_aln_str,
          aln.aln_str_size, arange, j, 0)
      inc(aligned_seq_count)
    free_aln_range(arange)
    free_alignment(aln)
    free_kmer_match(kmer_match_ptr)
    inc(j)
  fprintf(stderr, "XX mydelta:%llfus for align()\x0A",
          cast[clongdouble]((mydelta_us) div 1000000.0))
  clock_gettime(CLOCK_MONOTONIC_RAW, addr(`end`))
  delta_us = (`end`.tv_sec - start.tv_sec) * 1000000 +
      (`end`.tv_nsec - start.tv_nsec) div 1000
  fprintf(stderr, "XX delta:%llfus\x0A",
          cast[clongdouble]((delta_us) div 1000000.0))
  if aligned_seq_count > 0:
    consensus = get_cns_from_align_tags(tags_list, aligned_seq_count,
                                      strlen(input_seq[0]), min_cov)
  else:
    ## # allocate an empty consensus sequence
    consensus = calloc(1, sizeof((consensus_data)))
    consensus.sequence = calloc(1, sizeof((char)))
    consensus.eqv = calloc(1, sizeof(cuint))
  ## #free(consensus);
  free_seq_addr_array(sda_ptr)
  free_seq_array(sa_ptr)
  free_kmer_lookup(lk_ptr)
  j = 0
  while j < aligned_seq_count:
    free_align_tags(tags_list[j])
    inc(j)
  free(tags_list)
  clock_gettime(CLOCK_MONOTONIC_RAW, addr(`end`))
  delta_us = (`end`.tv_sec - start.tv_sec) * 1000000 +
      (`end`.tv_nsec - start.tv_nsec) div 1000
  fprintf(stderr, "XX len(cons) %d %llfus\x0A", strlen(consensus.sequence),
          cast[clongdouble]((delta_us) div 1000000.0))
  return consensus

proc generate_utg_consensus*(input_seq: cstringArray; offset: ptr seq_coor_t;
                            n_seq: cuint; min_cov: cuint; K: cuint; min_idt: cdouble): ptr consensus_data =
  var j: cuint
  var seq_count: cuint
  var aligned_seq_count: cuint
  var arange: ptr aln_range
  var aln: ptr alignment
  var tags_list: ptr ptr align_tags_t
  ## #char * consensus;
  var consensus: ptr consensus_data
  var max_diff: cdouble
  var utg_len: seq_coor_t
  var r_len: seq_coor_t
  max_diff = 1.0 - min_idt
  seq_count = n_seq
  ## #**
  ## #    for (j=0; j < seq_count; j++) {
  ## #        printf("seq_len: %u %u\n", j, strlen(input_seq[j]));
  ## #    };
  ## #    fflush(stdout);
  ## #    *
  tags_list = calloc(seq_count + 1, sizeof(ptr align_tags_t))
  utg_len = strlen(input_seq[0])
  aligned_seq_count = 0
  arange = calloc(1, sizeof((aln_range)))
  arange.s1 = 0
  arange.e1 = strlen(input_seq[0])
  arange.s2 = 0
  arange.e2 = strlen(input_seq[0])
  tags_list[aligned_seq_count] = get_align_tags(input_seq[0], input_seq[0],
      strlen(input_seq[0]), arange, 0, 0)
  inc(aligned_seq_count, 1)
  j = 1
  while j < seq_count:
    arange.s1 = 0
    arange.e1 = strlen(input_seq[j]) - 1
    arange.s2 = 0
    arange.e2 = strlen(input_seq[j]) - 1
    r_len = strlen(input_seq[j])
    ## #printf("seq_len: %u %u\n", j, r_len);
    if offset[j] < 0:
      if (r_len + offset[j]) < 128:
        continue
      if r_len + offset[j] < utg_len:
        ## #printf("1: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = align(input_seq[j] - offset[j], r_len + offset[j], input_seq[0],
                  r_len + offset[j], 500, 1)
      else:
        ## #printf("2: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = align(input_seq[j] - offset[j], utg_len, input_seq[0], utg_len, 500, 1)
      offset[j] = 0
    else:
      if offset[j] > utg_len - 128:
        continue
      if offset[j] + r_len > utg_len:
        ## #printf("3: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = align(input_seq[j], utg_len - offset[j], input_seq[0] + offset[j],
                  utg_len - offset[j], 500, 1)
      else:
        ## #printf("4: %ld %u %u\n", offset[j], r_len, utg_len);
        aln = align(input_seq[j], r_len, input_seq[0] + offset[j], r_len, 500, 1)
    if aln.aln_str_size > 500 and
        (cast[cdouble](aln.dist div cast[cdouble](aln.aln_str_size))) < max_diff:
      tags_list[aligned_seq_count] = get_align_tags(aln.q_aln_str, aln.t_aln_str,
          aln.aln_str_size, arange, j, offset[j])
      inc(aligned_seq_count)
    free_alignment(aln)
    inc(j)
  free_aln_range(arange)
  if aligned_seq_count > 0:
    consensus = get_cns_from_align_tags(tags_list, aligned_seq_count, utg_len, 0)
  else:
    ## # allocate an empty consensus sequence
    consensus = calloc(1, sizeof((consensus_data)))
    consensus.sequence = calloc(1, sizeof((char)))
    consensus.eqv = calloc(1, sizeof(cuint))
  ## #free(consensus);
  j = 0
  while j < aligned_seq_count:
    free_align_tags(tags_list[j])
    inc(j)
  free(tags_list)
  return consensus

proc free_consensus_data*(consensus: ptr consensus_data) =
  free(consensus.sequence)
  free(consensus.eqv)
  free(consensus)
