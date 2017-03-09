## #
## #  =====================================================================================
## # 
## #        Filename:  kmer_count.c
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

var KMERMATCHINC*: cuint = 10000

proc compare_seq_coor*(a: pointer; b: pointer): cint =
  var arg1: ptr seq_coor_t = a
  var arg2: ptr seq_coor_t = b
  return (arg1[]) - (arg2[])

proc allocate_kmer_lookup*(size: seq_coor_t): ptr kmer_lookup =
  var kl: ptr kmer_lookup
  ## #printf("%lu is allocated for kmer lookup\n", size);
  kl = cast[ptr kmer_lookup](malloc(size * sizeof((kmer_lookup))))
  init_kmer_lookup(kl, size)
  return kl

proc init_kmer_lookup*(kl: ptr kmer_lookup; size: seq_coor_t) =
  var i: seq_coor_t
  ## #printf("%lu is allocated for kmer lookup\n", size);
  i = 0
  while i < size:
    kl[i].start = INT_MAX
    kl[i].last = INT_MAX
    kl[i].count = 0
    inc(i)

proc free_kmer_lookup*(`ptr`: ptr kmer_lookup) =
  free(`ptr`)

proc allocate_seq*(size: seq_coor_t): seq_array =
  var sa: seq_array
  sa = cast[seq_array](malloc(size * sizeof((base))))
  init_seq_array(sa, size)
  return sa

proc init_seq_array*(sa: seq_array; size: seq_coor_t) =
  var i: seq_coor_t
  i = 0
  while i < size:
    sa[i] = 0x000000FF
    inc(i)

proc free_seq_array*(sa: seq_array) =
  free(sa)

proc allocate_seq_addr*(size: seq_coor_t): seq_addr_array =
  return cast[seq_addr_array](calloc(size, sizeof((seq_addr))))

proc free_seq_addr_array*(sda: seq_addr_array) =
  free(sda)

proc get_kmer_bitvector*(sa: seq_array; K: cuint): seq_coor_t =
  var i: cuint
  var kmer_bv: seq_coor_t = 0
  var kmer_mask: seq_coor_t
  kmer_mask = 0
  i = 0
  while i < K:
    kmer_mask = kmer_mask shl 2
    kmer_mask = kmer_mask or 0x00000003
    inc(i)
  i = 0
  while i < K:
    kmer_bv = kmer_bv shl 2
    kmer_bv = kmer_bv or ((cast[cuint](sa[i])) and 0x00000003)
    inc(i)
  return kmer_bv

proc add_sequence*(start: seq_coor_t; K: cuint; seq: cstring; seq_len: seq_coor_t;
                  sda: seq_addr_array; sa: seq_array; lk: ptr kmer_lookup) =
  var i: seq_coor_t
  var kmer_bv: seq_coor_t
  var kmer_mask: seq_coor_t
  kmer_mask = 0
  i = 0
  while i < K:
    kmer_mask = kmer_mask shl 2
    kmer_mask = kmer_mask or 0x00000003
    inc(i)
  i = 0
  while i < seq_len:
    case seq[i]
    of 'A':
      sa[start + i] = 0
    of 'C':
      sa[start + i] = 1
    of 'G':
      sa[start + i] = 2
    of 'T':
      sa[start + i] = 3
    inc(i)
  kmer_bv = get_kmer_bitvector(sa + start, K)
  i = 0
  while i < seq_len - K:
    ## #printf("%lu %lu\n", i, kmer_bv);
    ## #printf("lk before init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
    if lk[kmer_bv].start == INT_MAX:
      lk[kmer_bv].start = start + i
      lk[kmer_bv].last = start + i
      inc(lk[kmer_bv].count, 1)
      ## #printf("lk init: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
    else:
      sda[lk[kmer_bv].last] = start + i
      inc(lk[kmer_bv].count, 1)
      lk[kmer_bv].last = start + i
      ## #printf("lk change: %lu %lu %lu\n", kmer_bv, lk[kmer_bv].start, lk[kmer_bv].last);
    kmer_bv = kmer_bv shl 2
    kmer_bv = kmer_bv or sa[start + i + K]
    kmer_bv = kmer_bv and kmer_mask
    inc(i)

proc mask_k_mer*(size: seq_coor_t; kl: ptr kmer_lookup; threshold: seq_coor_t) =
  var i: seq_coor_t
  i = 0
  while i < size:
    if kl[i].count > threshold:
      kl[i].start = INT_MAX
      kl[i].last = INT_MAX
      ## #kl[i].count = 0;
    inc(i)

proc find_kmer_pos_for_seq*(seq: cstring; seq_len: seq_coor_t; K: cuint;
                           sda: seq_addr_array; lk: ptr kmer_lookup): ptr kmer_match =
  var i: seq_coor_t
  var kmer_bv: seq_coor_t
  var kmer_mask: seq_coor_t
  var kmer_pos: seq_coor_t
  var next_kmer_pos: seq_coor_t
  var half_K: cuint
  var kmer_match_rtn_allocation_size: seq_coor_t = KMERMATCHINC
  var kmer_match_rtn: ptr kmer_match
  var sa: ptr base
  kmer_match_rtn = cast[ptr kmer_match](malloc(sizeof((kmer_match))))
  kmer_match_rtn.count = 0
  kmer_match_rtn.query_pos = cast[ptr seq_coor_t](calloc(
      kmer_match_rtn_allocation_size, sizeof((seq_coor_t))))
  kmer_match_rtn.target_pos = cast[ptr seq_coor_t](calloc(
      kmer_match_rtn_allocation_size, sizeof((seq_coor_t))))
  sa = calloc(seq_len, sizeof((base)))
  kmer_mask = 0
  i = 0
  while i < K:
    kmer_mask = kmer_mask shl 2
    kmer_mask = kmer_mask or 0x00000003
    inc(i)
  i = 0
  while i < seq_len:
    case seq[i]
    of 'A':
      sa[i] = 0
    of 'C':
      sa[i] = 1
    of 'G':
      sa[i] = 2
    of 'T':
      sa[i] = 3
    inc(i)
  kmer_bv = get_kmer_bitvector(sa, K)
  half_K = K shr 1
  i = 0
  while i < seq_len - K:
    kmer_bv = get_kmer_bitvector(sa + i, K)
    if lk[kmer_bv].start == INT_MAX:
      ## #for high count k-mers
      continue
    kmer_pos = lk[kmer_bv].start
    next_kmer_pos = sda[kmer_pos]
    kmer_match_rtn.query_pos[kmer_match_rtn.count] = i
    kmer_match_rtn.target_pos[kmer_match_rtn.count] = kmer_pos
    inc(kmer_match_rtn.count, 1)
    if kmer_match_rtn.count > kmer_match_rtn_allocation_size - 1000:
      inc(kmer_match_rtn_allocation_size, KMERMATCHINC)
      kmer_match_rtn.query_pos = cast[ptr seq_coor_t](realloc(
          kmer_match_rtn.query_pos,
          kmer_match_rtn_allocation_size * sizeof((seq_coor_t))))
      kmer_match_rtn.target_pos = cast[ptr seq_coor_t](realloc(
          kmer_match_rtn.target_pos,
          kmer_match_rtn_allocation_size * sizeof((seq_coor_t))))
    while next_kmer_pos > kmer_pos:
      kmer_pos = next_kmer_pos
      next_kmer_pos = sda[kmer_pos]
      kmer_match_rtn.query_pos[kmer_match_rtn.count] = i
      kmer_match_rtn.target_pos[kmer_match_rtn.count] = kmer_pos
      inc(kmer_match_rtn.count, 1)
      if kmer_match_rtn.count > kmer_match_rtn_allocation_size - 1000:
        inc(kmer_match_rtn_allocation_size, KMERMATCHINC)
        kmer_match_rtn.query_pos = cast[ptr seq_coor_t](realloc(
            kmer_match_rtn.query_pos,
            kmer_match_rtn_allocation_size * sizeof((seq_coor_t))))
        kmer_match_rtn.target_pos = cast[ptr seq_coor_t](realloc(
            kmer_match_rtn.target_pos,
            kmer_match_rtn_allocation_size * sizeof((seq_coor_t))))
    inc(i, half_K)
  free(sa)
  return kmer_match_rtn

proc free_kmer_match*(`ptr`: ptr kmer_match) =
  free(`ptr`.query_pos)
  free(`ptr`.target_pos)
  free(`ptr`)

proc find_best_aln_range*(km_ptr: ptr kmer_match; K: seq_coor_t; bin_size: seq_coor_t;
                         count_th: seq_coor_t): ptr aln_range =
  var i: seq_coor_t
  var j: seq_coor_t
  var
    q_min: seq_coor_t
    q_max: seq_coor_t
    t_min: seq_coor_t
    t_max: seq_coor_t
  var d_count: ptr seq_coor_t
  var q_coor: ptr seq_coor_t
  var t_coor: ptr seq_coor_t
  var arange: ptr aln_range
  var
    d: clong
    d_min: clong
    d_max: clong
  var cur_score: clong
  var max_score: clong
  var max_k_mer_count: clong
  var max_k_mer_bin: clong
  var cur_start: seq_coor_t
  arange = calloc(1, sizeof((aln_range)))
  q_min = INT_MAX
  q_max = 0
  t_min = INT_MAX
  t_max = 0
  d_min = INT_MAX
  d_max = LONG_MIN
  i = 0
  while i < km_ptr.count:
    if km_ptr.query_pos[i] < q_min:
      q_min = km_ptr.query_pos[i]
    if km_ptr.query_pos[i] > q_max:
      q_max = km_ptr.query_pos[i]
    if km_ptr.target_pos[i] < t_min:
      t_min = km_ptr.target_pos[i]
    if km_ptr.query_pos[i] > t_max:
      t_max = km_ptr.target_pos[i]
    d = cast[clong](km_ptr.query_pos[i]) - cast[clong](km_ptr.target_pos[i])
    if d < d_min:
      d_min = d
    if d > d_max:
      d_max = d
    inc(i)
  ## #printf("%lu %ld %ld\n" , km_ptr->count, d_min, d_max);
  d_count = calloc((d_max - d_min) div bin_size + 1, sizeof((seq_coor_t)))
  q_coor = calloc(km_ptr.count, sizeof((seq_coor_t)))
  t_coor = calloc(km_ptr.count, sizeof((seq_coor_t)))
  i = 0
  while i < km_ptr.count:
    d = cast[clong]((km_ptr.query_pos[i])) - cast[clong]((km_ptr.target_pos[i]))
    inc(d_count[(d - d_min) div cast[clong](bin_size)], 1)
    q_coor[i] = INT_MAX
    t_coor[i] = INT_MAX
    inc(i)
  j = 0
  max_k_mer_count = 0
  max_k_mer_bin = INT_MAX
  i = 0
  while i < km_ptr.count:
    d = cast[clong]((km_ptr.query_pos[i])) - cast[clong]((km_ptr.target_pos[i]))
    if d_count[(d - d_min) div cast[clong](bin_size)] > max_k_mer_count:
      max_k_mer_count = d_count[(d - d_min) div cast[clong](bin_size)]
      max_k_mer_bin = (d - d_min) div cast[clong](bin_size)
    inc(i)
  ## #printf("k_mer: %lu %lu\n" , max_k_mer_count, max_k_mer_bin);
  if max_k_mer_bin != INT_MAX and max_k_mer_count > count_th:
    i = 0
    while i < km_ptr.count:
      d = cast[clong]((km_ptr.query_pos[i])) - cast[clong]((km_ptr.target_pos[i]))
      if labs(((d - d_min) div cast[clong](bin_size)) - max_k_mer_bin) > 5:
        continue
      if d_count[(d - d_min) div cast[clong](bin_size)] > count_th:
        q_coor[j] = km_ptr.query_pos[i]
        t_coor[j] = km_ptr.target_pos[i]
        ## #printf("d_count: %lu %lu\n" ,i, d_count[(d - d_min)/ (long int) bin_size]);
        ## #printf("coor: %lu %lu\n" , q_coor[j], t_coor[j]);
        inc(j)
      inc(i)
  if j > 1:
    arange.s1 = q_coor[0]
    arange.e1 = q_coor[0]
    arange.s2 = t_coor[0]
    arange.e2 = t_coor[0]
    arange.score = 0
    max_score = 0
    cur_score = 0
    cur_start = 0
    i = 1
    while i < j:
      inc(cur_score, 32 - (q_coor[i] - q_coor[i - 1]))
      ## #printf("deltaD, %lu %ld\n", q_coor[i] - q_coor[i-1], cur_score);
      if cur_score < 0:
        cur_score = 0
        cur_start = i
      elif cur_score > max_score:
        arange.s1 = q_coor[cur_start]
        arange.s2 = t_coor[cur_start]
        arange.e1 = q_coor[i]
        arange.e2 = t_coor[i]
        max_score = cur_score
        arange.score = max_score
        ## #printf("%lu %lu %lu %lu\n", arange.s1, arange.e1, arange.s2, arange.e2);
      inc(i)
  else:
    arange.s1 = 0
    arange.e1 = 0
    arange.s2 = 0
    arange.e2 = 0
    arange.score = 0
  ## # printf("free\n");
  free(d_count)
  free(q_coor)
  free(t_coor)
  return arange

proc find_best_aln_range2*(km_ptr: ptr kmer_match; K: seq_coor_t;
                          bin_width: seq_coor_t; count_th: seq_coor_t): ptr aln_range =
  var d_coor: ptr seq_coor_t
  var hit_score: ptr seq_coor_t
  var hit_count: ptr seq_coor_t
  var last_hit: ptr seq_coor_t
  var
    max_q: seq_coor_t
    max_t: seq_coor_t
  var
    s: seq_coor_t
    e: seq_coor_t
    max_s: seq_coor_t
    max_e: seq_coor_t
    max_span: seq_coor_t
    d_s: seq_coor_t
    d_e: seq_coor_t
    delta: seq_coor_t
    d_len: seq_coor_t
  var
    px: seq_coor_t
    py: seq_coor_t
    cx: seq_coor_t
    cy: seq_coor_t
  var max_hit_idx: seq_coor_t
  var
    max_hit_score: seq_coor_t
    max_hit_count: seq_coor_t
  var
    i: seq_coor_t
    j: seq_coor_t
  var
    candidate_idx: seq_coor_t
    max_d: seq_coor_t
    d: seq_coor_t
  var arange: ptr aln_range
  arange = calloc(1, sizeof((aln_range)))
  d_coor = calloc(km_ptr.count, sizeof((seq_coor_t)))
  max_q = - 1
  max_t = - 1
  i = 0
  while i < km_ptr.count:
    d_coor[i] = km_ptr.query_pos[i] - km_ptr.target_pos[i]
    max_q = if max_q > km_ptr.query_pos[i]: max_q else: km_ptr.query_pos[i]
    max_t = if max_t > km_ptr.target_pos[i]: max_q else: km_ptr.target_pos[i]
    inc(i)
  qsort(d_coor, km_ptr.count, sizeof((seq_coor_t)), compare_seq_coor)
  s = 0
  e = 0
  max_s = - 1
  max_e = - 1
  max_span = - 1
  delta = cast[clong]((0.05 * (max_q + max_t)))
  d_len = km_ptr.count
  d_s = - 1
  d_e = - 1
  while 1:
    d_s = d_coor[s]
    d_e = d_coor[e]
    while d_e < d_s + delta and e < d_len - 1:
      inc(e, 1)
      d_e = d_coor[e]
    if max_span == - 1 or e - s > max_span:
      max_span = e - s
      max_s = s
      max_e = e
    inc(s, 1)
    if s == d_len or e == d_len:
      break
  if max_s == - 1 or max_e == - 1 or max_e - max_s < 32:
    arange.s1 = 0
    arange.e1 = 0
    arange.s2 = 0
    arange.e2 = 0
    arange.score = 0
    free(d_coor)
    return arange
  last_hit = calloc(km_ptr.count, sizeof((seq_coor_t)))
  hit_score = calloc(km_ptr.count, sizeof((seq_coor_t)))
  hit_count = calloc(km_ptr.count, sizeof((seq_coor_t)))
  i = 0
  while i < km_ptr.count:
    last_hit[i] = - 1
    hit_score[i] = 0
    hit_count[i] = 0
    inc(i)
  max_hit_idx = - 1
  max_hit_score = 0
  i = 0
  while i < km_ptr.count:
    cx = km_ptr.query_pos[i]
    cy = km_ptr.target_pos[i]
    d = cx - cy
    if d < d_coor[max_s] or d > d_coor[max_e]: continue
    j = i - 1
    candidate_idx = - 1
    max_d = 65535
    while 1:
      if j < 0: break
      px = km_ptr.query_pos[j]
      py = km_ptr.target_pos[j]
      d = px - py
      if d < d_coor[max_s] or d > d_coor[max_e]:
        dec(j)
        continue
      if cx - px > 320:
        break
      if cy > py and cx - px + cy - py < max_d and cy - py <= 320:
        max_d = cx - px + cy - py
        candidate_idx = j
      dec(j)
    if candidate_idx != - 1:
      last_hit[i] = candidate_idx
      hit_score[i] = hit_score[candidate_idx] + (64 - max_d)
      hit_count[i] = hit_count[candidate_idx] + 1
      if hit_score[i] < 0:
        hit_score[i] = 0
        hit_count[i] = 0
    else:
      hit_score[i] = 0
      hit_count[i] = 0
    if hit_score[i] > max_hit_score:
      max_hit_score = hit_score[i]
      max_hit_count = hit_count[i]
      max_hit_idx = i
    inc(i)
  if max_hit_idx == - 1:
    arange.s1 = 0
    arange.e1 = 0
    arange.s2 = 0
    arange.e2 = 0
    arange.score = 0
    free(d_coor)
    free(last_hit)
    free(hit_score)
    free(hit_count)
    return arange
  arange.score = max_hit_count + 1
  arange.e1 = km_ptr.query_pos[max_hit_idx]
  arange.e2 = km_ptr.target_pos[max_hit_idx]
  i = max_hit_idx
  while last_hit[i] != - 1:
    i = last_hit[i]
  arange.s1 = km_ptr.query_pos[i]
  arange.s2 = km_ptr.target_pos[i]
  free(d_coor)
  free(last_hit)
  free(hit_score)
  free(hit_count)
  return arange

proc free_aln_range*(arange: ptr aln_range) =
  free(arange)
