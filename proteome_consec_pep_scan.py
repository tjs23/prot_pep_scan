import os, sys, io, subprocess, textwrap, time

try:
  from numba import njit, uint8, int64

except ModuleNotFoundError as err:

  print('* * Numba module not installed or accessible ** ')
  print('  Install via pip:\n     python3 -m pip install numba')
  print('  Install via conda:\n     conda install numba\n')
  raise(err)

import numpy as np

READ_BUFFER = 2**16

def open_file(file_path, mode=None, gzip_exts=('.gz','.gzip'), buffer_size=READ_BUFFER, partial=False):
  """
  GZIP and Python 2/3 agnostic file opening
  """
  import io
  
  if os.path.splitext(file_path)[1].lower() in gzip_exts:
    if mode and 'w' in mode:
      file_obj = io.BufferedWriter(gzip.open(file_path, mode), buffer_size)
      
    else:
      if partial:
        file_obj = io.BufferedReader(gzip.open(file_path, mode or 'rb'), buffer_size)
        
      else:
        try:
          file_obj = subprocess.Popen(['zcat', file_path], stdout=subprocess.PIPE).stdout
        except OSError:
          file_obj = io.BufferedReader(gzip.open(file_path, mode or 'rb'), buffer_size)
    
    if sys.version_info.major > 2:
      file_obj = io.TextIOWrapper(file_obj, encoding="utf-8")
 
  else:
    if sys.version_info.major > 2:
      file_obj = open(file_path, mode or 'rU', buffer_size, encoding='utf-8')
      
    else:
      file_obj = open(file_path, mode or 'rU', buffer_size)
  
  return file_obj


def iter_fasta(fasta_path):
  
  with open_file(fasta_path) as file_obj:
    name = None
    seq = []
 
    for line in file_obj:
      line = line.strip()
 
      if not line:
        continue
 
      if line[0] == '>':
        if name:
          yield name, ''.join(seq)

        seq  = []
        name = line[1:]
      else:
        seq.append(line)

    if name:
      yield name, ''.join(seq)


# Just-in-time compile with static typing and no Py objects; fast
@njit(int64[:](uint8[:], uint8[:,:], int64))
def scan_seq(seq, pep_array, max_gaps):

  n = len(seq)
  p, d = pep_array.shape
  hits = np.zeros(n, dtype=np.int64)
  out = np.empty(2, dtype=np.int64) # Max consecutive, position
  max_anonymous = 3
  n_anonymous = 0
  
  for i in range(n-d): # Start pos
   
    for j in range(p): # Possible peptides
      m = 0 # Matching residues
       
      for k in range(d): # Position in peptide
        if pep_array[j,k] == 88: # 'X' - any amino acid
          m += 1
        elif seq[i+k] == pep_array[j,k]:
          m += 1
        else:
          break  

      if m == d: # Whole peptide matched
        hits[i] = 1
        break
   
  max_consec = 0
  best_pos = 0
  gaps_remain = max_gaps
  
  for k in range(d): # Offset
    i = k
    start = 0
    n_consec = 0
    
    while i < n: # Scan whole length for any runs
      if hits[i] == 1:
        n_consec += 1
        if n_consec == 1:
          start = i
          
        if n_consec > max_consec: # Record best run
          max_consec = n_consec
          best_pos = start
      
      elif n_consec > 0 and gaps_remain > 0: # Gap: continue the run but don't end on a gap 
        n_consec += 1
        gaps_remain -= 1
      
      else:
        n_consec = 0
        start = 0
        gaps_remain = max_gaps
        
      i += d
  
  out[0] = max_consec
  out[1] = best_pos
       
  return out
        

def consec_pep_scan(fasta_paths, sub_peptides, req_consec=3, max_gaps=0, out_file=None, report_dt=0.2):
  
  start_t = time.time()
  rep_t = start_t
  n_hits = 0
  n_seqs = 0
  
  fmt1 = '\r    Inspected {:,} sequences found {:,} hits'.format
  fmt2 = '    Hit {}\n        \u2BA1 max_consec:{} num_diff:{} pos:{:4d} sub_seq:{}'.format
  fmt3 = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format
  
  if out_file:
    out_file_obj = open(out_file, 'w')
    write = out_file_obj.write
    write(fmt3('#seq_source','seq_name','max_consec','num_diff','start_pos','end_pos','sub_seq'))
    
  for fasta_path in fasta_paths:
  
    print(f'\nSearching {fasta_path}')
    n_peps = len(sub_peptides)
    p_size = len(sub_peptides[0])
    hits = []
    h = 0
    
    pep_array = np.fromstring(''.join(sub_peptides), dtype='uint8').reshape(n_peps, p_size)
 
    for i, (name, seq) in enumerate(iter_fasta(fasta_path)):
 
      seq_array = np.fromstring(seq, dtype='uint8')
 
      max_consec, best_pos = scan_seq(seq_array, pep_array, max_gaps)
      
      t = time.time()      
      if (t-rep_t) > report_dt or max_consec >= req_consec:
        msg = fmt1(i, h)
        sys.stdout.write(msg) # Must have enough columns to cover previous msg
        sys.stdout.flush()
        rep_t = t
        
      if max_consec >= req_consec:
        end_pos = best_pos+max_consec*p_size
        sub_seq = textwrap.wrap(seq[best_pos:end_pos], p_size)
        num_diff = len(set(sub_seq))
        sub_seq = ','.join(sub_seq)
        hits.append((sub_seq, max_consec, num_diff, best_pos, end_pos, name))
        h += 1  

    print(fmt1(i, h))
    n_hits += h
    n_seqs += i+1
        
    for sub_seq, max_consec, num_diff, start_pos, end_pos, name in hits:
      print(fmt2(name, max_consec, num_diff, start_pos, sub_seq))
            
      if out_file:
        write(fmt3(fasta_path, name, max_consec, num_diff, start_pos, end_pos, sub_seq))

  print(f'\nDone {n_seqs:,} sequences. Time taken:{time.time()-start_t:.2f} s')

  if out_file:
    print(f'Written {n_hits:,} hits to {out_file}')
    out_file_obj.close()
  
def main(argv=None):

  from argparse import ArgumentParser

  if argv is None:
    argv = sys.argv[1:]
  
  DEFAULT_N_CONSEC = 3
  DEFAULT_N_GAPS = 0
  
  epilog = 'For further help email tstevens@mrc-lmb.cam.ac.uk'

  arg_parse = ArgumentParser(prog='proteome_consec_pep_scan.py', description='Find sequences with the longest subsequences that match consecutive arrays of sub-peptide sequences',
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument(metavar='FASTA_FILE', nargs='+', dest='i',
                         help='One or more FASTA format sequence file paths (separated by spaces). Wildcards accepted.')

  arg_parse.add_argument('-p', metavar='PEP_SEQ', nargs='+', default=None,
                         help='Peptide amino acid sub-sequences to search for. May include "X" to match any amino acid. Space-separated, without quotes. For example: RKL PGG STQ')

  arg_parse.add_argument('-n', '--min-num-consec', default=DEFAULT_N_CONSEC, metavar='PEP_COUNT', type=int, dest="n",
                         help=f'Minimum number of consecutive sub-peptide sequences. Default {DEFAULT_N_CONSEC}')

  arg_parse.add_argument('-g', '--max-num-gaps', default=DEFAULT_N_GAPS, metavar='GAP_COUNT', type=int, dest="g",
                         help=f'Maximum number of unspecified, internal sub-peptides; gaps of sub-peptide length. Default {DEFAULT_N_GAPS}')

  arg_parse.add_argument('-o', '--out-file', default=None, metavar='OUT_TSV_FILE', dest="o",
                         help=f'Optional output file to write results as a tab-separated table')

  args = vars(arg_parse.parse_args(argv))

  fasta_paths = args['i']
  peps = args['p']
  n_consec = args['n']
  max_gaps = args['g']
  out_file = args['o']
  
  for fasta_path in fasta_paths:
    if not os.path.exists(fasta_path):
      print(f'ERROR: File {fasta_path} not found')
      sys.exit(1)     
  
  if not peps:
    print('ERROR: No peptide sub-sequences specified')
    arg_parse.print_help()
    sys.exit(1)
  
  size = len(peps[0])
  for i, pep in enumerate(peps):
    pep = pep.upper()
    peps[i] = pep
    
    if len(pep) != size:
      print('ERROR: Peptide sub-sequences not of the same length')
      sys.exit(1)
    
    if set(pep) == set(['X']):
      print('ERROR: Peptide sub-sequences cannot be all "X"s')
      sys.exit(1)
      
    
  consec_pep_scan(fasta_paths, peps, n_consec, max_gaps, out_file)


if __name__ == "__main__":
  main()



