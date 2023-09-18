This is a log of updates I made. This should be short and
  I will not log changes in graphs and stats.

# Sep 18th, 2023

- Updated the benchAll.sh script to not run tests that have
  already been run. It will still subsample and mutate
  genomes that have been run, however, this is a small
  amount of the run time.
  - For a restart You will have to clear everything, but
    the prefix-stats.tsv in 09.
  - You will also have to move backup.xxx.reference.fasta
    to its correct directory (this is the unmutated 
    reference).
- Ended up with some mutated consensus in 06-schemes
  - This was not a problem with the 07 directory (current
    graphed results).
