This is a log of updates I made. This should be short and
  I will not log changes in graphs and stats.

# Todo

- Need to add auto install in makefile for ivar
- Need to fix infinite loop error in stich
  - I have not been able to replicate this, but it has
    poped up in benchmarking.

# Sep 20th, 2023

- Fixed a couple errors in stich that caused crashes.
- added ivar to the testing (00-scripts/runIvar.sh)
  - ivar needs to be in 00-programs

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
