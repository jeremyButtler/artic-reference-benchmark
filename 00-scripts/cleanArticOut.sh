#!/usr/bin/bash

# This script removes all files, except the consensus
# from the artic output
articDirStr="$1";

rm -r "$articDirStr.a"*;
rm -r "$articDirStr.1"*;
rm -r "$articDirStr.2"*;
rm -r "$articDirStr.cov"*;
rm -r "$articDirStr.fail"*;
rm -r "$articDirStr.m"*;
rm -r "$articDirStr.p"*;
rm -r "$articDirStr.s"*;
rm -r "$articDirStr.t"*;

exit;
