#!/usr/bin/perl

  $src = @ARGV[0];
  $dest = @ARGV[1];
  $command = "cp $src/\*.include $dest";
  printf "$command \n";
  system "$command";
  $command = "cp $src/\*.py $dest";
  system "$command";
  $command = "ln -s $src/baseline $dest";
  system "$command";

