#!/usr/bin/perl

  $src = @ARGV[0];
  $dest = @ARGV[1];
  if($src != $dest) {
    $command = "cp $src/\*.include $dest";
    system "$command";
    $command = "cp $src/\*.py $dest";
    system "$command";
    $command = "ln -s $src/baseline $dest";
    system "$command";
  }  

