#!/usr/bin/perl

  $hostfile = $ENV{'PBS_NODEFILE'};
  print "hostfile = $hostfile\n";
  system("rm host.*");
  @machs = split(/\n/,`cat $hostfile`);
  $len = @machs;
  $j = 0;
  for($i = 0; $i < $len; $i++) {
    if(($i % 4) == 0) {
      $filename = ">host.".$j;
      open(OFILE,$filename);
      $j++;
    }
    print OFILE "@machs[$i] \n";
    if(($i % 4) == 3) { close OFILE; }
    
  }
