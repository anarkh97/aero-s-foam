#!/usr/bin/perl

   @lst = split("\n",`cat lst`);
   foreach (@lst) {
     @words = split(/\s+/,$_);
     $filename = @words[1];
     $pathhead = "~/FEM_test2/";
     $file2 = $pathhead.$filename;
 #   print "$filename $file2\n";
     $command = "diff $filename $file2\n";
     $ret = `$command`;
     if($ret != "") {print "$command\n$ret\n"; }
   }
