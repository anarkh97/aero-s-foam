#set terminal postscript color
set terminal postscript landscape color
set output "Discrepancies.ps"
set lmargin 10
set rmargin 5
set tmargin 5
set bmargin 5
set key bottom
set key box
plot "dsvm15/dsvm15_FETI_DP_displacy_2.dat" every ::2 using 1:2 title "test1", "baseline/dsvm15/dsvm15_FETI_DP_displacy_2.dat" every ::2 using 1:2 title "test2"
set lmargin 10
set rmargin 5
set tmargin 5
set bmargin 5
plot "dsvm15/dsvm15_FETI_DPH_displacy_2.dat" every ::2 using 1:2, "baseline/dsvm15/dsvm15_FETI_DPH_displacy_2.dat" every ::2 using 1:2
