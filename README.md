# rga
Scritps used for the study of RGAs

Command lines to run each of the used softwares:

1) Phobius
$ perl phobius.pl -short file.fas > file.fas.phobius 

2) TMHMM2
$ perl tmhmm file.fasta > file.fas.tmhmm

Note: Line 33 from the file "tmhmm" was commented, and full path was added to model on line 56.

3) Coils2
$ export COILSDIR=/path/to/coils
$ ./ncoils-osf -min_seg 1 -f < file.fas > file.fas.coils

4) PfamScan
perl pfam_scan.pl -fasta file.fas -dir ./PfamDATA/ -outfile file.fas.pfamscan

5) InterproScan
./interproscan.sh -i file.fas -f tsv
