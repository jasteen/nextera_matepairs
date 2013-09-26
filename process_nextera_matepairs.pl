#!/usr/bin/env perl

###############################################################################
#
# process_nextera_matepairs.pl
#
# This script will take raw data from an illumina matepair run and output reads 
# suitable for mapping or assembly
#
# Copyright (C) Jason Steen, 2013
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#modules and pragmas
use strict;
use warnings;
use String::Approx 'aindex';
use Getopt::Long;

##### Collect the options from the command line

#initialise default params
my $forward_read = 0;
my $reverse_read = 0;
my $prefix = "processed";
my $min_read_length = 50;
my $help = 0;

#actually get options
GetOptions ( 
	    'f|forward_read=s' => \$forward_read,
	    'r|reverse_read=s' => \$reverse_read,
	    'p|prefix=s' => \$prefix,
	    'm|min_read_length=i' => \$min_read_length,
	    'h|help' => \$help) or die "Invalid Parameters\n";

#give verbose usage on -h
exec("pod2usage $0 -v 2") if ($help);
#bail if either scaffold or contig file are unset
exec("pod2usage $0") if($forward_read eq 0 || $reverse_read eq 0);

#if all this passes, inform the user that we have started processing reads
print STDOUT "working...\n\n";

##### open input and output files

#open input files for reading. 
if ($forward_read =~ /.gz$/){
	open(R1, "gunzip -c $forward_read |") or die "cant open $_, are you sure it exists?";
	$forward_read =~ s/.gz//;
}else{
	open(R1, $forward_read) or die "cant open $_, are you sure it exists?";
}


if ($reverse_read =~ /.gz$/){
	open(R2, "gunzip -c $reverse_read |") or die "cant open $_, are you sure it exists?";
	$reverse_read =~ s/.gz//;
}else{
	open(R2, $reverse_read) or die "cant open $_, are you sure it exists?";
}

#deal with the path information on any input
if($forward_read =~ /\//){
	my @temp1 = split(/\//, $forward_read);
	my @temp2 = split(/\//, $reverse_read);
	$forward_read = pop(@temp1);
	$reverse_read = pop(@temp2);
}

#open our output files
open(OUTa, ">$prefix.$forward_read") or die "cant open $_ for writing, do you have write access to this directory?";
open(OUTb, ">$prefix.$reverse_read") or die "cant open $_ for writing, do you have write access to this directory?";
open(OUTc, ">$prefix.failed_filter_reads_interleaved.fq");

open(OUT2, ">$prefix.munted_reads.fq");
open(OUT3, ">$prefix.MatrixOfCounts.txt");

##### Main script

#define adaptors
my $dupJunctionAdaptor = 'CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG';
my $startDupJuncAdaptor10 = 'CTGTCTCTTA';
my $endDupJuncAdaptor10 = 'TAAGAGACAG';
my $singleJunctionAdaptor = 'CTGTCTCTTATACACATCT';
my $revSingleJunctionAdaptor = 'AGATGTGTATAAGAGACAG';

#declare global counters
my($CtrNoAdap, $CtrPairProc, $CtrFullHit, $CtrFiltOutLength, $CtrFwdSingle, $CtrRevSingle, $CtrWtf);


#required variables for the fastq reader
my @aux = undef;
my @aux2= undef;
my ($name, $seq, $qual, $name2, $seq2, $qual2);

#header for the matrix of values
print OUT3 join("\t", "R1_full", "R2_full", "R1_single", "R2_single", "R1_rev_single", "R2_rev_single", "R1_end", "R2_end", "R1_start", "R2_start"), "\n";

my($R1_full, $R2_full, $R1_single, $R2_single, $R1_rev_single, $R2_rev_single, $R1_end, $R2_end, $R1_start, $R2_start);
#Actual Guts of the script.  generate string matches and filter the reads based on them
while ((($name, $seq, $qual) = readfq(\*R1, \@aux))&&(($name2, $seq2, $qual2) = readfq(\*R2, \@aux2))) {
         $R1_full =             aindex($dupJunctionAdaptor, $seq);
         $R2_full =             aindex($dupJunctionAdaptor, $seq2);
         $R1_single =           aindex($singleJunctionAdaptor, $seq);
         $R2_single =           aindex($singleJunctionAdaptor, $seq2);
         $R1_rev_single =       aindex($revSingleJunctionAdaptor, $seq);
         $R2_rev_single =       aindex($revSingleJunctionAdaptor, $seq2);
         $R1_end =              aindex($startDupJuncAdaptor10, $seq);
         $R2_end =              aindex($startDupJuncAdaptor10, $seq2);
         $R1_start =            aindex($endDupJuncAdaptor10, $seq);
         $R2_start =            aindex($endDupJuncAdaptor10, $seq2);

        #Output matrix of aindex values
        print OUT3 join("\t", $R1_full, $R2_full, $R1_single, $R2_single, $R1_rev_single, $R2_rev_single, $R1_end, $R2_end, $R1_start, $R2_start), "\n";

 	#Start checking to see what we need to do with this pair
        $CtrPairProc++;
        if(($R1_full == -1)&&($R2_full == -1)&&($R1_single == -1)&&($R2_single == -1)&&($R1_rev_single == -1)&&($R2_rev_single == -1)&&($R1_end == -1)&&($R2_end == -1)&&($R1_start == -1)&&($R2_start)){
                #read-pair seems to contain no adaptors of any type
                $CtrNoAdap++;
                print OUTa "\@$name\n$seq\n+\n$qual\n";
                print OUTb "\@$name2\n$seq2\n+\n$qual2\n";
        }elsif(($R1_full != -1) || ($R2_full != -1)){
                $CtrFullHit++;
                if((($R1_full >= $min_read_length) || ($R1_full  == -1)) && (($R2_full >= $min_read_length) || ($R2_full  == -1))){
                        if($R1_full == -1){$R1_full = 250}; if ($R2_full == -1){$R2_full = 250};
                        my $R1_seq = substr($seq, 0, $R1_full);
                        my $R2_seq = substr($seq2, 0, $R2_full);
                        my $R1_qual = substr($qual, 0, $R1_full);
                        my $R2_qual = substr($qual2, 0, $R2_full);
                        #this could be included if you want interleaved pairs
                        #print OUT "\@$name\n$R1_seq\n+\n$R1_qual\n\@$name2\n$R2_seq\n+\n$R2_qual\n";
                        print OUTa "\@$name\n$R1_seq\n+\n$R1_qual\n";
                        print OUTb "\@$name2\n$R2_seq\n+\n$R2_qual\n";
                }else{
                        $CtrFiltOutLength++;
                        print OUTc "\@$name\n$seq\n+\n$qual\n\@$name2\n$seq2\n+\n$qual2\n";
                        next;
                }
        }elsif(($R1_single != -1) || ($R2_single!= -1)){
                $CtrFwdSingle++;
                #these reads probably contain forward single  adaptors.
                if($R1_single == -1){$R1_single = 250};if($R2_single == -1){$R2_single = 250};
                if((($R1_single >= $min_read_length) || ($R1_single == -1)) && (($R2_single >= $min_read_length) || ($R2_single == -1))){
                        my $R1_seq = substr($seq, 0, $R1_single);
                        my $R2_seq = substr($seq2, 0, $R2_single);
                        my $R1_qual = substr($qual, 0, $R1_single);
                        my $R2_qual = substr($qual2, 0, $R2_single);
                        print OUTa "\@$name\n$R1_seq\n+\n$R1_qual\n";
                        print OUTb "\@$name2\n$R2_seq\n+\n$R2_qual\n";
                }else{
                        $CtrFiltOutLength++;
                        print OUTc "\@$name\n$seq\n+\n$qual\n\@$name2\n$seq2\n+\n$qual2\n";
                        next;
                }
        }elsif(($R1_rev_single != -1) || ($R2_rev_single != -1)){
                $CtrRevSingle++;
                #these reads probably contain reverse single adaptors
        if($R1_rev_single == -1){$R1_rev_single = 250};if($R2_rev_single == -1){$R2_rev_single = 250};
                if((($R1_rev_single >= $min_read_length) || ($R1_rev_single == -1)) && (($R2_rev_single >= $min_read_length) || ($R2_rev_single == -1))){
                        my $R1_seq = substr($seq, 0, $R1_rev_single);
                        my $R2_seq = substr($seq2, 0, $R2_rev_single);
                        my $R1_qual = substr($qual, 0, $R1_rev_single);
                        my $R2_qual = substr($qual2, 0, $R2_rev_single);
                        print OUTa "\@$name\n$R1_seq\n+\n$R1_qual\n";
                        print OUTb "\@$name2\n$R2_seq\n+\n$R2_qual\n";
                }else{
                        $CtrFiltOutLength++;
                        print OUTc "\@$name\n$seq\n+\n$qual\n\@$name2\n$seq2\n+\n$qual2\n";
                        next;
                }
        }else{
                $CtrWtf++;
                print OUT2 "\@$name\n$seq\n+\n$qual\n\@$name2\n$seq2\n+\n$qual2\n";
        }
}

close(OUTa);
close(OUTb);
close(OUTc);
close(OUT2);
close(OUT3);
close(R1);
close(R2);

system("gzip $prefix.$forward_read");
system("gzip $prefix.$reverse_read");
system("$prefix.failed_filter_reads_interleaved.fq");
system("$prefix.munted_reads.fq.gz");

my $Remainder=$CtrPairProc-$CtrNoAdap-$CtrFullHit-$CtrFwdSingle-$CtrRevSingle-$CtrWtf;
print STDOUT join("\n",
		  "$CtrPairProc\tpairs were processed, of which:",
		  "$CtrNoAdap\thad no adaptors,",
		  "$CtrFullHit\thad a hit to the full dbl adaptor",
 		  "$CtrFwdSingle\thad a hit to the forward single adaptor",
		  "$CtrRevSingle\thas a hit to the reverse single adaptor",
		  "$CtrFiltOutLength\twere filtered because they didnt meet the length requirements of $min_read_length",
		  "$CtrWtf\tI had no idea what to do with","","","",
		  "The Folllowing files have been written:",
		  "$prefix.$forward_read.gz",
		  "$prefix.$reverse_read.gz",
		  "$prefix.failed_filter_reads_interleaved.fq.gz",
		  "$prefix.munted_reads.fq.gz",
  		  "$prefix.MatrixOfCounts.txt","","",
		 );
		
exit(0);


####heng_li fastq reading subroutine
sub readfq {
        my ($fh, $aux) = @_;
        @$aux = [undef, 0] if (!defined(@$aux));
        return if ($aux->[1]);
        if (!defined($aux->[0])) {
                while (<$fh>) {
                        chomp;
                        if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                                $aux->[0] = $_;
                                last;
                        }
                }
                if (!defined($aux->[0])) {
                        $aux->[1] = 1;
                        return;
                }
        }
        my $name = /^.(\S+)/? $1 : '';
        my $seq = '';
        my $c;
        $aux->[0] = undef;
        while (<$fh>) {
                chomp;
                $c = substr($_, 0, 1);
                last if ($c eq '>' || $c eq '@' || $c eq '+');
                $seq .= $_;
        }
        $aux->[0] = $_;
        $aux->[1] = 1 if (!defined($aux->[0]));
        return ($name, $seq) if ($c ne '+');
        my $qual = '';
        while (<$fh>) {
                chomp;
                $qual .= $_;
                if (length($qual) >= length($seq)) {
                        $aux->[0] = undef;
                        return ($name, $seq, $qual);
                }
        }
        $aux->[1] = 1;
        return ($name, $seq);
}


__DATA__

=head1 NAME

Process_nextera_matepairs.pl

=head1 COPYRIGHT

copyright (C) Jason Steen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

this script will take a raw illumina matepair run, and process each pair to remove internal adaptors.  Processed reads output by this script suitable for mapping.

=head1 SYNOPSIS

Process_nextera_matepairs [-help|h] -f -r [-m -p]
	
	-f, -fwd_read		File containing the first read of the pair (can be gzipped)
	-r, -rev_read		File containing the second read of the pair (can be gzipped)
	-p, -prefix		Prefix for the output files (optional, default = processed)
	-m, -min_read-length	Minimum read length of each pair to output reads (default = 50)
	[-h, -help] 		Displays this usage information



=cut

