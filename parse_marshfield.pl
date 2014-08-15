#!/usr/bin/perl

use strict;
use warnings;
use LWP;
use Carp;
use List::Util qw(shuffle min max);
use Readonly;
use Data::Dumper;
use Storable;
use Getopt::Long;
use File::Basename;

###################### Headers #####################



####################################################
main();

sub main {
    my $r_arg = arg_user();
    my $r_annotation = read_annotation($r_arg);
    generate_list_of_files($r_arg, $r_annotation);
    
};

sub arg_user {
    my %args;

    my $path_to_files;
    my $file;
    my $annotation;
    my $debug;
    my $help;


    if ($help) {
        help();
    }


    ### load parameters:
    GetOptions (
                'path_to_files=s'           => \$path_to_files,
                'file=s'                    => \$file,
                'annotation=s'              => \$annotation,
                'debug'                     => \$debug,
                'help'                      => \$help,
                );


    if ($help) {
        help();
    }

    if (not defined $path_to_files or
        not defined $annotation or
        0 ==1) {
        help();
    }

    if (not -d $path_to_files) {
        die "path to file: $path_to_files does not exist\n";
    }


    $args{'annotation'}         = $annotation;
    $args{'path_to_files'}      = $path_to_files;
    $args{'file'}               = $file;
    $args{'debug'}              = $debug;



    return \%args;
    
}


sub help {
    my $name = basename ($0);
    my $spad = (' ') x length $name;

    print <<EOF;
usage: 
$name --path_to_files <path> --annotation <file> [--file <file>] [--debug] [--help]
$spad Output is a bed like format with chr,start,stop at the begining.
Example:
perl $name  --path_to_files ./marshfield/par  --annotation Marshfield_annotations.txt

$spad Mandatory parameters:
$spad --path_to_files\tlocations of markers to parse
$spad --annotation\tfile with annotations of markers (Tom's folder)


$spad Optional parameters:
$spad --file\tThis is for debugging purposes. The program will only analze this file (no path is needed).
$spad --debug\tExtra output messages
$spad --help\tThis help messages


EOF
    exit(64);

}


sub read_annotation {
    my ($r_arg) = @_;

    my %hash;
    my $file_annotation = $r_arg->{annotation};
    my $FH;
    open $FH, '<', $file_annotation or
        confess "can't open $file_annotation :$!";

    while (<$FH>) {
        chomp;
        my ($chrom,
            $start,
            $end,
            $locus,
            $marker1,
            $marker2,
            $strand,
            $annotated,
            $product,
            $size) = split /\s+/;


        $hash{$marker2}{'location'} = "$chrom\t$start\t$end";
        $hash{$marker2}{'product'} = $product;
    }   
    close $FH;

    return \%hash;
}

sub generate_list_of_files {
    my ($r_arg, $r_annotation) = @_;

    my $path_to_files = $r_arg->{path_to_files};
    my $file_filter = $r_arg->{file}; #this is an optional parameter. if it is defined the program will only analyze this file
    opendir(DIR, $path_to_files) or die $!;
    print_header();

    while (my $file = readdir(DIR)) {

        next unless (-f "$path_to_files/$file");
        next if ($file =~m/\.zip/);
        next if (defined $file_filter and $file_filter ne $file);

        if ($r_arg->{debug}) {
            print "$path_to_files/$file\n";
        }

        my $marker = $file;
        $marker=~s/\..*//;
        $marker=uc($marker);

        next if (not exists $r_annotation->{$marker});
        analyze_file($r_arg, "$path_to_files/$file", $r_annotation->{$marker});

    }
    
}


sub analyze_file {
    my ($r_arg, $file, $r_marker) = @_;

    my $debug = $r_arg->{debug};
    my $FH;
    open $FH, '<', $file or
        confess "can't open $file: $!";


    my %dist;
   
    my $flag_freq = 0; #indicates if we are in the frequency region


    LINE:
    while (<$FH>) {
        
        chomp;
        s/^\s+//;

        if (m/Frequencies/) {#start of freq section
            $flag_freq = 1;
            next LINE;
        }
        if ($flag_freq == 1 and $_ eq '') {#end of freq section
            $flag_freq = 2;
            next LINE;
        }

        if ($flag_freq == 1 and m/\)/) {
            m/^(\d+)\)\s+(\d+)\s+-\s+(.*)/;
            my $number = $1;
            my $allele = $2;
            my $freq = $3;
            $freq =~ s/\x0D//;

            $dist{data}{$allele} = $freq;
            $dist{number} = $number;
            $dist{total} += $freq;
        }

        if (m/Het: (.+)/) {
            $dist{het} = $1;
            $dist{het} =~ s/\x0D//;
        }
        if (m/PIC: (.+)/) {
            $dist{pic} = $1;
            $dist{pic} =~ s/\x0D//;
        }
        if (m/Individuals Used:\s+(\d+)/i) {
            $dist{ind} = $1;
            $dist{ind} =~ s/\x0D//;
        }
        if (m/Marker:\s+(.+)/) {
            $dist{marker} = $1;
            $dist{marker} =~ s/\x0D//;
        }

    }
    close $FH;

    if (sanity_test($file, \%dist, $debug, $r_marker)) {
        print_hash($r_arg, \%dist, $file, $r_marker);
    }
    
}

sub sanity_test {
    my ($file, $r_hash, $debug, $r_marker) = @_;

    #test all variables are defined:
    foreach my $var ('ind', 'het', 'pic', 'marker') {
        if (not defined $r_hash->{$var}) {
            print "$var is not defined in $file\n" if ($debug);
            return 0;
        }
    }
    
    if (not defined $r_hash->{total} or $r_hash->{total}<0.99) {
        print "freq is not defined or not equal to 1 in $file\n" if ($debug);
        return 0;
    }

    #fill alleles with 0 freq:
    my $min = min keys %{$r_hash->{data}};
    my $max = max keys %{$r_hash->{data}};

    my $expected_product = $r_marker->{product}; #this is what the reference should be!
    my $flag_expect_product_seen = 0;
    foreach my $allele ($min..$max) {
        if (not defined $r_hash->{data}->{$allele}) {
            $r_hash->{data}->{$allele} = 0;
            $r_hash->{number}++;
        }
        if ($allele eq $expected_product) {
            $flag_expect_product_seen = 1;
        }
    }

    if ($flag_expect_product_seen == 0) {
        print "Expected product $expected_product has not seen in $file. Min: $min and Max: $max. \n" if ($debug);
        return 0;
    }

    return 1;
}

sub print_hash {
    my ($r_arg, $r_hash, $file, $r_marker) = @_;

    print $r_marker->{location},"\t";
    print "$file\t";
   

    print   $r_hash->{marker},"\t",
            $r_hash->{het},"\t",
            $r_hash->{pic},"\t",
            $r_hash->{ind},"\t";
            

     my $min = min keys %{$r_hash->{data}};
     my $max = max keys %{$r_hash->{data}};
     print "$min\t";
     print "$max\t";
     print $r_marker->{product},"\t";
     print $r_hash->{number},"\t";
 
     foreach my $allele ($min..$max) {
         print $allele,":",$r_hash->{data}->{$allele},"\t";
     }
 
 

   print "\n";
}

sub print_header {
    print "chr\t","start\t","stop\t","file\t","marker\t","het\t","pic\t","#individuals\t","min_allele\t","max_allele\t", "product\t","number_of_alleles\t","(allele:prop)\n";

}
