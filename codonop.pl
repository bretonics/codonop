#!/usr/bin/env perl
return 1 if caller(); # tests

use strict; use warnings; use diagnostics; use feature qw(say);
use Getopt::Long; use Pod::Usage;

use Bio::SeqIO; use Bio::Seq;

use FindBin; use lib "$FindBin::RealBin/lib";

# Own Modules (https://github.com/bretonics/Modules)
use MyConfig; use MyIO;
use Data::Dumper;
# ==============================================================================
#
#   CAPITAN:        Andres Breton, http://andresbreton.com
#   FILE:           codonop.pl
#   LICENSE:        GNU GPLv2
#   USAGE:
#   DEPENDENCIES:   - BioPerl
#                   - Modules
#
# ==============================================================================


#-------------------------------------------------------------------------------
# COMMAND LINE
my (@DNA, @RNA, $TABLE, $OUTFMT );

my $USAGE = "\n\n$0 [options]\n
Options:
  -dna          DNA sequence file(s)
  -rna          RNA sequence file(s)
  -table        Codon usage table
  -outfmt       Output format in DNA (default) or RNA [optional]
  -help         Shows this message
\n";

# OPTIONS
GetOptions(
  'dna=s{,}'      =>\@DNA,
  'rna=s{,}'      =>\@RNA,
  'table:s'       =>\$TABLE,
  'outfmt:s'      =>\$OUTFMT,
  # ''            =>\,
  help            =>sub{pod2usage($USAGE);}
)or pod2usage(2);

# checks(); # check CL arguments

#-------------------------------------------------------------------------------
# VARIABLES
my (%SEQS, %optimized);
$SEQS{'DNA'} = \@DNA;
$SEQS{'RNA'} = \@RNA;


#-------------------------------------------------------------------------------
# CALLS



foreach my $type (keys %SEQS) {
  my $codonUT = parseTable($TABLE);
  say "$type input sequence(s):";
# say join(",", keys $codonUT);

  my @SEQS = @{ $SEQS{$type} };
  foreach my $sequence (@SEQS){ # loop all sequence files
    my $inSeqIO   = Bio::SeqIO->new( -file => $sequence );

    while ( my $seqObj = $inSeqIO->next_seq) { # loop all sequences in file
      my $displayID = $seqObj->display_id;
      say "Optimizing $displayID";
      # my $protObj = $seq->translate;
      # my $prot = $protObj->seq;

      my $optSeq = optimize($codonUT, $seqObj, $OUTFMT);

  exit;
      $optimized{$displayID} = $optSeq;
    }
}
exit;



# say Dumper(\%optimized);
}

#-------------------------------------------------------------------------------
# SUBS
# sub checks {
#   unless ($TABLE) {
#     die "Need to provide both sequence and codon usage table files.";
#   }
#
# }

sub parseTable {
  my ($table) = @_;
  my $FH = getFH('<', $table);
  my %codons;

  while (my $line = <$FH>) {
    next if $. < 3; # skip headers
    last if($line =~ /^-/);
    $line =~ /^(\w{3})\s([\w\*])\s(\d\.\d+)\s+(\d+\.\d+)\s\((\d+)\)/;
    my $codon = $1;
    my $aa = $2;
    my $fraction = $3;
    my $frequency = $4;
    my $number = $5;

    $codons{$codon} = {
                    'aa'          => $aa,
                    'fraction'    => $fraction,
                    'frequency'   => $frequency,
                    'number'      => $number,
                  };
  }
  return(\%codons);
}

sub optimize {
  my ($codonUT, $seqObj, $OUTFMT) = @_;
  my @optSeq;
  my $seq     = uc $seqObj->seq;
  my $seqLen  = $seqObj->length;
  my $type    = $seqObj->alphabet;

  # Traverse through sequence &
  # Get best codon according to usage table
  for (my $i = 0; $i < $seqLen;) {
    my $codon = substr $seq, $i, 3;
# say "Using codon $codon ---";

    $codon = getCodon($codonUT, $codon, $type);
    if (!$OUTFMT) {
      $codon =~ s/U/T/g;
    } else {
      lc $OUTFMT eq 'dna' ? $codon =~ s/U/T/g : 0;
    }
    push @optSeq, $codon;
    $i = $i + 3;
  }

  return \@optSeq;
}

sub getCodon {
  my ($codonUT, $codon, $type) = @_;
  my %codonUT = %$codonUT;
  $codon =~ s/T/U/g if ($type eq 'dna');

  my @triplets = keys %codonUT;
  my $aa = $codonUT{$codon}{'aa'};


  my @AAcodons = grep { $codonUT{$_}{'aa'} eq $aa } @triplets;
# say "$codon -> $aa";
# say "\tFound: ", join(",", @AAcodons);
# print "\n\nFor AA '$aa' ($codon): ";
  my $best = 0;
  my ($fraction, $ultimate);
  foreach (@AAcodons) {
    $fraction = $codonUT{$_}{'fraction'};
# say "$codonUT{$codon}{'aa'} -> $codon";
# print "\nBest fraction is $best";
# print "\n\t'$_' --- ...found fraction: $fraction ";
    if($fraction > $best) {
      $best = $fraction;
# print "\tbest is now $best";
      $ultimate = $_; # best codon for specific AA
    }
  }
  # $optimized{$aa} = $ultimate;
# say "\nUltimate: $ultimate";
  return $ultimate;
}
