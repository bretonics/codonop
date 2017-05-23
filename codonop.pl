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
#   USAGE:          Optimize codons according to codon usage
#   DEPENDENCIES:   - BioPerl
#                   - Modules
#
# ==============================================================================


#-------------------------------------------------------------------------------
# COMMAND LINE
my (@SEQS, $TABLE, $OUTFMT );

my $USAGE = "\n\n$0 [options]\n
Options:
  -seqs           Sequence file(s)
  -table          Codon usage table
  -outfmt         Output format in DNA (default) or RNA [optional]
  -help           Shows this message
\n";

# OPTIONS
GetOptions(
  'seqs=s{,}'     =>\@SEQS,
  'table:s'       =>\$TABLE,
  'outfmt:s'      =>\$OUTFMT,
  # ''            =>\,
  help            =>sub{pod2usage($USAGE);}
)or pod2usage(2);

checks(); # check CL arguments

#-------------------------------------------------------------------------------
# VARIABLES
my (%optimizedSeqs);
# push @SEQS, @DNA, @RNA;
# $SEQS{'DNA'} = \@DNA;
# $SEQS{'RNA'} = \@RNA;


#-------------------------------------------------------------------------------
# CALLS
my $codonUT = parseTable($TABLE);

foreach my $sequence (@SEQS){ # loop all sequence files
  my $inSeqIO   = Bio::SeqIO->new( -file => $sequence );

  while ( my $seqObj = $inSeqIO->next_seq) { # loop all sequences in file
    my $displayID = $seqObj->display_id;
    print "Optimizing $displayID ";
    # my $protObj = $seq->translate;
    # my $prot = $protObj->seq;

    my $optSeq = optimize($codonUT, $seqObj, $OUTFMT);
    $optimizedSeqs{$displayID} = $optSeq;
  }
}
# say Dumper(\%optimizedSeqs);


#-------------------------------------------------------------------------------
# SUBS


sub checks {
  unless ($TABLE) {
    die "Need to provide both sequence and codon usage table files.";
  }

}



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

  my (%codons, @optSeq);
  my %codonUT = %$codonUT;
  my $seq     = uc $seqObj->seq;
  my $seqLen  = $seqObj->length;
  my $type    = uc $seqObj->alphabet;

  say $type, " sequence";

  # Traverse through sequence &
  # Get best codon according to usage table
  for (my $i = 0; $i < $seqLen;) {
    my $codon = substr $seq, $i, 3;

    # DNA input must be replaced to RNA's codon usage table
    $codon =~ s/T/U/g if ($type eq 'DNA');
    my $aa = $codonUT{$codon}{'aa'};

    # Performance Boost
    # Check if AA already stored with best codon identified
    if ( exists $codons{$aa} ) {
    } else {
      $codon = getCodon($codonUT, $aa);
      $codons{$aa} = $codon;
    }

    # Make codon correct alphabet (DNA/RNA) for output
    if (!$OUTFMT) {
      $codon =~ s/U/T/g; # DNA default
    } else {
      uc $OUTFMT eq 'DNA' ? $codon =~ s/U/T/g : 0; # DNA or leave as RNA?
    }

    push @optSeq, $codon;
    $i = $i + 3;
  }

  return(\@optSeq);
}



sub getCodon {
  my ($codonUT, $aa) = @_;

  my %codonUT = %$codonUT;
  my ($fraction, $ultimate);
  my $best = 0;

  # Get all codons in usage table for specified AA
  my @codonsAA = getCodonsAA($codonUT, $aa);

  # Get ultimate best codon for specific AA
  # For now, best == codon with greatest fraction used for AA
  foreach (@codonsAA) {
    $fraction = $codonUT{$_}{'fraction'};
    if($fraction > $best) {
      $best = $fraction;
      $ultimate = $_;
    }
  }

  return($aa, $ultimate);
}



sub getCodonsAA {
  my ($codonUT, $aa) = @_;

  my %codonUT   = %$codonUT;
  my @triplets  = keys %codonUT;
  my @codonsAA  = grep { $codonUT{$_}{'aa'} eq $aa } @triplets;

  return(@codonsAA);
}
