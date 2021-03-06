package algo::HxNyAdvanced;
use strict;
use warnings;
use English;

use Data::Dumper;
use algo::Algo;
use Genotype_newickBI;
use Genotype_util;

my $debug = 1;

#=========Inheriting Algo=============#
our @ISA = qw(algo::Algo);

# sub getClassification
sub getClassification {
  my $self = shift;
  my ( $treeFile, $accession ) = @_;
  my $subn = "getClassification";
  my $debug = 1 || Genotype_util::getDebugAll();

  my $nodelook = $self->Nodelook;
  my $lookup   = {};

  $debug && print STDERR "$subn: sequence file AND query :: $treeFile, $accession\n";
  my ( $X, $Sister ) = $self->Advance0( $treeFile, $accession );

  $debug && print STDERR "$subn: X :: " . join( '|', @{$X} ) . "\n";
  $debug && print STDERR "$subn: Sister :: " . join( '|', @{$Sister} ) . "\n";

  foreach ( @{$X} ) {
    chomp $_;
    $debug && print STDERR "$subn: X     ='$_'\n";
    if ( $_ =~ m/^[A-Za-z]+(\d+)$/ ) {
      $_ =~ s/\s//g;
      my $type = $self->getLookup($_);
      if ( defined($type) ) {
        $type =~ s/\s//g;
        $lookup->{$_} = $type;
        $debug && print STDERR "$subn: X     ='$_' type=$type\n";
      }
    }
  }
  foreach ( @{$Sister} ) {
    chomp $_;
    $debug && print STDERR "$subn: Sister='$_'\n";
    if ( $_ =~ m/^[A-Za-z]+(\d+)$/ ) {
      $_ =~ s/\s//g;
      my $type = $self->getLookup($_);
      if ( defined($type) ) {
        $type =~ s/\s//g;
        $lookup->{$_} = $type;
        $debug && print STDERR "$subn: Sister='$_' type=$type\n";
      }
    }
  }
  my $confirm = undef;
  my $final   = undef;

  $debug && print STDERR "$subn: lookup=\n".Dumper($self->{lookup})."End of lookup\n";
  $debug && print STDERR "$subn: lookup=\n".Dumper($lookup)."End of lookup\n";
  $debug && print STDERR "$subn:   Lookup:\n";
  foreach my $acc ( sort keys %{$lookup} ) {
    $debug && print STDERR "$subn:     $acc = " . $lookup->{$acc} . "\n";
  }

  my ( $needconfirmation, $Provision, $prob ) =
    $self->Classification( $X, $lookup, $nodelook, $accession );
  $debug && print STDERR "$subn: Provision=$Provision, prob=$prob\n";
  if ($needconfirmation) {

    $debug && print STDERR "$subn: Need Confirmation: Provision=$Provision, prob=$prob\n";

    $confirm = $self->Find_value( $Sister, $lookup, $nodelook );
    $final = $self->Final( $Provision, $confirm, $nodelook );

  }
  else {
    $confirm = ' no need for sister confirmation ';
    $final   = $Provision;
  }
  unless ( defined($final) ) { $final = 'ND'; }
  $debug && print STDERR "$subn: Provision, confirm, final ::: $Provision\t$confirm\t$final\n";

  return ( $final, $prob );
} # sub getClassification

############################################################
# sub callCladinator, call cladinator program by C Zmasek
sub callCladinator {
  my $self = shift;
  my ( $treeFile, $accession ) = @_;
  my $subn = "callCladinator";
  my $debug = 0 || Genotype_util::getDebugAll();

  my ( $final, $prob ) = ('NA', 0.000); # Default values
  my $cmd = "";
  my $fileCladinatorOutput = 'cladinator_out_1.tsv';
  if ($treeFile =~ m/^([^.]+)[.]/) {
      $fileCladinatorOutput = $accession . '_cladinator.tsv';
  }
  if (-e $fileCladinatorOutput) {
      `rm $fileCladinatorOutput`; # Remove cladinator output if already exists
  }
  $cmd = $self->config->getValue("cladinator");
  $cmd .= " $treeFile";
  $cmd .= " $fileCladinatorOutput";
  $debug && print STDERR "$subn: accession=$accession cmd='$cmd'\n";
  system("$cmd");
  my $status = $?;
  die "Failed to run cladinator for $accession ($status)\n" if ($status);

  # Read result from cladinator, and parse the result
  # V2.0.1, Nov 09, 2017: Only take 'Matching Clades'
  # V2.0.2, Nov 17, 2017: Take 'Matching Clades', further take 'Specific-hits' if prob>=0.7
  # V2.0.3, Nov 28, 2017: Changed parsing algorithm: BRC-11607
  my $PROB_MIN = 0.7;
  my ( $match, $matchProb ) = ('', 0.000);
  my ( $spec, $specProb ) = ('', 0.000);
  open my $infile, '<', "$fileCladinatorOutput"
       or croak("$subn: '$fileCladinatorOutput', couldn't open: $OS_ERROR");
  while (my $line = <$infile>) {
      next if ($line =~ /^#/ || $line !~ /Matching Clades|Specific-hits/);
      $debug && print STDERR "$subn: line=$line";
      my $words = [ split(/\t/, $line) ];
      if ($line =~ /Matching Clades/ && $words->[3] >= $PROB_MIN && $words->[3] > $matchProb) {
        ($match, $matchProb) = ($words->[2], $words->[3]);
      } elsif ($line =~ /Specific-hits/ && $words->[3] >= $PROB_MIN && $words->[3] > $specProb) {
        ($spec, $specProb) = ($words->[2], $words->[3]);
      }
  }
  close $infile
      or croak("$subn: '$fileCladinatorOutput', couldn't close: $OS_ERROR");

  # Changed parsing algorithm
  if ($specProb >= $PROB_MIN) {
      ( $final, $prob ) = ($spec, $specProb);
  } elsif ($matchProb >= $PROB_MIN) {
      ( $final, $prob ) = ($match, $matchProb);
      $final = 'NA' if ($final eq '?');
  } else {
      ( $final, $prob ) = ('NA', 1.000);
  }

  unless ( defined($final) ) { $final = 'ND'; }
  $debug && print STDERR "$subn: final, prob ::: \t$final\t$prob\n";

  return ( $final, $prob, $fileCladinatorOutput);
} # sub callCladinator

############################################################
# sub findGenotype, find genotype by 
sub findGenotype {
  my $self = shift;
  my ( $treeFile, $accession ) = @_;
  my $subn = "findGenotype";
  my $debug = 0 || Genotype_util::getDebugAll();

  my $nodelook = $self->Nodelook;
  my $lookup   = {};
  my $lookupSister   = {};

  $debug && print STDERR "$subn: sequence file AND query :: $treeFile, $accession\n";
  my ( $sibId, $extId ) = $self->Advance( $treeFile, $accession );

  $debug && print STDERR "$subn: lookup=\n".Dumper($self->{lookup})."End of lookup\n";
  $debug && print STDERR "$subn: \$sibId=\n".Dumper($sibId)."End of \$sibId\n";
  $debug && print STDERR "$subn: \$extId=\n".Dumper($extId)."End of \$extId\n";

  my $X = [];
  foreach (keys %{$sibId}) {
    push @{$X}, @{$sibId->{$_}};
  }
  my $Sister = [];
  foreach (keys %{$extId}) {
    push @{$Sister}, @{$extId->{$_}};
  }
  foreach ( @{$X} ) {
    my $type = '';
    chomp $_;
    if ( $_ =~ m/^[A-Za-z_]+(\d+)([a-zA-Z_]+(\d+)?)?$/ ) {
      $_ =~ s/\s//g;
      $type = $self->getLookup($_);
      if ( defined($type) ) {
        $type =~ s/\s//g;
        $lookup->{$_} = $type;
      }
    }
    $debug && print STDERR "$subn: X     ='$_'".(($type)?"\ttype=$type":'')."\n";
  }
  foreach ( @{$Sister} ) {
    my $type = '';
    chomp $_;
    if ( $_ =~ m/^[A-Za-z_]+(\d+)([a-zA-Z_]+(\d+)?)?$/ ) {
      $_ =~ s/\s//g;
      next if (exists($lookup->{$_})); # Skip node if it's already included in $X
      my $type = $self->getLookup($_);
      if ( defined($type) ) {
        $type =~ s/\s//g;
        $lookupSister->{$_} = $type;
      }
    }
    $debug && print STDERR "$subn: Sister='$_'".(($type)?"\ttype=$type":'')."\n";
  }
  my $confirm = undef;
  my $final   = undef;

  $debug && print STDERR "$subn: lookup=\n".Dumper($lookup)."End of lookup\n";
  $debug && print STDERR "$subn: lookupSister=\n".Dumper($lookupSister)."End of lookupSister\n";
  $debug && print STDERR "$subn:   Lookup:\n";
  foreach my $acc ( sort keys %{$lookup} ) {
    $debug && print STDERR "$subn:     $acc = " . $lookup->{$acc} . "\n";
  }

  my ( $needconfirmation, $Provision, $prob );
#  ( $needconfirmation, $Provision, $prob ) =
#    $self->Classification( $X, $lookup, $nodelook, $accession );
  ( $needconfirmation, $Provision, $prob ) =
    $self->summBranch( $sibId, $lookup, $nodelook, $accession );
  $debug && print STDERR "$subn: Provision=$Provision, prob1=$prob\n";
  my ( $needconfirmation1, $Provision1, $prob1 );
#  ( $needconfirmation1, $Provision1, $prob1 ) =
#    $self->Classification( $Sister, $lookupSister, $nodelook, $accession );
#  $debug && print STDERR "$subn: Provision1=$Provision1, prob1=$prob1\n";

  if ($needconfirmation) {
    $debug && print STDERR "$subn: Need Confirmation: Provision=$Provision, prob=$prob\n";

    $confirm = $self->Find_value( $Sister, $lookup, $nodelook );
    $final = $self->Final( $Provision, $confirm, $nodelook );

  }
  else {
    $confirm = ' no need for sister confirmation ';
    $final   = $Provision;
  }
  unless ( defined($final) ) { $final = 'ND'; }
  $debug && print STDERR "$subn: Provision, confirm, final ::: $Provision\t$confirm\t$final\n";

  return ( $final, $prob );
} # sub findGenotype

############################################################
#sub summBranch, looks through the sibling IDs, to determine the genotype
# sibId = $VAR1 = {
#          '_seq001_#0_M=1' => [ 'HQ537008','HQ537009', ]
#         };
#lookup=  $VAR1 = {
#          'EU246939' => '6t',
#          'EU246931' => '6e'
#         };
sub summBranch {
  my $self = shift;
  my ( $sibId, $lookup, $nodelook, $query ) = @_;
  my $subn = "algo::HxNyAdvanced::summBranch";
  my $debug = 0 || Genotype_util::getDebugAll();

  $debug && print STDERR "$subn: Start\n";
  $debug && print STDERR "$subn: query=$query\tnodelook=\n".Dumper($nodelook)."End of nodelook\n";
  $debug && print STDERR "$subn: query=$query\tlookup=\n".Dumper($lookup)."End of lookup\n";
  $debug && print STDERR "$subn: query=$query\tsibId=\n".Dumper($sibId)."End of sibId\n";
  my $provisions = {};
  my $provision;
  my $needsister = 0;

  my $probs      = {};
  for my $qnode (sort keys %{$sibId}) {
    my $sibs = $sibId->{$qnode};
    my $prob = 0;
    $prob = ($qnode =~ m/_M=((\d)+(\.+(\d+)+([eE][-+]?(\d+))?)?)$/) ? $1 : 0;
    if ($prob<0.001) {
      print STDERR "$subn: [WARNING] we might have a problem in: qnode=$qnode\tprob=$prob\n";
    }
    $debug && print STDERR "$subn: qnode=$qnode\tprob=$prob\n";

    my $counts = {}; # counts = { '1b' => 1 };
    my $total = 0;
    for my $sib (@{$sibs}) {
      next if (!exists($lookup->{$sib}) || !defined($lookup->{$sib}));
      $counts->{ $lookup->{$sib} }++;
      $total++;
    }
    $debug && print STDERR "$subn: qnode=$qnode\tcounts=\n".Dumper($counts)."End of counts\n";

    for my $type (keys %{$counts}) {
      if (0) {
        # Distribute according to ratio of counts
        $provisions->{$type} += $prob * $counts->{$type} / $total;
      } else {
        # Distribute equally between all genotypes
        $provisions->{$type} += $prob / scalar(keys %{$counts});
      }
    }
    $debug && print STDERR "$subn: qnode=$qnode\tcounts=\n".Dumper($counts)."End of counts\n";
    $debug && print STDERR "$subn: qnode=$qnode\tprovisions=\n".Dumper($provisions)."End of provisions\n";
  }

  $debug && print STDERR "$subn: provisions=\n".Dumper($provisions)."End of provisions\n";
  # If there are more genotypes, delete ND
  if (scalar(keys %$provisions)>1 && exists($provisions->{'ND'})) {
    print STDERR "$subn: removed 'ND' from provisions=\n".Dumper($provisions)."End of provisions\n";
    delete $provisions->{'ND'};
  }

  # Pick the genotype w/ highest probability
  $provision = (scalar(keys %$provisions)>=1) ? (sort {$provisions->{$b}<=>$provisions->{$a}} keys %$provisions)[0]
             : 'ND';
  my $prob = 0;
  $prob = $provisions->{$provision} if (exists($provisions->{$provision}) && defined($provisions->{$provision}));

  # Collect the genotypes (1) and subtypes (1a)
  # GU208812
  my $underOneSubtype = 1;
  my %genoProbs = (); # genoProbs = { '2' => '1.0000001' };
  my @subtypes = ();  # subtypes  = [ '2b' ];
  for  my $t (sort keys %{$provisions})  {
    next if ($t !~ m/^(\d+)/); # outgroup
      $genoProbs{$1} += $provisions->{$t};
      push @subtypes, $t;
  }
  $debug && print STDERR "$subn: genoProbs=\n".Dumper(\%genoProbs)."End of genoProbs\n";
  $debug && print STDERR "$subn: subtypes=\n".Dumper(\@subtypes)."End of subtypes\n";

  # Check the genotypes and subtypes
  my $THRES_PROBABILITY_MIN = 0.714;
  my $countGenotypes = scalar(keys %genoProbs);
  my @keysGenotypes = (keys %genoProbs);
  my $countSubtypes = scalar(@subtypes);
  $debug && print STDERR "$subn: Found genoProbs[$countGenotypes]='@keysGenotypes'; subtypes[$countSubtypes]='@subtypes'\n";
  if ( $countGenotypes >= 2 
      # If there is >1 genotype, ND
       && $prob < $THRES_PROBABILITY_MIN
     ) {
      if ( $genoProbs{(sort {$genoProbs{$b}<=>$genoProbs{$a}} keys %genoProbs)[0]} > 0.714 ) {
        $provision = (sort {$genoProbs{$b}<=>$genoProbs{$a}} keys %genoProbs)[0];
        $prob = $genoProbs{$provision};
      } else {
        $provision = 'ND';
        $prob = 0;
      }
  } elsif ( $countGenotypes == 1 && $countSubtypes > 1 ) {
      # If 1 genotype, but >1 subtypes, take only genotype, not subtype
    if ( $prob < $THRES_PROBABILITY_MIN ) {
      $debug && print STDERR "$subn: Subtype w/ max prob. is $provision:$prob\n";
      $provision = (keys %genoProbs)[0];
      $prob = $genoProbs{$provision};
      $debug && print STDERR "$subn: Subtype changed to genotype: $provision:$prob\n";
    }
  }

  $debug && print STDERR "$subn: provisions=\n".Dumper($provisions)."End of provisions\n";
  $debug && print STDERR "$subn: provision=$provision prob=$prob\n";
  return ( $needsister, $provision, $prob );
} # sub summBranch

############################################################
## sub Advance gets the subtree that is related to the query
sub Advance {
  my $self = shift;
  my ( $tree, $query ) = @_;
  my $subn = "algo::HxNyAdvanced::Advance";
  my $debug = 0 || Genotype_util::getDebugAll();

  my $X      = [];
  my $Sister = [];
  my $array  = $self->clean($tree);

  $debug && print STDERR "$subn: whole tree : " . join('|', @{$array}) . "\n";
  $debug && print STDERR "$subn: $query\n";
  my $len         = scalar @{$array};
  my @count_array = ();
  my $XStart = $len;
  my $XEnd   = -1;
  my $SStart = $len;
  my $SEnd   = -1;
  my $sibId   = {};
  my $extId   = {};
  for ( my $i = 0 ; $i < $len ; ++$i ) {
    $array->[$i] =~ s/\s+//g;
    if ( $array->[$i] !~ m/^$query/ ) {
      next;
    }    # end if($array->[$i]=~m/$query/)

    $debug && print STDERR "$subn:\n";
    $debug && print STDERR "$subn: $i :: $i :: $array->[$i-1] $array->[$i] $array->[$i+1] \tX='@{$X}'\n";

    if ( $array->[ $i - 1 ] =~ m/\(/ ) {
        push( @{$X}, $array->[ $i - 1 ] );
        $debug && print STDERR "$subn: $i :: $i-1 :: $array->[$i-1] \tX='@{$X}'\n";
        my $count = 1;
        my $k     = $i - 1;
        do {
          ++$k;
          push( @{$X}, $array->[$k] );
          $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tX='@{$X}'\n";
          $count += ( $array->[$k] =~ m/\(/ ) ? +1
                  : ( $array->[$k] =~ m/\)/ ) ? -1
                  : 0;

          $debug && print STDERR "$subn: X :: $k :: $array->[$k] ::: $count\n";
        } until ( ( $array->[$k] =~ m/\)/ ) && ( $count == 0 ) );
        for (my $j=$k; $j>$i; $j--) {
          next if ($array->[$j] !~ m/^[a-zA-Z_]+(\d+)([a-zA-Z_]+(\d+)?)?$/);
          push @{ $sibId->{ $array->[$i] } }, $array->[$j];
        }
        $XStart = $i - 1 if ($XStart > $i - 1);
        $XEnd   = $k if ($XEnd   < $k);

        $debug && print STDERR "$subn: final X :: " . join( '', @{$X} ) . "\n";
        if ( $array->[ $k + 1 ] =~ m/\,/ ) {
          my $count_nxt = 1;
          do {
            ++$k;
            push( @{$Sister}, $array->[$k] );
            $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tSister='@{$Sister}'\n";
            $count_nxt += ( $array->[$k] =~ m/\(/ ) ? +1
                        : ( $array->[$k] =~ m/\)/ ) ? -1
                        : 0;

            $debug && print STDERR "$subn: knext_forward :: $k :: $array->[$k] ::: $count_nxt\n";
          } until ( ( $array->[$k] =~ m/\)/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          shift( @{$Sister} );
          $SStart = $i - 2 if ($SStart > $i - 2);
          $SEnd   = $k if ($SEnd   < $k);

          $debug && print STDERR "$subn: sister_forward1 :: " . join( '', @{$Sister} ) . "\n";

        } elsif ( $array->[ $i - 2 ] =~ m/\,/ ) {
          my $count_nxt = 1;
          my $point     = $i - 1;
          do {
            --$point;
            push( @{$Sister}, $array->[$point] );
            $debug && print STDERR "$subn: $i :: $point :: $array->[$point] \tSister='@{$Sister}'\n";
            $count_nxt += ( $array->[$point] =~ m/\)/ ) ? +1
                        : ( $array->[$point] =~ m/\(/ ) ? -1
                        : 0;
            $debug && print STDERR "$subn: knext_backward1 :: $k :: $array->[$k] ::: $count_nxt\n";
          } until ( ( $array->[$point] =~ m/\(/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          @{$Sister} = reverse @{$Sister};
          pop( @{$Sister} );
          $SStart = $i-2 if ($SStart < $i-2);
          $SEnd   = $point - 1 if ($SEnd   < $point - 1);

          $debug && print STDERR "sister_backward1 :: " . join( '', @{$Sister} ) . "\n";
        }

    } elsif ( $array->[ $i + 1 ] =~ m/\)/ ) {
        my $point = $i + 1;
        push( @{$X}, $array->[ $i + 1 ] );
        $debug && print STDERR "$subn: $i :: $i+1 :: $array->[$i+1] \tX='@{$X}'\n";
        my $count = 1;
        my $k     = $i + 1;
        do {
          --$k;
          push( @{$X}, $array->[$k] );
          $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tX='@{$X}'\n";
          $count += ( $array->[$k] =~ m/\)/ ) ? +1
                  : ( $array->[$k] =~ m/\(/ ) ? -1
                  : 0;
        } until ( ( $array->[$k] =~ m/\(/ ) && ( $count == 0 ) );
        @{$X} = reverse @{$X};
        for (my $j=$k; $j<$i; $j++) {
          next if ($array->[$j] !~ m/^[a-zA-Z_]+(\d+)([a-zA-Z_]+(\d+)?)?$/);
          push @{ $sibId->{ $array->[$i] } }, $array->[$j];
        }
        $XStart = $k if ($XStart > $k);
        $XEnd   = $i + 1 if ($XEnd   < $i + 1);

        $debug && print STDERR "$subn: final X :: " . join( '', @{$X} ) . "\n";

        if ( $array->[ $i + 2 ] =~ m/\)/ ) {
          my $count_nxt = 1;
          do {
            --$k;
            push( @{$Sister}, $array->[$k] );
            $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tSister='@{$Sister}'\n";
            $count_nxt += ( $array->[$k] =~ m/\)/ ) ? +1
                        : ( $array->[$k] =~ m/\(/ ) ? -1
                        : 0;
          } until ( ( $array->[$k] =~ m/\(/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          @{$Sister} = reverse @{$Sister};
          pop( @{$Sister} );
          $debug && print STDERR "$subn: $i :: Sister='@$Sister'\n";
          $SStart = $k if ($SStart > $k);
          $SEnd   = $i + 2 if ($SEnd   < $i + 2);

          $debug && print STDERR "$subn: sister_backward2 :: " . join( '', @{$Sister} ) . "\n";

        } elsif ( $array->[ $k - 1 ] =~ m/\(/ ) {
          my $count_nxt = 1;
          do {
            ++$point;
            push( @{$Sister}, $array->[$point] );
            $debug && print STDERR "$subn: $i :: $point :: $array->[$point] \tSister='@{$Sister}'\n";
            $count_nxt += ( $array->[$point] =~ m/\(/ ) ? +1
                        : ( $array->[$point] =~ m/\)/ ) ? -1
                        : 0;
          } until ( ( $array->[$point] =~ m/\)/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          shift( @{$Sister} );
          $debug && print STDERR "$subn: $i :: Sister='@$Sister'\n";
          $SStart = $k - 1 if ($SStart > $k - 1);
          $SEnd   = $point if ($SEnd   < $point);

          $debug && print STDERR "$subn: sister_forward2 :: " . join( '', @{$Sister} ) . "\n";
        }
    }
  }

  $debug && print STDERR "$subn: X= " . join( '', @{$X} ) . "\n";
  $debug && print STDERR "$subn: Sister= " . join( '', @{$Sister} ) . "\n";
  $X = [ @{ $array }[$XStart .. $XEnd] ];
  $Sister = [ @{ $array }[$SStart.. $XStart-1, $XEnd+1 .. $SEnd] ];
  $debug && print STDERR "$subn: X= " . join( '', @{$X} ) . "\n";
  $debug && print STDERR "$subn: Sister= " . join( '', @{$Sister} ) . "\n";
  $X = [ @{ $array }[$SStart .. $SEnd] ];
  $debug && print STDERR "$subn: X= " . join( '', @{$X} ) . "\n";

  $debug && print STDERR "$subn: \$sibId=\n".Dumper($sibId)."End of \$sibId\n";
  $debug && print STDERR "$subn: \$extId=\n".Dumper($extId)."End of \$extId\n";

  return ( $sibId, $extId );
} # sub Advance

############################################################
## sub Advance gets the subtree that is related to the query
# Tom's version, more-or-less
sub Advance0 {
  my $self = shift;
  my ( $tree, $query ) = @_;
  my $subn = "algo::HxNyAdvanced::Advance0";
  my $debug = 0 || Genotype_util::getDebugAll();

  my $X      = [];
  my $Sister = [];
  my $array  = $self->clean($tree);

  print STDERR "$subn: whole tree : " . join('|', @{$array}) . "\n";
  print STDERR "$subn: $query\n";
  my $len         = scalar @{$array};
  my @count_array = ();
  for ( my $i = 0 ; $i < $len ; ++$i ) {
    $array->[$i] =~ s/\s+//g;
    if ( $array->[$i] =~ m/$query/ ) {

      $debug && print STDERR "$subn:\n";
      $debug && print STDERR "$subn: $i :: $i :: $array->[$i-1] $array->[$i] $array->[$i+1] \tX='@{$X}'\n";
      if ( $array->[ $i - 1 ] =~ m/\(/ ) {
        push( @{$X}, $array->[ $i - 1 ] );
        $debug && print STDERR "$subn: $i :: $i-1 :: $array->[$i-1] \tX='@{$X}'\n";
        my $count = 1;
        my $k     = $i - 1;
        do {
          ++$k;
          push( @{$X}, $array->[$k] );
          $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tX='@{$X}'\n";
          if ( $array->[$k] =~ m/\(/ ) {
            ++$count;
          }
          elsif ( $array->[$k] =~ m/\)/ ) {
            --$count;
          }

          # else{}
          $debug && print STDERR "$subn: X :: $k :: $array->[$k] ::: $count\n";
        } until ( ( $array->[$k] =~ m/\)/ ) && ( $count == 0 ) );

        print STDERR "$subn: final X :: " . join( '', @{$X} ) . "\n";
        if ( $array->[ $k + 1 ] =~ m/\,/ ) {
          my $count_nxt = 1;
          do {
            ++$k;
            push( @{$Sister}, $array->[$k] );
            $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tSister='@{$Sister}'\n";
            if ( $array->[$k] =~ m/\(/ ) {
              ++$count_nxt;
            }
            elsif ( $array->[$k] =~ m/\)/ ) {
              --$count_nxt;
            }
            else { }

            $debug && print STDERR "$subn: knext_forward :: $k :: $array->[$k] ::: $count_nxt\n";
          } until ( ( $array->[$k] =~ m/\)/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          shift( @{$Sister} );

          print STDERR "$subn: sister_forward1 :: " . join( '', @{$Sister} ) . "\n";
        }
        elsif ( $array->[ $i - 2 ] =~ m/\,/ ) {
          my $count_nxt = 1;
          my $point     = $i - 1;
          do {
            --$point;
            push( @{$Sister}, $array->[$point] );
            $debug && print STDERR "$subn: $i :: $point :: $array->[$point] \tSister='@{$Sister}'\n";
            if ( $array->[$point] =~ m/\)/ ) {
              ++$count_nxt;
            }
            elsif ( $array->[$point] =~ m/\(/ ) {
              --$count_nxt;
            }

            print STDERR "$subn: knext_backward1 :: $k :: $array->[$k] ::: $count_nxt\n";
          } until ( ( $array->[$point] =~ m/\(/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          @{$Sister} = reverse @{$Sister};
          pop( @{$Sister} );

          print STDERR "sister_backward1 :: " . join( '', @{$Sister} ) . "\n";
        }
      }
      elsif ( $array->[ $i + 1 ] =~ m/\)/ ) {
        my $point = $i + 1;
        push( @{$X}, $array->[ $i + 1 ] );
        $debug && print STDERR "$subn: $i :: $i+1 :: $array->[$i+1] \tX='@{$X}'\n";
        my $count = 1;
        my $k     = $i + 1;
        do {
          --$k;
          push( @{$X}, $array->[$k] );
          $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tX='@{$X}'\n";
          if ( $array->[$k] =~ m/\)/ ) {
            ++$count;
          } elsif ( $array->[$k] =~ m/\(/ ) {
            --$count;
          } else { }
        } until ( ( $array->[$k] =~ m/\(/ ) && ( $count == 0 ) );
        @{$X} = reverse @{$X};

        print STDERR "$subn: final X :: " . join( '', @{$X} ) . "\n";

        if ( $array->[ $i + 2 ] =~ m/\)/ ) {
          my $count_nxt = 1;
          do {
            --$k;
            push( @{$Sister}, $array->[$k] );
            $debug && print STDERR "$subn: $i :: $k :: $array->[$k] \tSister='@{$Sister}'\n";
            if ( $array->[$k] =~ m/\)/ ) {
              ++$count_nxt;
            } elsif ( $array->[$k] =~ m/\(/ ) {
              --$count_nxt;
            } else { }
          } until ( ( $array->[$k] =~ m/\(/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          @{$Sister} = reverse @{$Sister};
          pop( @{$Sister} );
          $debug && print STDERR "$subn: $i :: Sister='@$Sister'\n";

          print STDERR "$subn: sister_backward2 :: " . join( '', @{$Sister} ) . "\n";
        }
        elsif ( $array->[ $k - 1 ] =~ m/\(/ ) {
          my $count_nxt = 1;
          do {
            ++$point;
            push( @{$Sister}, $array->[$point] );
            $debug && print STDERR "$subn: $i :: $point :: $array->[$point] \tSister='@{$Sister}'\n";
            if ( $array->[$point] =~ m/\(/ ) {
              ++$count_nxt;
            }
            elsif ( $array->[$point] =~ m/\)/ ) {
              --$count_nxt;
            }
            else { }
          } until ( ( $array->[$point] =~ m/\)/ ) && ( $count_nxt == 0 ) );
          pop( @{$Sister} );
          shift( @{$Sister} );
          $debug && print STDERR "$subn: $i :: Sister='@$Sister'\n";

          print STDERR "$subn: sister_forward2 :: " . join( '', @{$Sister} ) . "\n";
        }
      }
    }    # end if($array->[$i]=~m/$query/)
  }

  return ( $X, $Sister );
} # sub Advance0

############################################################
sub clean {
  my $self = shift;
  my ($file) = @_;
  my $subn = "algo::HxNyAdvanced::clean";

  my $array = [];
  open IN, "$file" or die "cannot open file:$!";
  while (<IN>) {
    chomp $_;
    my @arr1 = ();
#    print STDERR "$subn: 0='$_'\n";
    $_ =~ s/:/ /g;
#    print STDERR "$subn: 1='$_'\n";
    $_ =~ s/,/ , /g;
#    print STDERR "$subn: 2='$_'\n";
#    $_ =~ s/(\d)+\.+(\d+)+([eE][-+]?(\d+))?//g; # Tom's code
    $_ =~ s/(?<!_M[=])(\d)+\.+(\d+)+([eE][-+]?(\d+))?//g;
#    print STDERR "$subn: 3='$_'\n";
#    $_ =~ s/(\d)[eE][-+]?\d+//g; # What is the purpose for this?
#    print STDERR "$subn: 4='$1' matched string\n";
#    print STDERR "$subn: 4='$_'\n";
    $_ =~ s/\)/ \) /g;
#    print STDERR "$subn: 5='$_'\n";
    $_ =~ s/\(/ \( /g;
#    print STDERR "$subn: 6='$_'\n";
    @arr1 = split( " ", $_ );
    push( @{$array}, @arr1 );
  }
  return $array;
}


############################################################
sub Classification {
  my $self = shift;
  my ( $X, $lookup, $nodelook, $query ) = @_;
  my $subn = "Classification";
  my $debug = 0 || Genotype_util::getDebugAll();

  print STDERR "$subn: Start\n";
  my $provisions = {};
  my $provision;
  my $len_X      = $#{$X};
  my $needsister = 0;
  my $thisx      = 0;
  my $thisxSet   = [];
  my $probs      = {};
  $debug && print STDERR "$subn: len_X :: " . ( $len_X - 1 ) . " @{$X}\n";
  for ( my $k = $len_X - 1 ; $k > 0 ; --$k ) {
    $debug && print STDERR "$subn: query=$query \tk=$k \t$X->[$k]\n";
#    if ( $X->[$k] =~ m/^$query$/ ) {
    if ( $X->[$k] =~ m/^$query/ ) {
      $debug && print STDERR "$subn: query $k \n";
      $thisx = $k;
      push @{$thisxSet}, $k;
      $probs->{$X->[$k]} = 1;
    }
  }
  $debug && print STDERR "$subn: thisx :: $thisx\n";
  $debug && print STDERR "$subn: thisxSet=\n".Dumper($thisxSet)."End of thisxSet\n";
  $debug && print STDERR "$subn: probs=\n".Dumper($probs)."End of probs\n";

for $thisx (@{$thisxSet}) {
  my $prob = 0;
  $prob = ($X->[ $thisx ] =~ m/_M=((\d)+\.+(\d+)+([eE][-+]?(\d+))?)$/) ? $1 : 1;
  $prob *= $probs->{$X->[ $thisx ]};
  $probs->{$X->[ $thisx ]} = 0;
  $debug && print STDERR "$subn: x=$thisx query=$X->[ $thisx ] \tprob=$prob\n";

  if ( $X->[ $thisx + 1 ] eq ',' ) {
#    if ( $X->[ $thisx + 2 ] =~ m/^[a-zA-Z]+(\d+)$/ ) {
#     #$provision  = $lookup->{ $X->[ $thisx + 2 ] };
#      $provisions->{ $lookup->{ $X->[ $thisx + 2 ] } } += $prob;
#      $needsister = $needsister || 0;
#    } elsif ( $X->[ $thisx + 2 ] eq ')' ) {
#    }
    my $k     = $thisx + 2;
    if ( $X->[ $k ] eq ')' ) {
        my $count = 1;
        do {
          ++$k;
          $count += ( $X->[$k] =~ m/\)/ ) ? -1
                  : ( $X->[$k] =~ m/\(/ ) ? +1
                  : 0;
          $debug && printf STDERR "$subn: thisx=$thisx count=$count $k=%-6s \tX='@{$X}'\n", $X->[$k];
        } until ( (( $X->[$k] =~ m/^[a-zA-Z]+(\d+)$/ ) && ( $count == 1 )) || $k==$#{$X} );
        
    }

    if ( $X->[ $k ] =~ m/^[a-zA-Z]+(\d+)$/ ) {
#      $provision  = $lookup->{ $X->[ $thisx - 2 ] };
      if ( exists($lookup->{ $X->[ $k ] }) ) {
        $provisions->{ $lookup->{ $X->[ $k ] } } += $prob;
        $needsister = $needsister || 0;
      }
    }

  }
  elsif ( $X->[ $thisx - 1 ] eq ',' ) {
    my $k     = $thisx - 2;
    if ( $X->[ $k ] eq ')' ) {
        my $count = 1;
        do {
          --$k;
          $count += ( $X->[$k] =~ m/\)/ ) ? +1
                  : ( $X->[$k] =~ m/\(/ ) ? -1
                  : 0;
          $debug && printf STDERR "$subn: thisx=$thisx count=$count $k=%-6s \tX='@{$X}'\n", $X->[$k];
        } until ( (( $X->[$k] =~ m/^[a-zA-Z]+(\d+)$/ ) && ( $count == 1 )) || $k==0 );
        
    }

    if ( $X->[ $k ] =~ m/^[a-zA-Z]+(\d+)$/ ) {
#      $provision  = $lookup->{ $X->[ $thisx - 2 ] };
      if ( exists($lookup->{ $X->[ $k ] }) ) {
        $provisions->{ $lookup->{ $X->[ $k ] } } += $prob;
        $needsister = $needsister || 0;
      }
    }

  } else {
    $needsister = $needsister || 1;
    if ( $thisx == 1 ) {
      splice( @{$X}, $len_X, 1 );
      splice( @{$X}, 0,      3 );

      $debug && print "$subn: begin :: " . join( '', @{$X} ) . "\n";
#      $provision = $self->sub_find_value( $X, $lookup, $nodelook );
      $provisions->{ $self->sub_find_value( $X, $lookup, $nodelook ) } += 1;
    }
    elsif ( $thisx == $len_X - 1 ) {
      splice( @{$X}, $len_X - 2, 3 );
      splice( @{$X}, 0,          1 );

      $debug && print STDERR "$subn: end :: " . join( '', @{$X} ) . "\n";
#      $provision = $self->sub_find_value( $X, $lookup, $nodelook );
      $provisions->{ $self->sub_find_value( $X, $lookup, $nodelook ) } += 1;
    }
    else {
      $debug && print STDERR "$subn: error in X :: " . join( '', @{$X} ) . "\n";
#      $provision = 'ND';
      $provisions->{ 'ND' } += 1;
    }
  }
}

  $debug && print STDERR "$subn: provisions=\n".Dumper($provisions)."End of provisions\n";
  # If there are more genotypes, delete ND
  delete $provisions->{'ND'} if (scalar(keys %$provisions)>1 && exists($provisions->{'ND'}));
  # Pick the genotype w/ highest probability
  $provision = (scalar(keys %$provisions)>=1) ? (sort {$provisions->{$b}<=>$provisions->{$a}} keys %$provisions)[0]
             : 'ND';
  my $prob = $provisions->{$provision};
  $prob = 0 if (!defined($prob));

  $debug && print STDERR "$subn: provisions=\n".Dumper($provisions)."End of provisions\n";
  print STDERR "$subn: provision=$provision prob=$prob\n";
  return ( $needsister, $provision, $prob );
} # sub Classification

############################################################
sub Find_value {
  my $self = shift;
  my ( $X, $hash_def, $nodelook ) = @_;

  print "Start Find_value\n";
  print "  Branch = " . join( '', @{$X} ) . "\n";
  my $value  = undef;
  my @values = ();
  print "  Assign Clades to Nodes in Branch\n";
  for ( my $i = 0 ; $i <= $#{$X} ; $i++ ) {
    if ( $X->[$i] =~ m/^[a-zA-Z]+(\d+)$/ ) {
      if ( defined( $hash_def->{ $X->[$i] } ) ) {
        $values[$i] = $hash_def->{ $X->[$i] };
        $value = $values[$i];
      }
      else {
        print "error :: $X->[$i] not defined\n";
        $value = 'ND';
        print "Find_value return :: $value\n";
        return $value;
      }
    }
    else { $values[$i] = $X->[$i]; }
  }
  print "  value = $value\n";
  print "  Branch Decorated with Clades = " . join( '', @values ) . "\n";
  my $cycle = 0;
  my $pair;
  print "  Cycling through the Branch\n";
  do {
    $cycle++;
    print "    Cycle = $cycle\n";
    $pair = 1;
    for ( my $i = 1 ; $i < $#values ; $i++ ) {
      if ( $values[ $i - 1 ] eq '('
        && $values[ $i + 1 ] eq ','
        && $values[ $i + 3 ] eq ')'
        && $values[$i] =~ m/^[a-zA-Z0-9.-]+/
        && $values[ $i + 2 ] =~ m/^[a-zA-Z0-9.-]+/ )
      {
        print "      Pair($pair) :: index = $i, "
          . $values[$i] . ", "
          . $values[ $i + 2 ] . "\n";
        $value = $nodelook->{ $values[$i] }->{ $values[ $i + 2 ] };
        print "      value = $value\n";
        unless ( defined($value) ) {
          print
            "error ::: no node lookup sister ::: $values[$i] , $values[$i+2]\n";
          $value = 'ND';
          print "Find_value return :: $value\n";
          return $value;
        }
        splice( @values, $i - 1, 4 );
        $values[ $i - 1 ] = $value;
        $pair++;

        print "      Sister Array :: " . join( '', @values ) . "\n\n";
      }
    }
  } until ( $#values <= 1 || $pair == 1 );

  print "Find_value return :: $value\n";
  return $value;
}

sub sub_find_value {
  my $self = shift;
  my ( $X, $hash_def, $nodelook ) = @_;

  print "Start sub_find_value\n";
  print "  Branch = " . join( '', @{$X} ) . "\n";
  my $value  = undef;
  my @values = ();
  print "  Assign Clades to Nodes in Branch\n";
  for ( my $i = 0 ; $i <= $#{$X} ; $i++ ) {
    if ( $X->[$i] =~ m/^[a-zA-Z]+(\d+)$/ ) {
      if ( defined( $hash_def->{ $X->[$i] } ) ) {
        $values[$i] = $hash_def->{ $X->[$i] };
        $value = $values[$i];
      }
      else {
        print "error :: $X->[$i] not defined\n";
        $value = 'ND';
        print "sub_find_value return :: $value\n";
        return $value;
      }
    }
    else { $values[$i] = $X->[$i]; }
  }
  print "  value = $value\n";
  print "  Branch Decorated with Clades = " . join( '', @values ) . "\n";
  my $cycle = 0;
  my $pair;
  print "  Cycling through the Branch\n";
  do {
    $cycle++;
    print "    Cycle = $cycle\n";
    $pair = 1;
    for ( my $i = 1 ; $i < $#values ; $i++ ) {
      if ( $values[ $i - 1 ] eq '('
        && $values[ $i + 1 ] eq ','
        && $values[ $i + 3 ] eq ')'
        && $values[$i] =~ m/^[a-zA-Z0-9.-]+/
        && $values[ $i + 2 ] =~ m/^[a-zA-Z0-9.-]+/ )
      {
        print "      Pair($pair) :: index = $i, "
          . $values[$i] . ", "
          . $values[ $i + 2 ] . "\n";
        $value = $nodelook->{ $values[$i] }->{ $values[ $i + 2 ] };
        print "      value = $value\n";
        unless ( defined($value) ) {
          print "error ::: no node lookup X ::: $values[$i] , $values[$i+2]\n";
          $value = 'ND';
          print "sub_find_value return :: $value\n";
          return $value;
        }
        splice( @values, $i - 1, 4 );
        $values[ $i - 1 ] = $value;
        $pair++;

        print "      Sister Array :: " . join( '', @values ) . "\n\n";
      }
    }
  } until ( $#values <= 1 || $pair == 1 );

  print "sub_find_value return :: $value\n";
  return $value;
}

sub Final {
  my $self = shift;
  my ( $Provisional, $Confirmation, $nodelook ) = @_;
  my $finaldef = $nodelook->{$Provisional}->{$Confirmation};
  print
"Final :: provision = $Provisional, confirm = $Confirmation, final = $finaldef\n";
  return $finaldef;
}

sub Nodelook {
  my $self = shift;
  my $subn = 'Nodelook';
  my $debug = 0;

  my $node_lookup = $self->{node_lookup};
  open IN, "$node_lookup" or die "cannot open file:$!";
  my $nodelook = {};
  while (<IN>) {
    chomp;
    my @split2;
    unless ($_) { next; }

    # print "node : $_\n";
    @split2 = split( "\t", $_ );
    $split2[0] =~ s/\s+//g;
    $split2[1] =~ s/\s+//g;
    $split2[2] =~ s/\s+//g;
    $nodelook->{ $split2[0] }->{ $split2[1] } = $split2[2];
    $nodelook->{ $split2[1] }->{ $split2[0] } = $split2[2];
  }
  close IN;
  $debug && print STDERR "$subn: nodelook=\n".Dumper($nodelook)."End of nodelook\n";
  return $nodelook;
}

1;
