package GBKUpdate::Configuration;

use strict;
use warnings;
use English;
our $AUTOLOAD;
use Carp;
my $debug = 1;

{
   my $_count = 0;
   sub get_count {
      return $_count; 
   }
   sub _incr_count {
      return ++$_count;
   }
   sub _decr_count {
      return --$_count;
   }

}

################################################################################
#
#class Configuration(object):
#  """ Reads and parses the configuration file. """
#
#  configPat = "\s*(\w+)\s*=\s*([^\s]+)\s*$"
#
#  def __init__(self,path):
#    self.pairs = {}
#    self._load(path)
#
#  def get(self,key):
#    return self.pairs[key]
#
#  def _load(self,path):
#    fd = file(path)
#    for line in fd:
#      mo = re.match(Configuration.configPat,line[:-1]);
#      if mo:
#        self.pairs[mo.group(1)] = mo.group(2)
#      else:
#        log.warn("skipping illegal line in config file '%s': '%s'" % (path,line[:-1]))
#    fd.close()
#

################################################################################
#
sub new {
   my ($class, $fn) = @_;

#   my $configPat = '/^\s*(\w+)\s*=\s*([^\s]+)\s*$/';

#   GBKUpdate::Message->info("Configuration: openning config file: $fn",1);
   my $configFILE;
   if (-r $fn) {
      open $configFILE, '<', $fn
#           or die $!;
           or croak("Configuration: Located $fn, but couldn't open");
#      print STDERR <$configFILE>;
   } else {
      croak("Couldn't located $fn");
   }

   my $data = {};
   while (<$configFILE>) {
#	print $_;
#	if ($_ =~ $configPat) {
	if ($_ =~ /^\s*(\w+)\s*=\s*([^\s]*)\s*$/) {
	   0 && print "\$1=\'$1\' \t\$2=\'$2\'\n";
	   my $key = '_'. $1;
	   my $value = $2;
#	   print "$key=$value\n";
	   $data->{$key} = $value;
#	   print "Configuration: $key='".$data->{$key}."'\n" if ($debug);
	} else {
	   print STDERR "Configuration: skipping illegal line in $fn: $_";
	   croak("Configuration: skipping illegal line in $fn: $_");
	}
   }
   close $configFILE or croak "Configuration: Couldn't close $fn: $OS_ERROR";

#   print keys %$data;
   my $self = bless $data, $class;

   $class->_incr_count();

#   GBKUpdate::Message->info("Configuration: finished reading config file: $fn",1);
   return $self;
}


sub AUTOLOAD {
   my ($self, $newvalue) = @_;
   my ($operation, $attribute) = ($AUTOLOAD =~ /(get|set)(_\w+)/);

   # Is this a legal method name?
   unless ($operation && $attribute) {
	croak "Configuration: Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
   }
   unless (exists $self->{$attribute}) {
	croak "Configuration: No such attribute '$attribute' exists in the class", ref($self);
   }

   # Turn off strict references to enable "magic" AUTOLOAD speedup
   no strict 'refs';

   # AUTOLOAD accessors
   if ($operation eq 'get') {
	# define subroutine
	*{$AUTOLOAD} = sub { shift->{$attribute} };

   # AUTOLOAD mutators
   } elsif ($operation eq 'set') {
	# define subroutine
	*{$AUTOLOAD} = sub { shift->{$attribute} = shift; };
	# Set the new attribute value
	$self->{$attribute} = $newvalue;
   }

   # Turn strict references back on
   use strict 'refs';

   # Return attribute value
   return $self->{$attribute};

}


sub DESTROY {
   my $self = @_;
   GBKUpdate::Configuration->_decr_count();
#   $self->_decr_count();
}

# Original python script
#
#  def __init__(self,path):
#    self.pairs = {}
#    self._load(path)
#
#  def get(self,key):
#    return self.pairs[key]
#
#  def _load(self,path):
#    fd = file(path)
#    for line in fd:
#      mo = re.match(Configuration.configPat,line[:-1]);
#      if mo:
#        self.pairs[mo.group(1)] = mo.group(2)
#      else:
#        log.warn("skipping illegal line in config file '%s': '%s'" % (path,line[:-1]))
#    fd.close()
#

1;
