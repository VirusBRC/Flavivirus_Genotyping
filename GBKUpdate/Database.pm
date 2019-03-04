package GBKUpdate::Database;

use strict;
use warnings;
our $AUTOLOAD;
use DBI;
use Exception::Class::DBI;
use Data::Dumper;
use Carp;
my $debug = 1;

sub new {
   my ($class, $conf) = @_;

   my $dsn = sprintf("DBI:mysql:%s:%s;max_allowed_packet=40MB",
                     $conf->get_dbname,
                     $conf->get_dbhost,
                    );
   my $db_user_name = $conf->get_dbuser;
   my $db_password = $conf->get_dbpasswd;
   my $dbh;
   $dbh = DBI->connect($dsn, $db_user_name, $db_password, 
     {PrintError  => 0, RaiseError  => 0, HandleError => Exception::Class::DBI->handler,}
    );
   $dbh->{'AutoCommit'} = 0;   # Turn on AutoCommit
   0 && $debug && print STDERR "Database: \$dbh is a '". ref($dbh) ."'\n";

   if (0) {
      my $sth = $dbh->prepare(qq{
        select id, name from family 
      });
      $sth->execute();

      my (@matrix) = ();
      while (my @ary = $sth->fetchrow_array())
      {
        push(@matrix, [@ary]);  # [@ary] is a reference
      }
      $sth->finish();
#      print Dumper(@matrix);
   }

   my $self = bless {dbh => $dbh}, $class;
   return $self;

}

# get genbank file by accession
sub get_genbank {
   my ($self, $accession) = @_;

   my $debug = 0;

   $debug && print STDERR "Database: get_genbank: accession=$accession\n";
   $debug && print STDERR "Database: get_genbank: \$dbh is a '". ref($self->{dbh}) ."'\n" if ($debug);

      my @columns = ('genbank');
      my $sth = join(',',@columns);
      $debug && print STDERR "Database: getFamily: sth=$sth\n";
      $sth = $self->{dbh}->prepare(qq{
         select $sth
           from genome
          where accession='$accession'
      });
      $sth->execute();

      my (@matrix) = ();
      while (my @ary = $sth->fetchrow_array()) {
        push(@matrix, [@ary]);  # [@ary] is a reference
      }
      $sth->finish();

#   print STDERR "Database: \$matrix[0]=\n". Dumper($matrix[0]) if ($debug);
   return $matrix[0];

}

sub getFamily {
   my ($self,$name) = @_;

   print STDERR "Database: getFamily: family name=$name\n";
   print STDERR "Database: getFamily: \$dbh is a '". ref($self->{dbh}) ."'\n" if ($debug);

      my @columns = ('id','name','taxid','min_length','max_length','query_string',
                'required_keywords','optional_keywords','reject_keywords');
      my $sth = join(',',@columns);
	print STDERR "Database: getFamily: sth=$sth\n";
      $sth = $self->{dbh}->prepare(qq{
         select $sth
           from family
          where name='$name'
      });
      $sth->execute();

      my (@matrix) = ();
      while (my @ary = $sth->fetchrow_array()) {
        push(@matrix, [@ary]);  # [@ary] is a reference
      }
      $sth->finish();

#   print STDERR "Database: \$matrix[0]=\n". Dumper($matrix[0]) if ($debug);
   return $matrix[0];

#    curs = self.connection.cursor()
#    stmt = """
#      select id,name,taxid,min_length,max_length,query_string,
#             required_keywords,optional_keywords,reject_keywords
#        from family
#       where name=%s"""
#    curs.execute(stmt,(name,))
#    if (curs.rowcount == 0):
#      return None
#    (id,na,tx,minl,maxl,qs,kw,op,re) = curs.fetchone()
#    curs.close()
#    return (int(id),na,tx,int(minl),int(maxl),qs,kw,op,re)
}

# returns a hashref containing all genomes with accession as key and
# array of version, gi and updatedate as value
sub getGenomes {
   my ($self,$family) = @_;

   print STDERR "Database: getGenomes: family name=$family\n";
   print STDERR "Database: getGenomes: \$dbh is a '". ref($self->{dbh}) ."'\n" if ($debug);

      my $sth = $self->{dbh}->prepare(qq{
      select accession,version,gi,updatedate
        from genome
        where family='$family'
      });
      $sth->execute();

      my $matrix;
      while (my @ary = $sth->fetchrow_array()) {
        $matrix->{$ary[0]} = [$ary[1], $ary[2], $ary[3]];  # [@ary] is a reference
      }
      $sth->finish();
#      print Dumper(@matrix);

#   print Dumper($matrix);
   return $matrix;

}

sub saveToDB {
    my ($self,$g,$exists) = @_;

#    curse = self.connection.cursor()
    my $stmt = "";
    $stmt .= "retrieved=now()";
    $stmt .= ",gi=?";
    $stmt .= ",version=?";
    $stmt .= ",title=?";
    $stmt .= ",taxid=?";
    $stmt .= ",genbank=?";
    $stmt .= ",fasta=?";
    $stmt .= ",createdate=?";
    $stmt .= ",updatedate=?";
    $stmt .= ",checksum=?";
    $stmt .= ",replaces=?";
    if ($exists) {
      $stmt .= ",status='changed'";
      $stmt = "
         update genome
            set $stmt
          where accession=?";
      print STDERR "Database: saveToDB: exists \$stmt=$stmt\n" if (0 && $debug);
#      curse.execute(stmt,())
      my $sth = $self->{dbh}->prepare_cached($stmt);
      $sth->execute(
              $g->getGI(),
              $g->getVersion(),
              $g->getTitle(),
              $g->getTaxon(),
              $g->getGBK(),
              $g->getFasta(),
              $g->getCreateDate(),
              $g->getUpdateDate(),
              $g->getChecksum(),
              $g->getReplaces(),
              $g->getAccession(),
              ) ;
      $self->recomputeIdentities($g->getAccession(),$g->getChecksum());
      $sth->finish();
    } else {
      $stmt .= ",accession=?";
      $stmt .= ",family=?";
      $stmt .= ",note=?";
      $stmt .= ",status='new'";
      $stmt = "
        insert into genome
          set $stmt
           ";
      print STDERR "Database: saveToDB: !exists \$stmt=$stmt\n" if (0 && $debug);
#      curse.execute(stmt,())

      my $sth = $self->{dbh}->prepare_cached($stmt);
      $sth->execute(
              $g->getGI(),
              $g->getVersion(),
              $g->getTitle(),
              $g->getTaxon(),
              $g->getGBK(),
              $g->getFasta(),
              $g->getCreateDate(),
              $g->getUpdateDate(),
              $g->getChecksum(),
              $g->getReplaces(),
              $g->getAccession(),
              $g->getFamily(),
              $g->getNote(),
              );
      $sth->finish();

      $self->recomputeIdentities($g->getAccession(),$g->getChecksum());
    }
#    curse.close()
    return;
}

sub markBad {
    my ($self,$g,$exists) = @_;

#    curse = self.connection.cursor()
    my $stmt = "";
    $stmt .= "gi=?";
    $stmt .= ",version=?";
    $stmt .= ",title=?";
    $stmt .= ",taxid=?";
    $stmt .= ",retrieved=now()";
    $stmt .= ",updatedate=?";
    $stmt .= ",createdate=?";

    if ($exists) {
      $stmt .= ",status='bad'";
      $stmt .= ",note=?";
      $stmt = "
        update genome
           set $stmt
         where accession = ?";
      print STDERR "Database: markBad: exists \$stmt=$stmt\n" if (0 && $debug);

      my $sth = $self->{dbh}->prepare_cached($stmt);
      $sth->execute(
              $g->getGI(),
              $g->getVersion(),
              $g->getTitle(),
              $g->getTaxon(),
              $g->getUpdateDate(),
              $g->getCreateDate(),
              $g->getNote(),
              $g->getAccession(),
              );
      $sth->finish();
#      curse.execute(stmt,(g.getGI(),g.getVersion(),g.getTitle(),g.getTaxon(),
#                          g.getUpdateDate(),g.getCreateDate(),
#                          g.getNote(),g.getAccession()))
    } else {
      $stmt .= ",status='bad'";
      $stmt .= ",family=?";
      $stmt .= ",note=?";
      $stmt .= ",accession=?";
      $stmt = "
        insert into genome
           set $stmt
          ";
      print STDERR "Database: markBad: !exists \$stmt=$stmt\n" if (0 && $debug);

      my $sth = $self->{dbh}->prepare_cached($stmt);
      $sth->execute(
              $g->getGI(),
              $g->getVersion(),
              $g->getTitle(),
              $g->getTaxon(),
              $g->getUpdateDate(),
              $g->getCreateDate(),
              $g->getFamily(),
              $g->getNote(),
              $g->getAccession(),
              );
      $sth->finish();
#      curse.execute(stmt,(g.getGI(),g.getVersion(),
#                          g.getTitle(),
#                          g.getTaxon(),g.getUpdateDate(),g.getCreateDate(),
#                          g.getFamily(),g.getNote(),g.getAccession()))
    }
#    curse.close()
    return;
}

sub recomputeIdentities {
    my ($self,$accession,$checksum) = @_;

    my $stmt = "delete from identical where src = ?";
    my $curse = $self->{dbh}->prepare_cached($stmt);
    $curse->execute($accession);
    my $src = $curse->rows;
    $curse->finish();
    $stmt = "delete from identical where dst = ?";
    $curse = $self->{dbh}->prepare_cached($stmt);
    $curse->execute($accession);
    my $dst = $curse->rows;
    $curse->finish();

    my $findstmt = "
       select accession from genome
       where checksum = ? and accession <> ?";
    $curse = $self->{dbh}->prepare_cached($findstmt);
    $curse->execute($checksum,$accession);
    my $matches = [];
    print STDERR sprintf("Database: recomputeIdentities: deleted src=$src dst=$dst\n");
    print STDERR sprintf("Database: recomputeIdentities: \$curse->rows=%s\n",$curse->rows);
    for my $i (1 .. $curse->rows) {
      my @ac = $curse->fetchrow_array();
      print STDERR "Database: recomputeIdentities: \@ac=@ac\n";
      push(@$matches,$ac[0]);
    }
    $curse->finish();

    my $addstmt = "insert into identical (src,dst) values (?,?)";
    $curse = $self->{dbh}->prepare_cached($addstmt);
    for my $ac (@$matches) {
      $curse->execute($ac,$accession);
      $curse->execute($accession,$ac);
    }
    $curse->finish();
    return;
}

sub saveTaxon {
   my ($self,$taxon) = @_;
    my $stmt = "
      insert into taxon
               set 
                id       =?
               ,speciesid=?
               ,genusid  =?
               ,familyid =?
               ,species  =?
               ,genus    =?
               ,family   =?
               ";
#    try:
#      curse = dbh.connection.cursor()
      my $curse = $self->{dbh}->prepare_cached($stmt);
      $curse->execute(
                      $taxon->{id},
                      $taxon->{speciesid},
                      $taxon->{genusid},
                      $taxon->{familyid},
                      $taxon->{species},
                      $taxon->{genus},
                      $taxon->{family},
                     );
      my $src = $curse->rows;
#      curse.close()
      $curse->finish();
#    except StandardError, zork:
#      log.error(zork)
    return;
}

sub commit {
    my ($self,$com) = @_;

#    self.connection.commit()
    if (!$com) {
       print STDERR "Database: commit: debugging, rollback\n";
       $self->{dbh}->rollback;
    } else {
       print STDERR "Database: commit: actual, commit\n";
       $self->{dbh}->commit;
    }
    return;
}

################################################################################
#
#class Database(object):
#
#  def __init__(self,conf):
#    self.connection = MySQLdb.connect(host=conf.get('dbhost'),
#                                      port=int(conf.get('dbport')),
#                                      user=conf.get('dbuser'),
#                                      passwd=conf.get('dbpasswd'),
#                                      db=conf.get('dbname'))
#
#  def getFamily(self,name):
#    curs = self.connection.cursor()
#    stmt = """
#      select id,name,taxid,min_length,max_length,query_string,
#             required_keywords,optional_keywords,reject_keywords
#        from family
#       where name=%s"""
#    curs.execute(stmt,(name,))
#    if (curs.rowcount == 0):
#      return None
#    (id,na,tx,minl,maxl,qs,kw,op,re) = curs.fetchone()
#    curs.close()
#    return (int(id),na,tx,int(minl),int(maxl),qs,kw,op,re)
#
#  def getGenomes(self,family):
#    """ retrieve genomes for a family, given family ID.  Result is a hash
#        with accessions as keys, and tuples as values: version, gi, date in
#        GenBank file, status """
#    curs = self.connection.cursor()
#    stmt = """
#      select accession,version,gi,updatedate
#        from genome
#        where family=%s"""
#    curs.execute(stmt,(family,))
#    data = {}
#    for i in range(0,curs.rowcount):
#      (ac,ve,gi,gb) = curs.fetchone()
#      data[ac] = (ve,gi,gb)
#    curs.close()
#    return data
#
#  def markBad(self,g,exists):
#    curse = self.connection.cursor()
#    if exists:
#      stmt = """
#        update genome
#          set status='bad',
#              gi=%s,
#              version=%s,
#              title=%s,
#              taxid=%s,
#              retrieved=now(),
#              updatedate=%s,
#              createdate=%s,
#              note=%s
#          where accession = %s"""
#      curse.execute(stmt,(g.getGI(),g.getVersion(),g.getTitle(),g.getTaxon(),
#                          g.getUpdateDate(),g.getCreateDate(),
#                          g.getNote(),g.getAccession()))
#    else:
#      stmt = """
#        insert
#          into genome (accession,gi,version,title,taxid,retrieved,
#                       updatedate,createdate,status,family,note)
#          values (%s,%s,%s,%s,%s,now(),%s,%s,'bad',%s,%s)"""
#      curse.execute(stmt,(g.getAccession(),g.getGI(),g.getVersion(),
#                          g.getTitle(),
#                          g.getTaxon(),g.getUpdateDate(),g.getCreateDate(),
#                          g.getFamily(),g.getNote()))
#    curse.close()
#
#  def recomputeIdentities(self,accession,checksum,curse):
#    curse.execute("delete from identical where src = %s",(accession,))
#    curse.execute("delete from identical where dst = %s",(accession,))
#    addstmt = "insert into identical (src,dst) values (%s,%s)"
#    findstmt = """
#      select accession from genome
#       where checksum = %s and accession <> %s"""
#    curse.execute(findstmt,(checksum,accession))
#    matches = []
#    for i in range(0,curse.rowcount):
#      (ac,) = curse.fetchone()
#      matches.append(ac)
#    for ac in matches:
#      curse.execute(addstmt,(ac,accession))
#      curse.execute(addstmt,(accession,ac))
#
#  def saveToDB(self,g,exists):
#    curse = self.connection.cursor()
#    if exists:
#      stmt = """
#        update genome
#          set status='changed',
#              gi=%s,
#              version=%s,
#              title=%s,
#              taxid=%s,
#              retrieved=now(),
#              createdate=%s,
#              updatedate=%s,
#              genbank=%s,
#              fasta=%s,
#              checksum=%s,
#              replaces=%s
#          where accession=%s"""
#      curse.execute(stmt,(g.getGI(),g.getVersion(),g.getTitle(),g.getTaxon(),
#                          g.getCreateDate(),g.getUpdateDate(),g.getGBK(),
#                          g.getFasta(),g.getChecksum(),
#                          g.getReplaces(),g.getAccession()))
#    else:
#      stmt = """
#        insert into genome
#          (accession,version,gi,family,title,taxid,retrieved,status,
#           note,genbank,fasta,updatedate,createdate,checksum,replaces)
#        values (%s,%s,%s,%s,%s,%s,now(),'new',%s,%s,%s,%s,%s,%s,%s)"""
#      curse.execute(stmt,(g.getAccession(),g.getVersion(),g.getGI(),
#                    g.getFamily(),g.getTitle(),g.getTaxon(),g.getNote(),
#                    g.getGBK(),
#                    g.getFasta(),g.getUpdateDate(),g.getCreateDate(),
#                    g.getChecksum(),g.getReplaces()))
#    self.recomputeIdentities(g.getAccession(),g.getChecksum(),curse)
#    curse.close()
#
#  def commit(self):
#    self.connection.commit()
#


1;
