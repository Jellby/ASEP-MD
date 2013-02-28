#!/usr/bin/perl -w

use strict;

# Use 2 or 3 arguments:
#  1 (read):  Gaussian output file
#  2 (read):  Gaussian fchk file
#  3 (write): generic output file (default: append ".gen" to Gaussian output)

my ($gau_output, $gau_fchk, $gen_output);
($gau_output, $gau_fchk, $gen_output) = @ARGV;

$gen_output = "$gau_output.gen" unless ($gen_output);

#===========================
# Read the data from the Gaussian output file

my $Q_natoms;
my $Q_charge;
my $Q_multiplicity;
my $Q_energy;
my $Q_energy_lower;
my $Q_selfenergy;
my @Q_dipole;
my @Q_mulliken;
my @Q_esp;
my @Q_gradient;
my @Q_gradient_lower;
my @Q_hessian;
my $Q_potfile;
my @Q_potential;

my $aux;

open GAU_OUTPUT, $gau_output or die "No se puede leer el fichero $gau_output";
open GAU_FCHK, $gau_fchk or die "No se puede leer el fichero $gau_fchk";

# Read the data available in the fchk file, in atomic units
while (<GAU_FCHK>) {
  if ($_ =~ /^Number of atoms/) {
    @_ = split;
    $Q_natoms = $_[4];
  }
  if ($_ =~ /^Charge/) {
    @_ = split;
    $Q_charge = $_[2];
  }
  if ($_ =~ /^Multiplicity/) {
    @_ = split;
    $Q_multiplicity = $_[2];
  }
  if ($_ =~ /^Total Energy/) {
    @_ = split;
    $Q_energy = $_[3];
  }
  if ($_ =~ /^Dipole Moment/) {
    $_ = <GAU_FCHK>;
    @Q_dipole = split;
  }
  if ($_ =~ /^Cartesian Gradient/) {
    @_ = split;
    $aux = int($_[4]/5)+($_[4]%5?1:0);
    $_ = "";
    for (my $i=0; $i<$aux; $i++) {
      $_ .= <GAU_FCHK>;
    }
    @Q_gradient = split;
  } 
  if ($_ =~ /^Cartesian Force Constants/) {
    @_ = split;
    $aux = int($_[5]/5)+($_[5]%5?1:0);
    $_ = "";
    for (my $i=0; $i<$aux; $i++) {
      $_ .= <GAU_FCHK>;
    }
    @Q_hessian = split;
  }
}

# Read the data in the Gaussian output file
while (<GAU_OUTPUT>) {
  if ($_ =~ /Self energy of the charges/) {
    @_ = split;
    $Q_selfenergy = $_[6];
    $Q_energy = $Q_energy-$Q_selfenergy;
  }
  # Energy difference and gradients for a conical intersection
  if ($_ =~ /Energy difference/) {
    @_ = split;
    $Q_energy_lower = $Q_energy+$_[2];
  }
  if ($_ =~ /^\s*Gradient of iOther State/) {
    undef @Q_gradient_lower;
    for (my $i=0; $i<$Q_natoms; $i++) {
      $_ = <GAU_OUTPUT>;
      push @Q_gradient_lower, split;
    }
  }
  if ($_ =~ /^\s*Gradient of iVec State/) {
    undef @Q_gradient;
    for (my $i=0; $i<$Q_natoms; $i++) {
      $_ = <GAU_OUTPUT>;
      push @Q_gradient, split;
    }
  }
  # Mulliken charges
  if (($_ =~ /^\s*Total atomic charges/) || ($_ =~ /^\s*Mulliken atomic charges/)) {
    $_ = <GAU_OUTPUT>;
    for (my $i=0; $i<$Q_natoms; $i++) {
      $_ = <GAU_OUTPUT>;
      @_ = split;
      $Q_mulliken[$i] = $_[2];
    }
  }
  # ESP charges
  if ($_ =~ /Charges from ESP fit,/) {
    $_ = <GAU_OUTPUT> for 1 .. 2;
    for (my $i=0; $i<$Q_natoms; $i++) {
      $_ = <GAU_OUTPUT>;
      @_ = split;
      $Q_esp[$i] = $_[2];
    }
  }
  if ($_ =~ /ESP Charges \(Molden\)/) {
    for (my $i=0; $i<$Q_natoms; $i++) {
      $_ = <GAU_OUTPUT>;
      @_ = split;
      $Q_esp[$i] = $_[1];
    }
  }
  # Fortran unit where the electrostatic potential is written
  if ($_ =~ /Compute potential derivative range/) {
    @_ = split;
    $Q_potfile = $_[11];
  }
}

close(GAU_OUTPUT);
close(GAU_FCHK);

# Read the potential
if ($Q_potfile) {
  open POTFILE, "fort.$Q_potfile" or die "No se puede leer el fichero fort.$Q_potfile";
  foreach (<POTFILE>) {
    @_ = split;
    push @Q_potential, $_[3];
  }
  close(POTFILE);
}

#===========================
# Write the generic output file, in atomic units

open GEN_OUTPUT, ">$gen_output"; 

printf GEN_OUTPUT "Number of atoms\n%4d\n\n",
  $Q_natoms;
printf GEN_OUTPUT "Charge\n%4d\n\n",
  $Q_charge;
printf GEN_OUTPUT "Multiplicity\n%4d\n\n",
  $Q_multiplicity;
printf GEN_OUTPUT "Energy\n%20.12E\n\n",
  $Q_energy;
printf GEN_OUTPUT "Energy (lower state)\n%20.12E\n\n",
  $Q_energy_lower if (defined $Q_energy_lower);
printf GEN_OUTPUT "Dipole moment\n%20.12f %20.12f %20.12f\n\n",
  @Q_dipole;
# These formats write 1 or 5 values per line...
printf GEN_OUTPUT "Mulliken charges";
  printf GEN_OUTPUT "%s%20.12f", (($_%1)?" ":"\n"), $Q_mulliken[$_] for 0 .. $#Q_mulliken;
  print GEN_OUTPUT "\n\n";
if (@Q_esp) {
  printf GEN_OUTPUT "ESP charges";
    printf GEN_OUTPUT "%s%20.12f", (($_%1)?" ":"\n"), $Q_esp[$_] for 0 .. $#Q_esp;
    print GEN_OUTPUT "\n\n";
}
if (@Q_gradient) {
  printf GEN_OUTPUT "Cartesian gradient";
    printf GEN_OUTPUT "%s%20.12E", (($_%5)?" ":"\n"), $Q_gradient[$_] for 0 .. $#Q_gradient;
    print GEN_OUTPUT "\n\n";
}
if (@Q_gradient_lower) {
  printf GEN_OUTPUT "Cartesian gradient (lower state)";
    printf GEN_OUTPUT "%s%20.12E", (($_%5)?" ":"\n"), $Q_gradient_lower[$_] for 0 .. $#Q_gradient_lower;
    print GEN_OUTPUT "\n\n";
}
if (@Q_hessian) {
  printf GEN_OUTPUT "Cartesian Hessian";
    printf GEN_OUTPUT "%s%20.12E", (($_%5)?" ":"\n"), $Q_hessian[$_] for 0 .. $#Q_hessian;
    print GEN_OUTPUT "\n\n";
}
if (@Q_potential) {
  printf GEN_OUTPUT "Electrostatic potential";
    printf GEN_OUTPUT "%s%20.12E", (($_%5)?" ":"\n"), $Q_potential[$_] for 0 .. $#Q_potential;
    print GEN_OUTPUT "\n\n";
}

close(GEN_OUTPUT);

exit 0;
