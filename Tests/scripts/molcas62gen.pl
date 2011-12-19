#!/usr/bin/perl -w

use strict;

# Use 1 or 2 arguments:
#  1 (read):  Molcas output file
#  2 (write): generic output file (default: append ".gen" to Molcas output)

my ($mol_output, $gen_output);
($mol_output, $gen_output) = @ARGV;

$gen_output = "$mol_output.gen" unless ($gen_output);

#===========================
# Read the data from the Molcas output file

# Root(s) to read
# Change as needed
my $Root=1;
my $Root_lower=0; #set to 0 if not an intersection

my $Q_natoms;
my $Q_charge;
my $Q_multiplicity;
my $Q_energy;
my $Q_energy_lower;
my @Q_dipole;
my @Q_mulliken;
my @Q_esp;
my @Q_gradient;
my @Q_gradient_lower;
my @Q_hessian;

my $aux;
my $module;
my $root;
my $alaskaroot;
my $debye=0.3934303073;

open MOL_OUTPUT, $mol_output or die "No se puede leer el fichero $mol_output";

while (<MOL_OUTPUT>) {
  # Update which Molcas module the output belongs to
  if ($_ =~ /executing module (\w*) /) {
    $module = $1;
    undef $root;
    $root = $Root if ($module eq "SCF");
  }
  # Select the root to which the values refer
  if ($_ =~ /for (root number:|state)\s*(\d*)/) {
    $root = $2;
  }
  if (($_ =~ /for root no\. =\s*(\d*)/) && ($module eq "MCLR")) {
    $alaskaroot = $1;
  }
  # Read the number of atoms (from SEWARD)
  if (($_ =~ /Cartesian Coordinates/) && ($module eq "SEWARD")) {
    $_ = <MOL_OUTPUT> for 1 .. 3;
    $Q_natoms = 0;
    while (1) {
      $_ = <MOL_OUTPUT>;
      last if ($_ =~ /^\s*$/);
      $Q_natoms++;
    }
  }

  # Read the charge (from SCF or RASSCF)
  if (($_ =~ /Molecular charge/) && ($module eq "SCF")) {
    @_ = split;
    $Q_charge = int($_[2]+0.5);
  }
  if (($_ =~ /Total            charge/) && ($module eq "RASSCF")) {
    @_ = split;
    $Q_charge = int($_[2]+0.5);
  }

  # Read the multiplicity (from SCF or RASSCF)
  if (($_ =~ /Fermi aufbau/) && ($module eq "SCF")) {
    $_ = <MOL_OUTPUT>;
    @_ = split;
    $Q_multiplicity = $_[1];
    if ($_[0] =~ /alpha/) {
      $_ = <MOL_OUTPUT>;
      @_ = split;
    }
    $Q_multiplicity = abs($Q_multiplicity-$_[1])+1;
  }
  if (($_ =~ /Spin quantum number/) && ($module eq "RASSCF")) {
    @_ = split;
    $Q_multiplicity = int(2*$_[3]+1.5);
  }

  # Read the energy (from SCF or RASSCF or CASPT2)
  if (($_ =~ /Total (SCF|KS-DFT) energy/) && ($module eq "SCF")) {
    @_ = split;
    $Q_energy = $_[3];
  }
  if (($_ =~ /Final state energy\(ies\)/) && ($module eq "RASSCF")) {
    $_ = <MOL_OUTPUT> for 1 .. 2;
    until ($_ =~ /^\s*$/) {
      $_ = <MOL_OUTPUT>;
      @_ = split;
      $Q_energy = $_[5] if ($#_ > 1 && $_[2] == $Root);
      $Q_energy_lower = $_[5] if ($#_ > 1 && $_[2] == $Root_lower);
    }
  }
  if (($_ =~ /^\s+Total energy/) && ($module eq "CASPT2")) {
    @_ = split;
    $Q_energy = $_[2] if ($root == $Root);
    $Q_energy_lower = $_[2] if ($root == $Root_lower);
  }

  # Read the Mulliken charges (from SCF, RASSC or CASPT2)
  if (($_ =~ /Mulliken charges/) && ($module eq "SCF") && ($root == $Root)) {
    undef @Q_mulliken;
    while (1) {
      $_ = <MOL_OUTPUT>;
      $_ = <MOL_OUTPUT>;
      last if ($_ =~ /^\s+Total electronic charge/);   
      until ($_ =~ /^\s*$/) {
        $_ = <MOL_OUTPUT>;
      }
      $_ = <MOL_OUTPUT>;
      @_ = split;
      shift @_;
      push @Q_mulliken, @_;
    }
  }

  # Read the ESP charges
  if ($_ =~ /ESP Charges \(Molden\)/) {
    for (my $i=0; $i<$Q_natoms; $i++) {
      $_ = <MOL_OUTPUT>;
      @_ = split;
      $Q_esp[$i] = $_[1];
    }
  }

  # Read the dipole moment, in Debye (from SCF, RASSCF or CASPT2)
  if (($_ =~ /^Dipole Moment/) && ($module eq "SCF") && ($root == $Root)) {
    $_ = <MOL_OUTPUT>;
    $_ = <MOL_OUTPUT>;
    @_ = split;
    @Q_dipole = map($_*$debye, @_[1,3,5]);
  }

  # Read the Cartesian gradient (from ALASKA)
  if (($_ =~ /Molecular gradients/) && ($module eq "ALASKA")) {
    if ((not defined($alaskaroot)) || ($alaskaroot == $Root)) {
      $_ = <MOL_OUTPUT> for 1 .. 5;
      $aux = $Q_natoms*3;
      for (my $i=0; $i<$aux; $i++) {
        $_ = <MOL_OUTPUT>;
        @_ = split;
        $Q_gradient[$i] = $_[2];
      }
    } elsif ($alaskaroot == $Root_lower) {
      $_ = <MOL_OUTPUT> for 1 .. 5;
      $aux = $Q_natoms*3;
      for (my $i=0; $i<$aux; $i++) {
        $_ = <MOL_OUTPUT>;
        @_ = split;
        $Q_gradient_lower[$i] = $_[2];
      }
    }
  } 

  # Read the Cartesian Hessian (from the UnSym file after MCKINLEY or MCLR)
  if (($_ =~ /BEGIN HESSIAN/) && (($module eq "MCKINLEY") || ($module eq "MCLR"))) {
    undef @Q_hessian;
    $_ = <MOL_OUTPUT>;
    $aux = int($Q_natoms*3/4)+($Q_natoms*3%4?1:0);
    for (my $i=0; $i<$Q_natoms*3; $i++) {
      $_ = <MOL_OUTPUT>;
      $_ = "";
      for (my $j=0; $j<$aux; $j++) {
        $_ .= <MOL_OUTPUT>;
      }
      @_ = split;
      push @Q_hessian, @_[0 .. $i];
    }
  }
}

#  # Fortran unit where the electrostatic potential is written
#  if ($_ =~ /Compute potential derivative range/) {
#    @_ = split;
#    $Q_potfile = $_[11];
#  }

close(MOL_OUTPUT);

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

#if ($Q_potfile) {
#  open POTFILE, "fort.$Q_potfile" or die "No se puede leer el fichero fort.$Q_potfile";
#  printf GEN_OUTPUT "Electrostatic potential\n";
#    print GEN_OUTPUT foreach (<POTFILE>);
#    print GEN_OUTPUT "\n";
#  close(POTFILE);
#}

close(GEN_OUTPUT);

exit 0;
