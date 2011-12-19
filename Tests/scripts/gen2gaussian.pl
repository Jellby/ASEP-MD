#!/usr/bin/perl -w

use strict;

# Use 1 or 2 arguments:
#  1 (read):  generic input file
#  2 (write): Gaussian input file (default: append ".gau" to generic input)

my ($gen_input, $gau_input);
($gen_input, $gau_input) = @ARGV;

$gau_input = "$gen_input.gau" unless ($gau_input);

#===========================
# Read the data from the generic input file

my $Q_potfile1 = 63; #default fortran unit for potential input
my $Q_potfile2 = 64; #default fortran unit for potential output
my $Q_natoms;
my @Q_molecule;
my $Q_derivative;
my $Q_potential;
my @Q_potential_points;
my $Q_chargesfile;
my $Q_charges;
my @Q_charges_points;

open GEN_INPUT, $gen_input or die "No se puede leer el fichero $gen_input";

while (<GEN_INPUT>) {
  # Read the geometry
  #  atom name, atomic symbol, atomic number, coordinates in Angstrom
  if ($_ =~ /^Geometry$/) {
    $_ = <GEN_INPUT>;
    @_ = split;
    $Q_natoms = $_[0];
    for (my $i=0; $i<$Q_natoms; $i++) {
      $_ = <GEN_INPUT>;
      @_ = split;
      $Q_molecule[$i]{name} = $_[0];
      $Q_molecule[$i]{type} = $_[1];
      $Q_molecule[$i]{number} = $_[2];
      $Q_molecule[$i]{x} = $_[3];
      $Q_molecule[$i]{y} = $_[4];
      $Q_molecule[$i]{z} = $_[5];
    }
  }
  # Read the order of the energy derivative needed
  if ($_ =~ /^Derivative$/) {
    $_ = <GEN_INPUT>;
    @_ = split;
    $Q_derivative = $_[0];
  }
  # Read the external charges positions and values
  if ($_ =~ /^External charges$/) {
    $_ = <GEN_INPUT>;
    @_ = split;
    $Q_charges = $_[0];
    for (my $i=0; $i<$Q_charges; $i++) {
      $_ = <GEN_INPUT>;
      @_ = split;
      $Q_charges_points[$i]{x} = $_[0];
      $Q_charges_points[$i]{y} = $_[1];
      $Q_charges_points[$i]{z} = $_[2];
      $Q_charges_points[$i]{q} = $_[3];
    }
  }
  # ... or read the filename containing external charges
  if ($_ =~ /^External charge file$/) {
    $_ = <GEN_INPUT>;
    @_ = split;
    $Q_chargesfile = $_[0];
  }
  # Read the points where the electrostatic potential is to be calculated
  if ($_ =~ /^Potential points$/) {
    $_ = <GEN_INPUT>;
    @_ = split;
    $Q_potential = $_[0];
    for (my $i=0; $i<$Q_potential; $i++) {
      $_ = <GEN_INPUT>;
      @_ = split;
      $Q_potential_points[$i]{x} = $_[0];
      $Q_potential_points[$i]{y} = $_[1];
      $Q_potential_points[$i]{z} = $_[2];
    }
  }
}

#===========================
# Write the Gaussian input file

open GAU_INPUT, ">$gau_input";

# Header
# Change as needed
print GAU_INPUT <<"EOF";
#P HF/6-31G
GFInput iop(6/7=3) Pop=CHELPG
FCheck=All NoSymm
EOF

printf GAU_INPUT "%s\n", "Force" if ($Q_derivative == 1);
printf GAU_INPUT "%s\n", "Freq " if ($Q_derivative == 2);
printf GAU_INPUT "%s\n", "Prop=(Grid,Potential)" if ($Q_potential);
printf GAU_INPUT "%s\n", "Charge=Angstrom" if ($Q_charges || $Q_chargesfile);

# Job name, charge and multiplicity
# Change as needed
print GAU_INPUT <<"EOF";

Job Name

0 1
EOF

# Molecular geometry, in Angstrom
for (my $i=0; $i<$Q_natoms; $i++) {
  printf GAU_INPUT "%16s %20.12f %20.12f %20.12f\n",
    $Q_molecule[$i]{type}, $Q_molecule[$i]{x}, $Q_molecule[$i]{y}, $Q_molecule[$i]{z};
}
print GAU_INPUT "\n";

# External charges, in Angstrom
if ($Q_charges) {
  for (my $i=0; $i<$Q_charges; $i++) {
    printf GAU_INPUT "%20.12f %20.12f %20.12f %20.12f\n",
      $Q_charges_points[$i]{x}, $Q_charges_points[$i]{y}, $Q_charges_points[$i]{z}, $Q_charges_points[$i]{q};
  }
  print GAU_INPUT "\n";
} elsif ($Q_chargesfile) {
  printf GAU_INPUT "@%s/N\n\n", $Q_chargesfile;
}

# Potential points
if ($Q_potential) {
  # Number of points and fortran units
  printf GAU_INPUT "%i 3 %i %i\n\n", $Q_potential, $Q_potfile1, $Q_potfile2;
  # Coordinates in the fort.XX file, in Angstrom
  open GAU_POT, ">fort.$Q_potfile1";
  for (my $i=0; $i<$Q_potential; $i++) {
    printf GAU_POT "%20.12f %20.12f %20.12f\n",
      $Q_potential_points[$i]{x}, $Q_potential_points[$i]{y}, $Q_potential_points[$i]{z};
  }
  close(GAU_POT);
}

close(GAU_INPUT);

exit 0;
