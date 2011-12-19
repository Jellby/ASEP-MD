#!/usr/bin/perl -w

use strict;

# Use 1 or 2 arguments:
#  1 (read):  generic input file
#  2 (write): Molcas input file (default: append ".mol" to generic input)

my ($gen_input, $mol_input);
($gen_input, $mol_input) = @ARGV;

$mol_input = "$gen_input.mol" unless ($mol_input);

#===========================
# Read the data from the generic input file

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
# Write the Molcas input file

my $angstrom = 1.88972613289;

open MOL_INPUT, ">$mol_input";

# Write the SEWARD INPUT
print MOL_INPUT <<"EOF";
&SEWARD &END
  Title
    Job Name
  >>> coord inline basis 6-31G
  NoSymm
EOF

# Write the molecular geometry, in Angstrom
print MOL_INPUT "    $Q_natoms\n\n";
for (my $i=0; $i<$Q_natoms; $i++) {
  printf MOL_INPUT "    %16s %20.12f %20.12f %20.12f\n",
    $Q_molecule[$i]{type}, $Q_molecule[$i]{x}, $Q_molecule[$i]{y}, $Q_molecule[$i]{z};
}
print MOL_INPUT "  End of input\n";

# Write the external charges, in Angstrom, or the charge file
if ($Q_charges) {
  print MOL_INPUT "  XField\n";
  print MOL_INPUT "    $Q_charges 1\n";
  for (my $i=0; $i<$Q_charges; $i++) {
    printf MOL_INPUT "%20.12f %20.12f %20.12f %20.12f %3.1f %3.1f %3.1f\n",
      $Q_charges_points[$i]{x}, $Q_charges_points[$i]{y}, $Q_charges_points[$i]{z}, $Q_charges_points[$i]{q}, 0, 0, 0;
  }
} elsif ($Q_chargesfile) {
  print MOL_INPUT "  XField\n";
  printf MOL_INPUT "\@\$MOLCAS_SUBMIT_PWD/%s\n", $Q_chargesfile;
}

# Write the potential points, in Bohr
if ($Q_potential) {
  print MOL_INPUT "  EPot\n";
  print MOL_INPUT "    $Q_potential\n";
  for (my $i=0; $i<$Q_potential; $i++) {
    printf MOL_INPUT "%20.12f %20.12f %20.12f\n",
      $Q_potential_points[$i]{x}*$angstrom, $Q_potential_points[$i]{y}*$angstrom, $Q_potential_points[$i]{z}*$angstrom;
  }
}
print MOL_INPUT "End of input\n";

# Write the SCF input
print MOL_INPUT <<"EOF";

&SCF &END
End of input
EOF

my $electrons = 0;
$electrons += $Q_molecule[$_]{number} foreach (0 .. $#Q_molecule);
my $orbitals = $electrons/2;
# Write the RASSCF input
# Change as needed
#	print MOL_INPUT <<"EOF";
#	
#	&RASSCF &END
#	  Inactive = $orbitals
#	  Ras2 = 0
#	  NActElectrons = 0 0 0
#	End of input
#	EOF

# Write the ALASKA/MCKINLEY input
if ($Q_derivative >= 1) {
  print MOL_INPUT <<"EOF";

&ALASKA &END
End of input
EOF
}
if ($Q_derivative >= 2) {
  print MOL_INPUT <<"EOF";

&MCKINLEY &END
End of input
!cat \$Project.UnSym
EOF
}

close(MOL_INPUT);

exit 0;
