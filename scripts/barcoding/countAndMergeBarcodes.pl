#!/bin/perl -w
#
# Need to go through ALL fastq files,
# look for bonafide barcode-containing reads, count them, and
# depth-normalize.  Then we can sum up the normalized reads for
# all barcodes across ALL samples, and plot the density distribution
#

use strict;

# Set up global variables for results
my %barcodes = ();
my %allBarcodes = ();
my %samples = ();


# Routine for Hamming distance
sub hd {
   return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

# Subroutine to check for expected library codes, allowing ONE mismatch
sub hasLibCode {

   my $lc1 = "CCAA";
   my $lc2 = "ACGT";
   my $lc3 = "TGGA";

   my $libCode = $_[0];
   if (hd($libCode,$lc1) <= 1) {
      return $lc1;
   } elsif (hd($libCode,$lc2) <= 1) {
      return $lc2;
   } elsif (hd($libCode,$lc3) <= 1) {
      return $lc3;
   } 

   return 0;
}

# Subroutine to check for expected vector sequence, allowing ONE mismatch
sub hasVector {

   if (hd($_[0],"ATCGATAC") <= 1) {
      return 1;
   }

   return 0;
}



# Routine for parsing a FASTQ file.
# IN:  filename to parse
# OUT:  Nothing.  Stores BC in %barcodes hashtable and %samples
sub readFile {

   my $numReads = 0;
   my $keptReads = 0;
   my $fName = $_[0];

   # Grab sample name from directory name just above filename:
   my @path = split("/",$fName);
   my $sample = $path[ $#path ];
   $sample =~ s/.fastq.gz//;

   # Store sample name in samples hashtable
   $samples{ $sample } = 1;

   open (IN,"zcat -c $fName | ") or die "Could not open input $fName";

   while (my $line1 = <IN>) {

      my $line2 = <IN>;
      my $line3 = <IN>;
      my $line4 = <IN>;

      $numReads++;

      my $libCode = substr($line2,14,4);
      my $vector = substr($line2,18,8);

      $libCode = hasLibCode($libCode);
      if ($libCode && hasVector($vector)) {

         $keptReads++;
         
         my $barcode = substr($line2,0,14).$libCode;

         $allBarcodes{$barcode} = 1;
         if (defined($barcodes{$sample}{$barcode})) {
            $barcodes{$sample}{$barcode} = $barcodes{$sample}{$barcode}+1;
         } else {
            $barcodes{$sample}{$barcode} = 1;
         }
      }

   }
   #print "Read in $numReads barcode reads...\n";

   #return %localBarcodes;
   print $fName."\t".$numReads."\t".$keptReads."\t".($keptReads/$numReads)."\n";
}


# Finally, need a function to write out the results
sub writeBarcodeMatrix {

   my $outFile = $_[0];
   open (OUT,">$outFile");

   my @keys = keys %samples;
   print OUT "Barcode\t".join("\t",keys %samples)."\n";
   
   foreach my $barcode (keys %allBarcodes) {
      print OUT $barcode;
      foreach my $sample (keys %samples) {

         print OUT "\t";
         if (defined($barcodes{$sample}{$barcode})) {
            print OUT $barcodes{$sample}{$barcode};
         } else {
            print OUT "0";
         }

      }
      print OUT "\n";
   } 
}

# Now, to write the main driver!
# List of gzipped FASTQ files to process is passed in as a comma-separated list
my @fileList = split(",",<STDIN>);
print "Length of file list:  ".scalar(@fileList)."\n";
foreach my $file (@fileList) {

#   print "Reading in $file\n";
   readFile($file);

}

print "Writing out merged matrix...";
writeBarcodeMatrix("mergedBarcodeMatrix.txt");
print "Done!\n";

