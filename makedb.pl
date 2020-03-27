#!/usr/bin/perl
use strict;
use lib "./";
use fasta;

# tcag
%::dna2aa = ('ttt' => 'F',
             'ttc' => 'F',
             'tta' => 'L',
             'ttg' => 'L',
             'tct' => 'S',
             'tcc' => 'S',
             'tca' => 'S',
             'tcg' => 'S',
             'tat' => 'Y',
             'tac' => 'Y',
             'taa' => '*',
             'tag' => '*',
             'tgt' => 'C',
             'tgc' => 'C',
             'tga' => '*',
             'tgg' => 'W',
             
             'ctt' => 'L',
             'ctc' => 'L',
             'cta' => 'L',
             'ctg' => 'L',
             'cct' => 'P',
             'ccc' => 'P',
             'cca' => 'P',
             'ccg' => 'P',
             'cat' => 'H',
             'cac' => 'H',
             'caa' => 'Q',
             'cag' => 'Q',
             'cgt' => 'R',
             'cgc' => 'R',
             'cga' => 'R',
             'cgg' => 'R',
           
             'att' => 'I',
             'atc' => 'I',
             'ata' => 'I',
             'atg' => 'M',
             'act' => 'T',
             'acc' => 'T',
             'aca' => 'T',
             'acg' => 'T',
             'aat' => 'N',
             'aac' => 'N',
             'aaa' => 'K',
             'aag' => 'K',
             'agt' => 'S',
             'agc' => 'S',
             'aga' => 'R',
             'agg' => 'R',
             
             'gtt' => 'V',
             'gtc' => 'V',
             'gta' => 'V',
             'gtg' => 'V',
             'gct' => 'A',
             'gcc' => 'A',
             'gca' => 'A',
             'gcg' => 'A',
             'gat' => 'D',
             'gac' => 'D',
             'gaa' => 'E',
             'gag' => 'E',
             'ggt' => 'G',
             'ggc' => 'G',
             'gga' => 'G',
             'ggg' => 'G');

my $dataDir = "./data";
my $outDir  = "./mydata";

`mkdir -p $outDir`;

my @files = GetFileList($dataDir);
ProcessFiles($outDir, 'heavy_v', @files);
ProcessFiles($outDir, 'heavy_d', @files);
ProcessFiles($outDir, 'heavy_j', @files);
ProcessFiles($outDir, 'heavy_c', @files);

ProcessFiles($outDir, 'light_v', @files);
ProcessFiles($outDir, 'light_j', @files);
ProcessFiles($outDir, 'light_c', @files);

sub ProcessFiles
{
    my ($outDir, $type, @files) = @_;
    my $outFilename = "${outDir}/${type}.dat";
    if(open(my $outFp, '>', $outFilename))
    {
        foreach my $file (@files)
        {
            if($file =~ /$type/)
            {
                ProcessAFile($outFp, $file)
            }
        }
        close $outFp;
    }
    else
    {
        print STDERR "Can't write $outFilename\n";
    }
}

sub ProcessAFile
{
    my($outFp, $file) = @_;
    if(open(my $inFp, '<', $file))
    {
        my($id, $info, $sequence);
        while(1)
        {
            ($id, $info, $sequence) = fasta::ReadFasta($inFp);
            last if ($id eq '');
            PrintTranslations($outFp, $info, $sequence);
        }
              
        close $inFp;
    }
    else
    {
        print STDERR "Can't read $file\n";
    }
}

sub PrintTranslations
{
    my($outFp, $info, $sequence) = @_;


    for(my $frame=0; $frame < 3; $frame++)
    {
        my $header = CreateHeader($info, $frame);
        my $aaSeq = Translate($sequence, $frame);
        print $outFp "$header\n$aaSeq\n";
    }
}

sub CreateHeader
{
    my($info, $frame) = @_;

    my @fields = split(/\|/, $info);
    my $species = $fields[2];
    $species =~ s/\_.*//;  # Remove strain
    my $gene    = $fields[1];
    my $domain  = $fields[4];
    my $header = ">${domain}_${gene}_F${frame}_$species";
    return($header);
}

sub Translate
{
    my ($dnaSeq, $frame) = @_;
    my $aaSeq = '';
    
    $dnaSeq = substr($dnaSeq, $frame);
    $dnaSeq = "\L$dnaSeq";

    while(length($dnaSeq))
    {
        my $codon  = substr($dnaSeq, 0, 3);
        $dnaSeq    = substr($dnaSeq, 3);
        $aaSeq    .= $::dna2aa{$codon};
    }

    return($aaSeq);
}


sub GetFileList
{
    my($dataDir) = @_;
    my @files = ();
    
    if(opendir(my $dp, $dataDir))
    {
        @files = readdir($dp);
        @files = grep (/\.faa/, @files);
        
        closedir($dp);
    }

    foreach my $file (@files)
    {
        $file = "$dataDir/$file";
    }
    
    return @files;
}
    
