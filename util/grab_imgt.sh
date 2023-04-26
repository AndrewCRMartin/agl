wget --mirror --no-parent https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/
find www.imgt.org -name 'index.html*' -exec rm -f {} \;
cp www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta ../data/human_heavy_d_dna.faa
cp www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Mus_musculus/IG/IGHD.fasta ../data/mouse_heavy_d_dna.faa
rm -rf mirror/www.imgt.org
mkdir -p mirror
mv www.imgt.org mirror
