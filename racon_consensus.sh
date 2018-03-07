fasta_nopath="${2##*/}";
mkdir $5;
python scripts/get_fastqs.py -fq $1 -dr $3/temp_alls -se "$3"/"$fasta_nopath"_reformat_out_COIpred -o $5/"$fasta_nopath"_dem_fastqs
python scripts/split_fasta_to_each.py -i $4 -o $5/"$fasta_nopath"_for_graphmap

for f in $5/"$fasta_nopath"_for_graphmap/*; 
do 
nops="${f##*/}";
echo $nops
nops="${nops%.fa*}";
echo $nops
graphmap align --max-error 0.05 -r $f -d $5/"$fasta_nopath"_dem_fastqs/"$nops".fastq -o $5/"$fasta_nopath"_dem_fastqs/"${f##*/}".sam;
racon --sam $5/"$fasta_nopath"_dem_fastqs/"$nops".fastq $5/"$fasta_nopath"_dem_fastqs/"${f##*/}".sam $f $5/"$fasta_nopath"_dem_fastqs/"${f##*/}"_racon.fasta;
cat $5/"$fasta_nopath"_dem_fastqs/*racon.fasta | sed 's/Consensus_//g' > "$5".fa;
done
