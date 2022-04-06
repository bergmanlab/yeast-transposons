#!/bin/bash
# This file aims to create TE lib for Scer, 
# documenting all source fasta seqs with stable paths or Genbank ID,
# and the coordinates for LTR/Internal

run_dir="/scratch/jc33471/create_te_lib"
mkdir -p $run_dir
script_dir=${run_dir}/script
full_dir=$run_dir/full_length
split_dir=$run_dir/split_te
mkdir -p ${full_dir} ${split_dir} ${script_dir}

conda activate create_te_lib

############# Full-length library #############
cd $full_dir

# ### tester elements
# # Canonical Ty1
# accession="M18706.1"
# esearch -db nucleotide -query ${accession} | efetch -format fasta > ${accession}.fasta
# sed -E "s/>${accession}.*/>TY1canonical_gb/" ${accession}.fasta | seqkit seq -w 0 > TY1canonical_gb.fasta

# # Ty1'
# accession="BK006936.2"
# esearch -db nucleotide -query ${accession} | efetch -format fasta > ${accession}.fasta
# seqkit subseq -w 0 --chr ${accession} -r 221037:226952 ${accession}.fasta | sed -E "s/>${accession}.*/>TY1prime_gb/" > TY1prime_gb.fasta


# # TY2_917
# accession="KT203716.1"
# esearch -db nucleotide -query ${accession} | efetch -format fasta > ${accession}.fasta
# sed -E "s/>${accession}.*/>TY2_917_gb/" ${accession}.fasta | seqkit seq -w 0 > TY2_917_gb.fasta



########## prepare ORFfinder program ##########
wget -O $script_dir/ORFfinder.gz "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz"
gunzip $script_dir/ORFfinder.gz
chmod a+x $script_dir/ORFfinder
LD_LIBRARY_PATH="/home/jc33471/miniconda/envs/create_te_lib/lib"

########## make annotation for LTRs from NCBI ##########
cd $split_dir

# ### tester elements
# # Canonical Ty1_gb
# printf "##gff-version 3\n" > TY1canonical_gb.gff
# # annotation ref https://www.ncbi.nlm.nih.gov/pmc/articles/PMC363300/
# printf "TY1canonical_gb\tGenbank\tLTR_retrotransposon\t1\t5918\t.\t+\t1\tName=M18706.1;organism=Saccharomyces cerevisiae;strain=JB84A\n" >> TY1canonical_gb.gff
# printf "TY1canonical_gb\tGenbank\tlong_terminal_repeat\t1\t334\t.\t+\t1\tName=5' long terminal repeat\n" >> TY1canonical_gb.gff
# printf "TY1canonical_gb\tGenbank\tlong_terminal_repeat\t5585\t5918\t.\t+\t1\tName=3' long terminal repeat\n" >> TY1canonical_gb.gff
# # with ORFfinder, start with ATG only, stop codon included
# printf "TY1canonical_gb\tORFfinder\tCDS\t294\t1616\t.\t+\t1\tName=GAG\n" >> TY1canonical_gb.gff
# # with ORFfinder, start with any sense codons, with stop codon included
# printf "TY1canonical_gb\tORFfinder\tCDS\t1576\t5562\t.\t+\t1\tName=POL\n" >> TY1canonical_gb.gff
# gt gff3 -sort TY1canonical_gb.gff > TY1canonical_gb.sorted.gff


# # Ty1prime_gb
# printf "##gff-version 3\n" > TY1prime_gb.gff
# # no publication found, annotation from https://www.yeastgenome.org/locus/S000006808
# printf "TY1prime_gb\tSGD\tLTR_retrotransposon\t1\t5916\t.\t+\t1\tID=SGD:S000006808;Name=YBLWTy1-1;organism=Saccharomyces cerevisiae;strain=S288C\n" >> TY1prime_gb.gff
# printf "TY1prime_gb\tSGD\tlong_terminal_repeat\t1\t334\t.\t+\t1\tName=5' long terminal repeat\n" >> TY1prime_gb.gff
# printf "TY1prime_gb\tSGD\tlong_terminal_repeat\t5585\t5916\t.\t+\t1\tName=3' long terminal repeat\n" >> TY1prime_gb.gff
# # with ORFfinder, start with ATG only, stop codon included
# printf "TY1prime_gb\tORFfinder\tCDS\t294\t1616\t.\t+\t1\tName=GAG\n" >> TY1prime_gb.gff
# # with ORFfinder, start with any sense codons, with stop codon included
# printf "TY1prime_gb\tORFfinder\tCDS\t1576\t5562\t.\t+\t1\tName=POL\n" >> TY1prime_gb.gff
# gt gff3 -sort TY1prime_gb.gff > TY1prime_gb.sorted.gff


# # TY2_917_gb
# printf "##gff-version 3\n" > TY2_917_gb.gff
# # frameshift confirmed https://pubmed.ncbi.nlm.nih.gov/2842793/ , annotation of LTRs from Genbank record
# # and correct 5' LTR from 333 to 332
# printf "TY2_917_gb\tGenbank\tLTR_retrotransposon\t1\t5960\t.\t+\t1\tID=KT203716.1;Name=KT203716.1;organism=Saccharomyces cerevisiae;strain=S288C\n" >> TY2_917_gb.gff
# printf "TY2_917_gb\tGenbank\tlong_terminal_repeat\t1\t332\t.\t+\t1\tName=5' long terminal repeat\n" >> TY2_917_gb.gff
# printf "TY2_917_gb\tGenbank\tlong_terminal_repeat\t5628\t5960\t.\t+\t1\tName=3' long terminal repeat\n" >> TY2_917_gb.gff
# # with ORFfinder, start with ATG only, stop codon included
# printf "TY2_917_gb\tORFfinder\tCDS\t292\t1608\t.\t+\t1\tName=GAG\n" >> TY2_917_gb.gff
# # with ORFfinder, start with any sense codons, with stop codon included
# printf "TY2_917_gb\tORFfinder\tCDS\t1562\t5605\t.\t+\t1\tName=POL\n" >> TY2_917_gb.gff
# gt gff3 -sort TY2_917_gb.gff > TY2_917_gb.sorted.gff

### Ty query sequences identified using phylogenetic method https://github.com/bergmanlab/jingxuan/issues/24#issuecomment-1059297528
cd $split_dir
export LD_LIBRARY_PATH="/home/jc33471/miniconda/envs/create_te_lib/lib:$LD_LIBRARY_PATH"
anno_dir="/scratch/cbergman/yeast_long_read_asms_12-24-2021/"

create_annotation(){
    ltrharvest_file=${anno_dir}/${dataset}/${strain}/ltrharvest/${strain}.ltrdigest.renamed.gff
    # grep coordinates for LTRs from LTRharvest
    ID=`awk -v accession="${accession}" -v start="${start}" -v end="${end}" 'BEGIN{FS="\t";OFS="\t";}{
        if($1==accession && $3=="LTR_retrotransposon" && $4==start && $5==end) {
            split($9,attr,";");
            sub(/ID=/,"",attr[1]); print attr[1]
        }
        }' ${ltrharvest_file}`
    if [[ -z $ID ]]; then
        echo "No same prediction in LTRharvest."
        return 1
    fi
    if [[ ${strand} == "+" ]]; then
        length_5=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$4==start){print $5-$4+1}}'`
        length_3=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$5==end){print $5-$4+1}}'`
    elif [[ ${strand} == "-" ]]; then
        length_3=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$4==start){print $5-$4+1}}'`
        length_5=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$5==end){print $5-$4+1}}'`   
    fi
    printf "##gff-version 3\n" > ${TE}.gff
    printf "${TE}\tLTRharvest\tLTR_retrotransposon\t1\t$((${end}-${start}+1))\t.\t+\t1\tName=${TE};organism=${organism};strain=${strain};Note=${accession}:${start}-${end}\n" >> ${TE}.gff
    printf "${TE}\tLTRharvest\tlong_terminal_repeat\t1\t${length_5}\t.\t+\t1\tName=5' long terminal repeat\n" >> ${TE}.gff
    printf "${TE}\tLTRharvest\tlong_terminal_repeat\t$((${end}-${start}+1-${length_3}+1))\t$((${end}-${start}+1))\t.\t+\t1\tName=3' long terminal repeat\n" >> ${TE}.gff
    # ORFfinder results are shown in a wired coordinating system, neither 0-based nor 1-based.
    # need to +1 to both start and end to calibrate to 1-based system
    ${script_dir}/ORFfinder -in ${full_dir}/${TE}.fasta -s 0 -ml 600 -out $split_dir/${TE}.orf1 -outfmt 3
    gag_start=`grep "CDS" $split_dir/${TE}.orf1 | sort -n -k1 | awk 'NR==1' | cut -f1`
    gag_end=`grep "CDS" $split_dir/${TE}.orf1 | sort -n -k1 | awk 'NR==1' | cut -f2`
    ${script_dir}/ORFfinder -in ${full_dir}/${TE}.fasta -s 2 -ml 600 -out $split_dir/${TE}.orf2 -outfmt 3
    pol_start=`grep "CDS" $split_dir/${TE}.orf2 | sort -n -k1 | awk 'NR==2' | cut -f1`
    pol_end=`grep "CDS" $split_dir/${TE}.orf2 | sort -n -k1 | awk 'NR==2' | cut -f2`
    # with ORFfinder, start with ATG only, stop codon included
    printf "${TE}\tORFfinder\tCDS\t$((${gag_start}+1))\t$((${gag_end}+1))\t.\t+\t1\tName=GAG\n" >> ${TE}.gff
    # with ORFfinder, start with any sense codons, with stop codon included
    printf "${TE}\tORFfinder\tCDS\t$((${pol_start}+1))\t$((${pol_end}+1))\t.\t+\t1\tName=POL\n" >> ${TE}.gff
    gt gff3 -sort ${TE}.gff > ${TE}.sorted.gff
}

print_ltr_length(){
    ltrharvest_file=${anno_dir}/${dataset}/${strain}/ltrharvest/${strain}.ltrdigest.renamed.gff
    # grep coordinates for LTRs from LTRharvest
    ID=`awk -v accession="${accession}" -v start="${start}" -v end="${end}" 'BEGIN{FS="\t";OFS="\t";}{
        if($1==accession && $3=="LTR_retrotransposon" && $4==start && $5==end) {
            split($9,attr,";");
            sub(/ID=/,"",attr[1]); print attr[1]
        }
        }' ${ltrharvest_file}`
    if [[ ${strand} == "+" ]]; then
        length_5=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$4==start){print $5-$4+1}}'`
        length_3=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$5==end){print $5-$4+1}}'`
    elif [[ ${strand} == "-" ]]; then
        length_3=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$4==start){print $5-$4+1}}'`
        length_5=`grep -E "${ID}(;|$)" ${ltrharvest_file} | awk -v accession=${accession} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";}{if($3=="long_terminal_repeat"&&$5==end){print $5-$4+1}}'`   
    fi
    echo $length_5 $length_3
}

# Ty1-canonical CZA_DBVPG6044_f590+_TY1|CABIKB010000018.1|42993-48917
# modify for each elements
dataset="czaja";accession="CABIKB010000018.1";
start=42994;end=48917;strand="+";
TE="Ty1-canonical";organism="Saccharomyces cerevisiae";strain="DBVPG6044"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f590+_TY1"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same

# Ty1-prime CZA_Y12_f294+_TY1|CABIJY010000010.1|606514-612428 # same
# modify for each elements
dataset="czaja";accession="CABIJY010000010.1";
start=606515;end=612428;strand="+";
TE="Ty1-prime";organism="Saccharomyces cerevisiae";strain="Y12"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f294+_TY1"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length 

# # TY2_917 CZA_S288c_f462+_TY2|CABIJX010000016.1|90664-96623
# # modify for each elements
# dataset="czaja";accession="CABIJX010000016.1";
# start=90665;end=96623;strand="+";
# TE="TY2_917";organism="Saccharomyces cerevisiae";strain="S288c"
# if [[ ${strand} == "+" ]]; then
#     esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
#     seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
# elif [[ ${strand} == "-" ]]; then
#     esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
#     seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
# fi
# create_annotation
# # check repeatmasker
# element_id="f462+_TY2"
# grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
# print_ltr_length # same

# Ty1-mosaic CZA_S288c_f87+_TY1|CABIJX010000002.1|536094-542020
# modify for each elements
dataset="czaja";accession="CABIJX010000002.1";
start=536095;end=542020;strand="+";
TE="Ty1-mosaic";organism="Saccharomyces cerevisiae";strain="S288c"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f87+_TY1"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same 

# Ty1p-ow yue_CBS432_f139+_TY1|CP020245.1|1300770-1306656
dataset="yue";accession="CP020245.1";
start=1300771;end=1306656;strand="+";
TE="Ty1p-ow";organism="Saccharomyces paradoxus";strain="CBS432"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f139+_TY1"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same

# Ty1p-nw yue_UFRJ50816_f400+_TY1|CP020300.1|115626-121510
dataset="yue";accession="CP020300.1";
start=115627;end=121510;strand="+";
TE="Ty1p-nw";organism="Saccharomyces paradoxus";strain="UFRJ50816"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f400+_TY1"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same

# Ty2 CZA_S288c_f339-_TY2|CABIJX010000011.1|138075-144034
dataset="czaja";accession="CABIJX010000011.1";
start=138076;end=144034;strand="-";
TE="Ty2";organism="Saccharomyces cerevisiae";strain="S288c"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f339-_TY2"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same


# Ty3 CZA_S288c_f95+_TY3|CABIJX010000002.1|707538-712889
dataset="czaja";accession="CABIJX010000002.1";
start=707539;end=712889;strand="+";
TE="Ty3";organism="Saccharomyces cerevisiae";strain="S288c"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f95+_TY3"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same

# Ty3p-ow yue_CBS432_f326+_TY3_1p|CP020250.1|169726-175086
dataset="yue";accession="CP020250.1";
start=169727;end=175086;strand="+";
TE="Ty3p-ow";organism="Saccharomyces paradoxus";strain="CBS432"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f326+_TY3_1p"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same


# Ty4 CZA_Y12_f354+_TY4|CABIJY010000012.1|501183-507406
dataset="czaja";accession="CABIJY010000012.1";
start=501184;end=507406;strand="+";
TE="Ty4";organism="Saccharomyces cerevisiae";strain="Y12"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f354+_TY4"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same

# Ty4p-ow yue_N44_f223-_TY4|CP020264.1|288028-294302
dataset="yue";accession="CP020264.1";
start=288029;end=294302;strand="-";
TE="Ty4p-ow";organism="Saccharomyces paradoxus";strain="N44"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f223-_TY4"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same, repeatmasker 371,1909,3615,371?

# Tsu4p-nw yue_UFRJ50816_f45-_TSU4|CP020294.1|554569-560566
dataset="yue";accession="CP020294.1";
start=554570;end=560566;strand="-";
TE="Tsu4p-nw";organism="Saccharomyces paradoxus";strain="UFRJ50816"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
create_annotation
# check repeatmasker
element_id="f45-_TSU4"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same

# Ty5p-ow yue_CBS432_f409+_TY5|CP020252.1|406512-411888
dataset="yue";accession="CP020252.1";
start=406513;end=411888;strand="+";
TE="Ty5p-ow";organism="Saccharomyces paradoxus";strain="CBS432"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
# create_annotation
ltrharvest_file=${anno_dir}/${dataset}/${strain}/ltrharvest/${strain}.ltrdigest.renamed.gff
# grep coordinates for LTRs from LTRharvest
ID=`awk -v accession="${accession}" -v start="${start}" -v end="${end}" 'BEGIN{FS="\t";OFS="\t";}{
    if($1==accession && $3=="LTR_retrotransposon" && $4==start && $5==end) {
        split($9,attr,";");
        sub(/ID=/,"",attr[1]); print attr[1]
    }
    }' ${ltrharvest_file}`
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
printf "##gff-version 3\n" > ${TE}.gff
printf "${TE}\tLTRharvest\tLTR_retrotransposon\t1\t$((${end}-${start}+1))\t.\t+\t1\tName=${TE};organism=${organism};strain=${strain};Note=${accession}:${start}-${end}\n" >> ${TE}.gff
printf "${TE}\tLTRharvest\tlong_terminal_repeat\t1\t${length_5}\t.\t+\t1\tName=5' long terminal repeat\n" >> ${TE}.gff
printf "${TE}\tLTRharvest\tlong_terminal_repeat\t$((${end}-${start}+1-${length_3}+1))\t$((${end}-${start}+1))\t.\t+\t1\tName=3' long terminal repeat\n" >> ${TE}.gff
# ORFfinder
${script_dir}/ORFfinder -in ${full_dir}/${TE}.fasta -s 0 -ml 600 -out $split_dir/${TE}.orf1 -outfmt 3
gag_start=`grep "CDS" $split_dir/${TE}.orf1 | sort -n -k1 | awk 'NR==1' | cut -f1`
gag_end=`grep "CDS" $split_dir/${TE}.orf1 | sort -n -k1 | awk 'NR==1' | cut -f2`
# with ORFfinder, start with ATG only, stop codon included
printf "${TE}\tORFfinder\tCDS\t$((${gag_start}+1))\t$((${gag_end}+1))\t.\t+\t1\tName=GAG-POL\n" >> ${TE}.gff
gt gff3 -sort ${TE}.gff > ${TE}.sorted.gff
# check repeatmasker
element_id="f409+_TY5"
grep "${element_id}" ${anno_dir}/${dataset}/${strain}/repeatmasker/${strain}.bed | cut -f11
print_ltr_length # same


# Ty5 bnd_EM14S01-3B_f230-_TY5|LR813555.1|572815-578164 # ltr harvest not working, treat differently
dataset="bendixsen";accession="LR813555.1";
start=572816;end=578164;strand="-";
TE="Ty5";organism="Saccharomyces cerevisiae";strain="EM14S01-3B"
if [[ ${strand} == "+" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
elif [[ ${strand} == "-" ]]; then
    esearch -db nucleotide -query ${accession} | efetch -format fasta > ${full_dir}/${accession}.fasta
    seqkit subseq -w 0 --chr ${accession} -r ${start}:${end} ${full_dir}/${accession}.fasta | seqkit seq -r -p -w 0 | sed -E "s/>${accession}.*/>${TE}/" > ${full_dir}/${TE}.fasta
fi
# create_annotation based on repeatmasker LTR annotation
length_5=251
length_3=251
printf "##gff-version 3\n" > ${TE}.gff
printf "${TE}\tLTRharvest\tLTR_retrotransposon\t1\t$((${end}-${start}+1))\t.\t+\t1\tName=${TE};organism=${organism};strain=${strain};Note=${accession}:${start}-${end}\n" >> ${TE}.gff
printf "${TE}\tLTRharvest\tlong_terminal_repeat\t1\t${length_5}\t.\t+\t1\tName=5' long terminal repeat\n" >> ${TE}.gff
printf "${TE}\tLTRharvest\tlong_terminal_repeat\t$((${end}-${start}+1-${length_3}+1))\t$((${end}-${start}+1))\t.\t+\t1\tName=3' long terminal repeat\n" >> ${TE}.gff
# ORFfinder
${script_dir}/ORFfinder -in ${full_dir}/${TE}.fasta -s 0 -ml 600 -out $split_dir/${TE}.orf1 -outfmt 3
gag_start=`grep "CDS" $split_dir/${TE}.orf1 | sort -n -k1 | awk 'NR==1' | cut -f1`
gag_end=`grep "CDS" $split_dir/${TE}.orf1 | sort -n -k1 | awk 'NR==1' | cut -f2`
# with ORFfinder, start with ATG only, stop codon included
printf "${TE}\tORFfinder\tCDS\t$((${gag_start}+1))\t$((${gag_end}+1))\t.\t+\t1\tName=GAG-POL\n" >> ${TE}.gff
gt gff3 -sort ${TE}.gff > ${TE}.sorted.gff

# combine all seqs
# cat ${full_dir}/TY1canonical.fasta ${full_dir}/TY1prime.fasta ${full_dir}/TY1.fasta ${full_dir}/TY1p_ow.fasta ${full_dir}/TY1p_nw.fasta ${full_dir}/TY2.fasta ${full_dir}/TY3.fasta ${full_dir}/TY3_1p.fasta ${full_dir}/TY4.fasta ${full_dir}/TY4p_ow.fasta ${full_dir}/TSU4.fasta ${full_dir}/TY5p.fasta ${full_dir}/TY5c.fasta | seqkit seq -w 0 --upper-case > ${full_dir}/ty_elements_from_genome.fasta
cat ${full_dir}/Ty*.fasta ${full_dir}/Tsu*.fasta | seqkit seq -w 0 --upper-case > ${full_dir}/ty_elements_from_genome.fasta
# combine all annotations
# gt gff3 -sort -tidy -addids -retainids -checkids ${split_dir}/TY1canonical.sorted.gff ${split_dir}/TY1prime.sorted.gff ${split_dir}/TY1.sorted.gff ${split_dir}/TY1p_ow.sorted.gff ${split_dir}/TY1p_nw.sorted.gff ${split_dir}/TY2.sorted.gff ${split_dir}/TY3.sorted.gff ${split_dir}/TY3_1p.sorted.gff ${split_dir}/TY4.sorted.gff ${split_dir}/TY4p_ow.sorted.gff ${split_dir}/TSU4.sorted.gff ${split_dir}/TY5p.sorted.gff ${split_dir}/TY5c.sorted.gff > ${split_dir}/full_length_from_genome.sorted.gff
gt gff3 -sort -tidy -addids -retainids -checkids ${split_dir}/Ty*.sorted.gff ${split_dir}/Tsu*.sorted.gff > ${split_dir}/full_length_from_genome.sorted.gff
gt gff3validator ${split_dir}/full_length_from_genome.sorted.gff

### split LTR and internal
cp /home/jc33471/jingxuan/src/python/split_cns_te_lib_yeast.py $script_dir/split_cns_te_lib_yeast.py
python ${script_dir}/split_cns_te_lib_yeast.py $full_dir/ty_elements_from_genome.fasta $split_dir/full_length_from_genome.sorted.gff | seqkit sort --by-name -w 0 > $split_dir/split_full_from_genome.fasta

