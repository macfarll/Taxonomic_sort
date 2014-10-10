#!/bin/bash
while test 1 -gt 0; do
        case "$1" in
                -h)
                        echo "options:"
                        echo "-h,       Show brief help."
                        echo "-f,      Input fasta."
                        echo "-dir,      Working directory, defualts to pwd."
                        echo "-nproc,      Number of processors to be used in rapsearch."
			echo "-tax_name,	Taxonomic name to be searched for, must be spelled properly with first letter capitalized (ex: Cyanobacteria)."
			echo "-tax_level,	Taxonomic level to be searched at, must be spelled properly, all lower case (ex: phylum)."
			echo "-s,	Turns off summary statistics."
                        exit 0
                        ;;
                -f)
                        shift
                        if test $# -gt 0; then
                                export FASTA=$1
				echo "Fasta to be used:"
				echo " "$FASTA
                        else
                                echo "Fasta not provided."
                                exit 1
                        fi
                        shift
                        ;;

                -tax_name)
                        shift
                        if test $# -gt 0; then
                                export LABEL=$1
				echo "\nTaxonomic label to be used:"
				echo " "$LABEL
                        else
				echo "Taxonomic name to search for not specified"
				exit 1
			fi
			shift
			;;
                -tax_level)
                        shift
                        if test $# -gt 0; then
                                export LEVEL=$1
				echo "\nTaxonomic level to search at:"
				echo " "$LEVEL
                        else
                                echo "Taxonomic level to search at not specified."
				exit 1
                        fi
                        shift
                        ;;

                -dir)
                        shift
                        if test $# -gt 0; then
                                export DIR=$1
				echo "Directory to be used:"
				echo " "$DIR
                        else
                                echo "No directory specified."
				exit 1
                        fi
                        shift
                        ;;
                -nproc)
                        shift
                        if test $# -gt 0; then
                                export NPROC=$1
				echo "\nNumber of processors to be used in rapsearch:"
				echo " "$NPROC
                        fi
                        shift
                        ;;
                -s)
                        shift
                        if test 1 -gt 0; then
                                export SUMMARY="1"
			else
				break
                        fi
                        ;;
                *)
                        break
                        ;;
        esac
done
if [ -z $NPROC ]; then
        echo "\nNumber of processors not specified, using 1."
	NPROC="1"
        break
fi
if [ -z $FASTA ]; then
        echo "Fasta not specified, exiting."
        exit
fi
if [ -z $LEVEL ]; then
        echo "Taxonomic level not specified, exiting."
        exit
fi
if [ -z $DIR ]; then
        echo "\nDirectory not specified, using pwd:"
	echo " "`pwd`
        DIR=`pwd`
fi
if [ -z $LABEL ]; then
        echo "Target taxonomic label not specified, exiting."
        exit
fi
#fasta2single_line.pl ${FASTA} > ${DIR}/1l_fasta_temp.fa
#FASTA=${DIR}/1l_fasta_temp.fa
#grep ">" ${FASTA} | awk '{printf ">%d\n", NR}' > ${DIR}/headers
#grep -v ">" ${FASTA} > ${DIR}/seqs
#paste -d"\n" headers seqs > ${DIR}/renumbered.fa
#rm ${DIR}/seqs ${DIR}/headers
#echo running prodigal
#prodigal -m -q -a ${DIR}/prodigal.faa -i ${DIR}/renumbered.fa -o /dev/null
#sed 's/\# [[:print:]]*$//g' prodigal.faa > ${DIR}/formatted_prodigal.faa
#rapsearch -q ${DIR}/formatted_prodigal.faa -d /usr/local/Programs/ncbi-blast-2.2.28+/db/rap2.16_nr -o ${DIR}/rapsearch_out -z ${NPROC} -v 1 -b 1 -t a -p t -a t
grep "^[0-9]" rapsearch_out.aln | sed 's/\t/|/g' | sed 's/ vs gi//g' | awk -F"|" '{print $1"\t"$2}' > ${DIR}/gi_mapping
sed 's/_/\t/g' ${DIR}/gi_mapping > ${DIR}/3col_gi_mapping
grep -v ">" ${DIR}/rapsearch_out.m8 | awk -F"|" '{print $2}' | sed -e '1,5d' > ${DIR}/gi.txt
#ruby etaxa.rb ${DIR}/gi.txt >/dev/null 2>/dev/null
#sh etaxa.sh >/dev/null 2>/dev/null
grep -A1 "^[0-9]" ${DIR}/taxa.txt | sed '/[0-9]/ a\HOLDER' | sed ':a;N;$!ba;s/\nHOLDER\n/;/g' | sed 's/ //g' | sed -e 's/;/-/7' > ${DIR}/gi_taxa.txt
echo "#!/usr/bin/Rscript" > created_r_merger
echo "f <- file(\"stdin\")" >> created_r_merger
echo "open(f)" >> created_r_merger
echo "while(length(line <- readLines(f,n=1)) > 0) {" >> created_r_merger
echo "  write(line, stderr())" >> created_r_merger
echo "  GI<-paste(line[1],\"/gi_taxa.txt\", sep=\"\")" >> created_r_merger
echo "  MAPPING<-paste(line[1],\"/3col_gi_mapping\", sep=\"\")" >> created_r_merger
echo "  OUT<-paste(line[1],\"/merged_r_output\", sep=\"\")" >> created_r_merger
echo "}" >> created_r_merger
echo "print(\"Merging by gi in r\")" >> created_r_merger
echo "print(GI)" >> created_r_merger
echo "lineage <- read.csv2(file = {GI}, header = F, sep = \";\")" >> created_r_merger
echo "lineage <- lineage[-7]" >> created_r_merger
echo "lineage <- lineage[-8]" >> created_r_merger
echo "lineage <- lineage[-9] " >> created_r_merger
echo "colnames(lineage) <- c(\"gi\", \"kingdom\", \"phylum\", \"class\", \"order\", \"family\")" >> created_r_merger
echo "mapping <- read.csv2(file = {MAPPING}, header = F, sep =\"\\\t\")" >> created_r_merger
echo "colnames(mapping) <- c(\"contig\", \"gene\", \"gi\")" >> created_r_merger
echo "merged <- merge(x = lineage, y = mapping, by = \"gi\", all.y = T)" >> created_r_merger
echo "write.csv(merged, quote = F,file = {OUT})" >> created_r_merger
echo ${DIR} | Rscript ${DIR}/created_r_merger
# >/dev/null 2>/dev/null
#rm created_r_merger
sed 's/^[0-9]*[0-9],//g' ${DIR}/merged_r_output | sed -e '1s/^.//' | awk '!a[$0]++' | sed 's/ /_/g' > ${DIR}/no_duplicate_csv
echo ${DIR}"\n"${LEVEL} > tax_search_param
echo "#!/usr/bin/Rscript" > created_r_mode
echo "f <- file(\"stdin\")" >> created_r_mode
echo "open(f)" >> created_r_mode
echo "while(length(line <- readLines(f,n=2)) > 0) {" >> created_r_mode
echo "  write(line, stderr())" >> created_r_mode
echo "  CSV<-paste(line[1],\"/no_duplicate_csv\", sep=\"\")" >> created_r_mode
echo "  TABLE<-paste(line[1],\"/table_of_modes\",sep=\"\")" >> created_r_mode
echo "  LEVEL<-print(line[2])" >> created_r_mode
echo "}" >> created_r_mode
echo "Mode <- function(x) {" >> created_r_mode
echo "  ux <- unique(x)" >> created_r_mode
echo "  ux[which.max(tabulate(match(x, ux)))]" >> created_r_mode
echo "}" >> created_r_mode
echo "data <- read.csv2(file ={CSV}, header = T, sep = \",\")" >> created_r_mode
echo "Table <- sapply(1:max(as.numeric(unique(data\$contig))), FUN = function(x) {Mode(data[data\$contig == x,{LEVEL}])})" >> created_r_mode
echo "write.csv(Table, file = {TABLE})" >> created_r_mode
cat tax_search_param | Rscript ${DIR}/created_r_mode
#>/dev/null 2>/dev/null
rm created_r_mode tax_search_param
grep "${LABEL}" ${DIR}/table_of_modes | sed -e "s/$LABEL//" | sed 's/"//g'| sed ':a;N;$!ba;s/\n/\$\n\>/g' | sed 's/,//g' | sed '1s/^/>/' | sed "\$a$" | sed ':a;N;$!ba;s/\n\$/$/g' > ${DIR}/line_sep_grep_commands
while read p; do
  grep -A1 "$p" ${DIR}/renumbered.fa >> ${DIR}/target_contigs.fa
done <${DIR}/line_sep_grep_commands
rm line_sep_grep_commands 1l_fasta_temp.fa
echo "\nThe output file, 'target_contigs.fa' has been generated.\n"
if [ -z $SUMMARY ]; then
	echo "#!/usr/bin/Rscript" > created_r_summary
	echo "f <- file(\"stdin\")" >> created_r_summary
	echo "open(f)" >> created_r_summary
	echo "while(length(line <- readLines(f,n=1)) > 0) {" >> created_r_summary
	echo "  FASTA<-paste(line[1],\"/target_contigs.fa\", sep=\"\")" >> created_r_summary
	echo "}" >> created_r_summary
	echo "suppressMessages(library(Biostrings))" >> created_r_summary
	echo "target <- readDNAStringSet({FASTA})" >> created_r_summary
	echo "print(\"Total length of assembled contigs:\")" >> created_r_summary
	echo "paste(\" \", sum(width(target)), sep =\"\")" >> created_r_summary
	echo "print(\" \")" >> created_r_summary
	echo "print(\"Number of target contigs:\")" >> created_r_summary
	echo "paste(\" \", length(target), sep=\"\")" >> created_r_summary
	echo "print(\" \")" >> created_r_summary
	echo "print(\"Average target contig length:\")" >> created_r_summary
	echo "paste(\" \", sum(width(target))/length(target), sep=\"\")" >> created_r_summary
	echo ${DIR} | Rscript ${DIR}/created_r_summary > ${DIR}/CLEANER
	cat ${DIR}/CLEANER | sed 's/\[1\]//g' | sed 's/"//g' | sed 's/^ //g'
	rm ${DIR}/CLEANER created_r_summary
	break
fi
