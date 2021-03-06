#!/bin/bash
while test 1 -gt 0; do
        case "$1" in
                -h)
                        echo "options:"
                        echo "-h,       Show brief help."
                        echo "-f,      Input fasta."
			echo "-d	Data prefix to be prepended to outputs."
			echo "-r	Start from rapsearch out, must be followed by <rapsearch_out>.aln file. Reduces runtime if program has been previously run."
			echo "-q	Quickly re-run the program with new parameters, requires the script have been ran previously in this directory with this data."
                        echo "-dir,      Working directory, defualts to pwd."
                        echo "-nproc,      Number of processors to be used in rapsearch."
			echo "-tax_name,	Taxonomic name to be searched for, spelling sensitive. (ex: Cyanobacteria)"
			echo "-tax_level,	Taxonomic level to be searched at, spelling sensitive. (ex: phylum)"
			echo "-s,	Turns off summary statistics."
                        exit 0
                        ;;
		-r)
			shift
			if test $# -gt 0; then
				export RAPSEARCHOUT=$1
				echo "\nUsing rapsearch output:"
				echo " "$RAPSEARCHOUT
			fi
			shift
			;;
                -q)
                        shift
                        if test 1 -gt 0; then
                                export QUICKPARAM="1"
                                echo "\nUsing quick option."
                        fi
                        ;;
                -f)
                        shift
                        if test $# -gt 0; then
                                export FASTA=$1
				echo "\nFasta to be used:"
				echo " "$FASTA
                        else
                                echo "Fasta not provided."
                                exit 1
                        fi
                        shift
                        ;;
                -d)
                        shift
                        if test $# -gt 0; then
                                export DATA=$1
                                echo "\nOutput will have data prefix:"
                                echo " "$DATA
                        else
                                echo "\nPlease provide a dataname."
                                exit
                        fi
                        shift
                        ;;

                -tax_name)
                        shift
                        if test $# -gt 0; then
                                export LABEL=`echo $1 | tr [A-Z] [a-z] | sed -e 's/^./\U&/g; s/ ./\U&/g'`
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
                                export LEVEL=`echo $1 | tr [A-Z] [a-z] | sed -e 's/^./\U&/g; s/ ./\U&/g'`
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
				echo "\nDirectory to be used:"
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
rm created_r* 1l_fasta_temp.fa 3col_gi_mapping formatted_prodigal.faa gi_mapping gi_taxa.txt gi.txt no_duplicate_csv merged_r_output line_sep_grep_commands target_contigs.fa Esearcher_subscript.sh 2>/dev/null
if [ -z $NPROC ]; then
        echo "\nNumber of processors not specified, using 1."
	NPROC="1"
fi
if [ -z $DATA ]; then
	echo "\nData prefix not provided, please provide a data prefix."
	exit
fi
if [ -n "$QUICKPARAM" ]; then
	if [ -s "${DATA}_taxa.txt" ]; then
		echo "\nUsing previous output: "${DATA}"_taxa.txt"
	else
		echo "Program has not been run on this data before in this directory, "${DATA}"_taxa.txt is missing."
	fi
elif [ -n "$RAPSEARCHOUT" ]; then
        echo "\nUsing rapsearch output: "${RAPSEARCHOUT}
elif [ -n "$FASTA" ]; then
        echo "\nUsing fasta:"
	echo " "${FASTA}
else
	echo "No input provided, exiting."
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
if [ -z $RAPSEARCHOUT ]; then
	if [ -z $QUICKPARAM ]; then
		cat ${FASTA} | awk '{if (substr($0,1,1)==">"){if (p) {print "\n";} print $0} else printf("%s",$0);p++;}END{print"\n"}' | sed '/^\s*$/d' > ${DIR}/1l_fasta_temp.fa
		FASTA=${DIR}/1l_fasta_temp.fa
		grep ">" ${FASTA} | awk '{printf ">%d\n", NR}' > ${DIR}/headers
		grep -v ">" ${FASTA} > ${DIR}/seqs
		paste -d"\n" headers seqs > ${DIR}/${DATA}_renumbered.fa
		rm ${DIR}/seqs ${DIR}/headers
		echo "running prodigal"
		prodigal -m -q -a ${DIR}/${DATA}_prodigal.faa -i ${DIR}/${DATA}_renumbered.fa -o /dev/null
		sed 's/\# [[:print:]]*$//g' ${DATA}_prodigal.faa > ${DIR}/${DATA}_formatted_prodigal.faa
		rapsearch -q ${DIR}/${DATA}_formatted_prodigal.faa -d /usr/local/Programs/ncbi-blast-2.2.28+/db/rap2.16_nr -o ${DIR}/${DATA}_rapsearch_out -z ${NPROC} -v 1 -b 1 -t a -p t -a t
		RAPSEARCHOUT=${DIR}/${DATA}_rapsearch_out.aln
		rm prodigal.faa
	fi
fi
grep "^[0-9]" ${RAPSEARCHOUT} | sed 's/\t/|/g' | sed 's/ vs gi//g' | awk -F"|" '{print $1"\t"$2}' > ${DIR}/gi_mapping
sed 's/_/\t/g' ${DIR}/gi_mapping > ${DIR}/3col_gi_mapping
grep -v ">" ${DIR}/${DATA}_rapsearch_out.m8 | awk -F"|" '{print $2}' | sed -e '1,5d' > ${DIR}/gi.txt
GI=`cat ${DIR}/gi.txt`
COUNTER=0
for i in $GI; do
        COUNTER=$((COUNTER+1))
        if [ $COUNTER -lt 250 ]; then
                TEMP=$TEMP$i","
        else
                COUNTER=0
                LIST=`echo ${TEMP} | sed 's/,$//g'`
                echo "esearch -db protein -query \""${LIST}"\" | efetch -format gpc | grep --color -E \"INSDSeq_taxonomy|<INSDSeqid>gi\" | sed -e 's/^[ \\\t]*//' -e 's/<INSDSeq_taxonomy>//g' -e 's/<\\/INSDSeq_taxonomy>//g' -e 's/<INSDSeqid>gi|//g' -e 's/<\/INSDSeqid>//g' >> taxa.txt" >> Esearcher_subscript.sh
                TEMP=`echo " "`
        fi
done
LIST=`echo ${TEMP} | sed 's/,$/ [gi]/g'`
echo "esearch -db protein -query \""${LIST}"\" | efetch -format gpc | grep --color -E \"INSDSeq_taxonomy|<INSDSeqid>gi\" | sed -e 's/^[ \\\t]*//' -e 's/<INSDSeq_taxonomy>//g' -e 's/<\\/INSDSeq_taxonomy>//g' -e 's/<INSDSeqid>gi|//g' -e 's/<\/INSDSeqid>//g' >> taxa.txt" >> Esearcher_subscript.sh
if [ -z ${QUICKPARAM} ]; then
	echo "Finding taxanomy for gene ids."
	rm ${DIR}/taxa.txt 2>/dev/null
	sh Esearcher_subscript.sh >/dev/null 2>/dev/null
	if [ -n "$DATA" ]; then
		mv ${DIR}/taxa.txt ${DIR}/${DATA}_taxa.txt
	fi
fi
rm Esearcher_subscript.sh
grep -A1 "^[0-9]" ${DIR}/${DATA}_taxa.txt | sed '/[0-9]/ a\HOLDER' | sed ':a;N;$!ba;s/\nHOLDER\n/;/g' | sed 's/ //g' | sed -e 's/;/-/7' > ${DIR}/gi_taxa.txt
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
echo "colnames(lineage) <- c(\"gi\", \"Kingdom\", \"Phylum\", \"Class\", \"Order\", \"Family\")" >> created_r_merger
echo "mapping <- read.csv2(file = {MAPPING}, header = F, sep =\"\\\t\")" >> created_r_merger
echo "colnames(mapping) <- c(\"contig\", \"gene\", \"gi\")" >> created_r_merger
echo "merged <- merge(x = lineage, y = mapping, by = \"gi\", all.y = T)" >> created_r_merger
echo "write.csv(merged, quote = F,file = {OUT})" >> created_r_merger
echo ${DIR} | Rscript ${DIR}/created_r_merger >/dev/null 2>/dev/null
rm created_r_merger
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
cat tax_search_param | Rscript ${DIR}/created_r_mode >/dev/null 2>/dev/null
rm created_r_mode tax_search_param
grep "${LABEL}" ${DIR}/table_of_modes | sed -e "s/$LABEL//" | sed 's/"//g'| sed ':a;N;$!ba;s/\n/\$\n\>/g' | sed 's/,//g' | sed '1s/^/>/' | sed "\$a$" | sed ':a;N;$!ba;s/\n\$/$/g' > ${DIR}/line_sep_grep_commands
while read p; do
  grep -A1 "$p" ${DIR}/renumbered.fa >> ${DIR}/target_contigs.fa
done <${DIR}/line_sep_grep_commands
rm line_sep_grep_commands
if [ -z $SUMMARY ]; then
	echo "#!/usr/bin/Rscript" > created_r_summary
	echo "f <- file(\"stdin\")" >> created_r_summary
	echo "open(f)" >> created_r_summary
	echo "while(length(line <- readLines(f,n=1)) > 0) {" >> created_r_summary
	echo "  FASTA<-paste(line[1],\"/target_contigs.fa\", sep=\"\")" >> created_r_summary
	echo "}" >> created_r_summary
	echo "suppressMessages(library(Biostrings))" >> created_r_summary
	echo "target <- readDNAStringSet({FASTA})" >> created_r_summary
	echo "target_frame <- data.frame(name=\"target\",length=width(target))" >> created_r_summary
	echo "print(\"Total length of assembled contigs:\")" >> created_r_summary
	echo "paste(\" \", sum(width(target)), sep =\"\")" >> created_r_summary
	echo "print(\" \")" >> created_r_summary
	echo "print(\"Number of target contigs:\")" >> created_r_summary
	echo "paste(\" \", length(target), sep=\"\")" >> created_r_summary
	echo "print(\" \")" >> created_r_summary
	echo "print(\"Average target contig length:\")" >> created_r_summary
	echo "paste(\" \", sum(width(target))/length(target), sep=\"\")" >> created_r_summary
        echo "print(\" \")" >> created_r_summary
	echo "print(\"N50:\")" >> created_r_summary
	echo "length_range <- unlist(tapply(target_frame\$length, target_frame\$length, function(x) rep(x[1], sum(x))))" >> created_r_summary
	echo "paste(\" \", median(length_range), sep=\"\")" >> created_r_summary
	echo ${DIR} | Rscript ${DIR}/created_r_summary > ${DIR}/CLEANER
	cat ${DIR}/CLEANER | sed 's/\[1\]//g' | sed 's/"//g' | sed 's/^ //g'
	rm ${DIR}/CLEANER created_r_summary
fi
if [ -n "$DATA" ]; then
        mv ${DIR}/target_contigs.fa ${DIR}/${DATA}_${LABEL}_contigs.fa
        mv ${DIR}/table_of_modes ${DIR}/${DATA}_table_of_modes
	echo "\nOutput has been written to: "${DIR}"/"${DATA}"_table_of_modes."
fi
rm no_duplicate_csv merged_r_output gi_taxa.txt gi.txt 3col_gi_mapping gi_mapping 2>/dev/null
