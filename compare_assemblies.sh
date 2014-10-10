#!/bin/bash
while test 1 -gt 0; do
	case "$1" in
		-h)
			echo "options:"
			echo "-h,	show brief help"
			echo "-f1,	first fasta to be compared"
			echo "-l1,	first fasta label"
                        echo "-f2,      second fasta to be compared(optional)"
                        echo "-l2,      second fasta label(only required with -f2 flag)"
			echo "-f3, 	third fasta to be compared(optional)"
			echo "-l3,	third fasta label(only required with -l3 flag)"
			echo "-o,	output name for graph"
                        exit 0
                        ;;
                -f1)
                        shift
                        if test $# -gt 0; then
                                export FASTA1=$1
                        else
                                echo "First fasta not provided"
                                exit 1
                        fi
                        shift
                        ;;
                -l1)
                        shift
                        if test $# -gt 0; then
                                export LABEL1=$1
                        else
                                echo "First label not provided"
                                exit 1
                        fi
                        shift
                        ;;
                -f2)
                        shift
                        if test $# -gt 0; then
                                export FASTA2=$1
                        else
                                echo "Second fasta not provided"
                                exit 1
                        fi
                        shift
                        ;;
                -l2)
                        shift
                        if test $# -gt 0; then
                                export LABEL2=$1
                        else
                                echo "Second label not provided"
                                exit 1
                        fi
                        shift
                        ;;
                -f3)
                        shift
                        if test $# -gt 0; then
                                export FASTA3=$1
                        else
                                echo "Third fasta not provided"
                                exit 1
                        fi
                        shift
                        ;;
                -l3)
                        shift
                        if test $# -gt 0; then
                                export LABEL3=$1
                        else
                                echo "Third label not provided"
                                exit 1
                        fi
                        shift
                        ;;
                -o)
                        shift
                        if test $# -gt 0; then
                                export OUTPUT=$1
                        else
                                echo "Output label not provided"
                                exit 1
                        fi
                        shift
                        ;;
                *)
                        break
                        ;;
        esac
done
if [ -z $FASTA1 ]; then
        echo "First fasta not specified, exiting."
        exit
fi
if [ -z $LABEL1 ]; then
        echo "Fasta Label 1 not specified, exiting."
        exit
fi
if [ -z $OUTPUT ]; then
        echo "Graph output not specified, writing to output_graph."
        OUTPUT="output_graph"
fi
if [ -z $FASTA2 ]; then
  echo using 1 fasta
  echo $FASTA1"\n"$LABEL1"\n"$OUTPUT > input_temp
  echo "#!/usr/bin/Rscript" > created_r_script
  echo "f <- file('stdin')" >> created_r_script
  echo "open(f)" >> created_r_script
  echo "while(length(line <- readLines(f,n=3)) > 0) {" >> created_r_script
  echo "  FASTA1<-print(line[1])" >> created_r_script
  echo "  FASTA1label<-print(as.character(line[2]))" >> created_r_script
  echo "  OUTPUT_NAME<-print(line[3])" >> created_r_script
  echo "}" >> created_r_script
  echo "library(Biostrings, quietly = TRUE)" >> created_r_script
  echo "library(ggplot2, quietly = TRUE)" >> created_r_script
  echo "print(FASTA1label)" >> created_r_script
  echo "fasta1 <- readDNAStringSet({FASTA1})" >> created_r_script
  echo "fasta1_frame <- data.frame(name=\"fasta1\",length=width(fasta1))" >> created_r_script
  echo "all_frame <- rbind(fasta1_frame)" >> created_r_script
  echo "output_graph <- ggplot(all_frame, aes(x=length, fill=name)) + geom_histogram(data=fasta1_frame, alpha=0.2, binwidth=0.01) + scale_x_log10()" >> created_r_script
  echo "jpeg(filename = {OUTPUT_NAME},width = 1024, height = 768, units = \"px\")" >> created_r_script
  echo "output_graph + scale_fill_discrete(name=\"Assemblies\", breaks=c(\"fasta1\"),labels=c(FASTA1label))" >> created_r_script
  echo "dev.off()" >> created_r_script
  cat input_temp | Rscript created_r_script >/dev/null 2>/dev/null
  echo "Graph has been written to:"
  echo $OUTPUT
  rm input_temp
  rm created_r_script
  exit 1
fi
if [ -z $FASTA3 ]; then
  echo using 2 fastas
  echo $FASTA1"\n"$LABEL1"\n"$OUTPUT"\n"$FASTA2"\n"$LABEL2 > input_temp
  echo "#!/usr/bin/Rscript" > created_r_script
  echo "f <- file('stdin')" >> created_r_script
  echo "open(f)" >> created_r_script
  echo "while(length(line <- readLines(f,n=5)) > 0) {" >> created_r_script
  echo "  FASTA1<-print(line[1])" >> created_r_script
  echo "  FASTA1label<-print(as.character(line[2]))" >> created_r_script
  echo "  OUTPUT_NAME<-print(line[3])" >> created_r_script
  echo "  FASTA2<-print(line[4])" >> created_r_script
  echo "  FASTA2label<-print(line[5])" >> created_r_script
  echo "}" >> created_r_script
  echo "library(Biostrings, quietly = TRUE)" >> created_r_script
  echo "library(ggplot2, quietly = TRUE)" >> created_r_script
  echo "print(FASTA1label)" >> created_r_script
  echo "print(FASTA2label)" >> created_r_script
  echo "fasta1 <- readDNAStringSet({FASTA1})" >> created_r_script
  echo "fasta2 <- readDNAStringSet({FASTA2})" >> created_r_script
  echo "fasta1_frame <- data.frame(name=\"fasta1\",length=width(fasta1))" >> created_r_script
  echo "fasta2_frame <- data.frame(name=\"fasta2\",length=width(fasta2))" >> created_r_script
  echo "all_frame <- rbind(fasta1_frame, fasta2_frame)" >> created_r_script
  echo "output_graph <- ggplot(all_frame, aes(x=length, fill=name)) + geom_histogram(data=fasta1_frame, alpha=0.2, binwidth=0.01) + geom_histogram(data=fasta2_frame, alpha=0.2, binwidth=0.01) + scale_x_log10()" >> created_r_script
  echo "jpeg(filename = {OUTPUT_NAME},width = 1024, height = 768, units = \"px\")" >> created_r_script
  echo "output_graph + scale_fill_discrete(name=\"Assemblies\", breaks=c(\"fasta1\", \"fasta2\"),labels=c(FASTA1label,FASTA2label))" >> created_r_script
  echo "dev.off()" >> created_r_script
  cat input_temp | Rscript created_r_script >/dev/null 2>/dev/null
  echo "Graph has been written to:"
  echo $OUTPUT
  rm input_temp
  rm created_r_script
  exit 1
fi
if [ -z $FASTA4 ]; then
  echo Generating plot with 3 fastas
  echo $FASTA1"\n"$LABEL1"\n"$OUTPUT"\n"$FASTA2"\n"$LABEL2"\n"$FASTA3"\n"$LABEL3 > input_temp
  echo "#!/usr/bin/Rscript" > created_r_script
  echo "f <- file('stdin')" >> created_r_script
  echo "open(f)" >> created_r_script
  echo "while(length(line <- readLines(f,n=7)) > 0) {" >> created_r_script
  echo "  FASTA1<-print(line[1])" >> created_r_script
  echo "  FASTA1label<-print(as.character(line[2]))" >> created_r_script
  echo "  OUTPUT_NAME<-print(line[3])" >> created_r_script
  echo "  FASTA2<-print(line[4])" >> created_r_script
  echo "  FASTA2label<-print(line[5])" >> created_r_script
  echo "  FASTA3<-print(line[6])" >> created_r_script
  echo "  FASTA3label<-print(line[7])" >> created_r_script
  echo "}" >> created_r_script
  echo "library(Biostrings, quietly = TRUE)" >> created_r_script
  echo "library(ggplot2, quietly = TRUE)" >> created_r_script
  echo "print(FASTA1label)" >> created_r_script
  echo "print(FASTA2label)" >> created_r_script
  echo "print(FASTA3label)" >> created_r_script
  echo "fasta1 <- readDNAStringSet({FASTA1})" >> created_r_script
  echo "fasta2 <- readDNAStringSet({FASTA2})" >> created_r_script
  echo "fasta3 <- readDNAStringSet({FASTA3})" >> created_r_script
  echo "fasta1_frame <- data.frame(name=\"fasta1\",length=width(fasta1))" >> created_r_script
  echo "fasta2_frame <- data.frame(name=\"fasta2\",length=width(fasta2))" >> created_r_script
  echo "fasta3_frame <- data.frame(name=\"fasta3\",length=width(fasta3))" >> created_r_script
  echo "all_frame <- rbind(fasta1_frame, fasta2_frame, fasta3_frame)" >> created_r_script
  echo "output_graph <- ggplot(all_frame, aes(x=length, fill=name)) + geom_histogram(data=fasta1_frame, alpha=0.2, binwidth=0.01) + geom_histogram(data=fasta2_frame, alpha=0.2, binwidth=0.01) + geom_histogram(data=fasta3_frame, alpha=0.2, binwidth=0.01) + scale_x_log10()" >> created_r_script
  echo "jpeg(filename = {OUTPUT_NAME},width = 1024, height = 768, units = \"px\")" >> created_r_script
  echo "output_graph + scale_fill_discrete(name=\"Assemblies\", breaks=c(\"fasta1\", \"fasta2\", \"fasta3\"),labels=c(FASTA1label,FASTA2label,FASTA3label))" >> created_r_script
  echo "dev.off()" >> created_r_script
  cat input_temp | Rscript created_r_script >/dev/null 2>/dev/null
  echo "Graph has been written to:"
  echo $OUTPUT
  rm input_temp
  rm created_r_script
  exit 1
fi
#output_graph <- ggplot(all_frame, aes(x=length, fill=name)) + geom_histogram(data=singletons_frame, alpha=0.2, binwidth=0.01, colour = "darkgreen", fill = "blue") + scale_x_log10(limits=c(1000,10000)) + xlab("Contig Length") + ylab("Number Of Contigs")
