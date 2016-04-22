#!/bin/bash
set -e						

if [ $# -eq 0 ]; then
    echo "Usage: $0 DataFile"
	echo "It is better to move your data to separet folder before running scripts. It produce a lot of mid files."
	exit 1
fi

SELF_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"	# DO NOT EDIT

################################### EDITABLE PART ############################################
PET_COUNT_CUTOFF=4			              # It is the smallest value for which the interaction will be treated as a cluster, not a singleton
KARYOTYPE="$SELF_PATH/karyotype.hg19.txt" # File with karyotype in Circos format
CELL_LINE="MyExperiment"				  # Prefix in some files. Can not be empty!
RESOLUTION=10000				          # Signal file resolution. Lower the better 100-10000 should be ok. Lower values are time consuming. Larger may couse rounding errors.
##############################################################################################

# Enviroment vars. Do not edit!
CURRENT_PATH=`pwd`
cd $CURRENT_PATH
GETSIGNAL=$SELF_PATH/getSignal.py
SEGMENTATE=$SELF_PATH/segmentate.py
INPUT_PATH=$1
INPUT_FILENAME=$(basename "$INPUT_PATH")
INPUT_DIR=$(dirname "$INPUT_PATH")
FILENAME_1=PET$PET_COUNT_CUTOFF.$INPUT_FILENAME

cd $INPUT_DIR

# Trim according to PET-Count
echo -n "Trimming according to PET-Count $PET_COUNT_CUTOFF... "
awk ' $7 >= 4 ' $INPUT_FILENAME > $FILENAME_1
echo "OK"

ORIGINAL_INTERACTIONS_NUMBER=$(wc -l $INPUT_FILENAME | cut -f1 -d ' ')
TRIMMED_INTERACTIONS_NUMBER=$(wc -l $FILENAME_1 | cut -f1 -d ' ')
ORIGINAL_FILESIZE=$(du -h "$INPUT_FILENAME" | cut -f1)
TRIMMED_FILESIZE=$(du -h "$FILENAME_1" | cut -f1)
SIZE_REDUCTION=$(python -c "print '{:0.2f} %'.format(-($TRIMMED_INTERACTIONS_NUMBER - $ORIGINAL_INTERACTIONS_NUMBER)/float($ORIGINAL_INTERACTIONS_NUMBER)*100) ")

echo "Original file contain $ORIGINAL_INTERACTIONS_NUMBER interactions ($ORIGINAL_FILESIZE)"
echo "Trimmed file contain  $TRIMMED_INTERACTIONS_NUMBER interactions ($TRIMMED_FILESIZE)"
echo "Interaction number reduction: $SIZE_REDUCTION"

# Extract the list of chromosomes
echo  "Exracting chromosomes list... "
LIST_OF_CHROMOSOMES=$(cut -f1,4 "$FILENAME_1" | awk '{print $1"\n"$2}' | sort -Vu)
LIST_OF_CHROMOSOMES_IN_KARYOTYPE=$(cat $KARYOTYPE | grep ^chr | cut -f7 -d' ' | sort -V)
echo $LIST_OF_CHROMOSOMES
echo $LIST_OF_CHROMOSOMES_IN_KARYOTYPE
for i in $LIST_OF_CHROMOSOMES; do
	if [[ $LIST_OF_CHROMOSOMES_IN_KARYOTYPE =~ $i ]]; then
		:
	else
		 echo "[ERROR] $i was not found in karyotype file ($KARYOTYPE). This inconsistency must be corected before further processing. Remeber that chromosomes are case sensitive," 
		exit 1
	fi
done
echo "OK"

# Splitting file by chromosomes
echo "Splitting file by chromosomes... "
for i in $LIST_OF_CHROMOSOMES; do
	FILENAME_2=$i.$FILENAME_1
	set +e
	grep -P "^$i\t" $FILENAME_1 | grep -P "\t$i\t" > $FILENAME_2
	set -e
	echo -e "$i\t$(wc -l $FILENAME_2 | cut -f1 -d ' ')"
done
echo "OK"

echo "Converting signal into np.uint16..."
for i in $LIST_OF_CHROMOSOMES; do
	FILENAME_2=$i.$FILENAME_1
	$GETSIGNAL -i -p $CELL_LINE -r $RESOLUTION $KARYOTYPE $FILENAME_2 $i
done

echo "Performing segmentation..."
for i in $LIST_OF_CHROMOSOMES; do
	FILENAME_3=$CELL_LINE.$i.signal.np.uint16 # signal file name
	set +e
	$SEGMENTATE $KARYOTYPE $FILENAME_3 $i 
	set -e
done

OUT_FILE=$CELL_LINE.segments.txt

#cat *.seg.txt | sort -V > $OUT_FILE
find . -name '*.segments.txt' ! -name $OUT_FILE | xargs cat | sort -V > $OUT_FILE

echo "OK"
echo "Results saved in $OUT_FILE"
echo "Everything is finished."

cd $CURRENT_PATH
