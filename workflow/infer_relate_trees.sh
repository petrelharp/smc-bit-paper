#binaries from website
RELATE_DIR=`pwd`/../../../relate
#git clone https://github.com/leospeidel/relate_lib.git
#cd relate_lib && mkdir -p build && cd build
#cmake .. && make && cd ../..
RELATE_LIB=`pwd`/../../../relate_lib 

#hapmap
printf '
cl <- %s
cat("Position(bp) Rate(cM/Mb) Map(cM)\n")
cat("0 1.0 0.0\n")
cat(as.integer(cl), "0.0", cl/1e6, "\n")
' $SEQLEN | R --slave >$ID.hapmap

NE="20000" #haploid
MU="1.25e-8"
SEED="1024"
ID="chr1"
SEQLEN=1000000
SAMPLES=100

WORK_DIR=`pwd`/sim_$SEED
IN_DIR=$WORK_DIR/relate_inputs
OUT_DIR=$WORK_DIR/relate_outputs
mkdir -p $WORK_DIR
mkdir -p $IN_DIR
mkdir -p $OUT_DIR

#dump into haplotypes
$RELATE_DIR/bin/RelateFileFormats \
  --mode ConvertFromVcf --haps $ID.haps \
  --sample samples -i $ID

#ancestral states
awk 'NR > 1 { printf(">%s\n%s\n", $3, $4); }' ID.vcf >$ID.anc.fa

#infer trees
$RELATE_DIR/scripts/PrepareInputFiles/PrepareInputFiles.sh \
  --haps $ID.haps \
  --sample samples \
  --ancestor $ID.anc.fa \
  -o $ID
$RELATE_DIR/bin/Relate --mode All -m $MU -N $NE \
  --haps $ID.haps.gz --sample sample.gz --annot $ID.annot \
  --dist $ID.dist.gz --map $ID.hapmap -o $ID

# convert to tree sequence with edges linked across trees
$RELATE_LIB/bin/Convert --mode ConvertToTreeSequence \
  --anc $OUT_DIR/$ID.anc.gz --mut $OUT_DIR/$ID.mut.gz --iterations 1000 --compress -o $OUT_DIR/$ID