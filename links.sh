mkdir /home/asmyasnikova83/SCR/links

prefix="/net/server/data/Archive/prob_learn/experiment/"
DIRS=$(ls -d ${prefix}P0*/*/)

for DIR in $DIRS; do
    SUBJ=${DIR#"$prefix"}
    SUBJ=${SUBJ:0:4}
    echo "ln -s $DIR links/$SUBJ"
    ln -s $DIR links/$SUBJ
done
