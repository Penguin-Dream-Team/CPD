#!/bin/bash

TESTS=(./tests/*.tree)
BALLQUERIES="ballqueries"

# replace with EXE=echo to not run the ballAlg
EXE=./ballAlg
#EXE=echo

# output directory of ballalg
OUTPUTS="outputs"
BALLQUERYOUT="ballqueryOut"

# set to false if not meant to run python verifier
VERIFY=false

# create output dir
mkdir -p $OUTPUTS
mkdir -p $BALLQUERYOUT

NUMARGS=${#@}

for i in "${!TESTS[@]}";
do
    # filename preprocessing (extract dims nodes and seed)
    FILE=${TESTS[$i]}
    IFS="." read -a f <<< $FILE
    IFS="-" read -a f_split <<< ${f[1]}
    DIMS=${f_split[1]}
    NODES=${f_split[2]}
    SEED=${f_split[3]}

    # run ballAlg algorithm
    echo "Tree #$i [${FILE}] => dims: $DIMS | nodes: $NODES | seed: $SEED"
    OUTPUT="$OUTPUTS/$(basename $FILE)"
    if [ $NUMARGS == 0 ]; then
        echo "Running '$EXE $DIMS $NODES $SEED > $OUTPUT'"
        $EXE $DIMS $NODES $SEED > $OUTPUT
    else
        echo "Running '$EXE $DIMS $NODES $SEED > $OUTPUT' with $1 threads"
        OMP_NUM_THREADS=$1 && $EXE $DIMS $NODES $SEED > $OUTPUT
    fi

    # if success run python verifier
    RESULT=$?
    if [ $RESULT == 0 ]; then
        echo "Finished successfully"

        if [ "$VERIFY" = true ]; then
            python3 verifier.py $FILE $OUTPUT
        fi
    fi

    # print new line if not the last tree
    if [ ${#TESTS[@]} != $(($i + 1)) ]; then
        echo ""
    fi
done

echo ""
echo "Running ballQuery"

# Ball Query stuff
function runBallQuery() {
    echo "Querying for ($2) in $OUTPUTS/$1"
    ./ballQuery $OUTPUTS/$1 $2 > $BALLQUERYOUT/$1
    diff -Z $BALLQUERIES/$1 $BALLQUERYOUT/$1
    RESULT=$?
    if [ $RESULT == 0 ]; then
        echo "No difference found"
    fi
    echo ""
}

runBallQuery ex-2-5-0.tree "3 1"
runBallQuery ex-2-8-0.tree "8 8"
runBallQuery ex-2-6-8.tree "3 7"
runBallQuery ex-3-3-3.tree "5 5 5"
runBallQuery ex-5-4-3.tree "9 7 5 3 1"
runBallQuery ex-3-5000000-0.tree "4 5 6"
runBallQuery ex-20-1000000-0.tree "1 2 3 4 5 6 7 8 9 1 2 3 4 5 6 7 8 9 1 2"
runBallQuery ex-3-20000000-0.tree "1 5 9"
