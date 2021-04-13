#!/bin/bash

TESTS=(./tests/*.tree)

# replace with EXE=echo to not run the ballAlg
EXE=./ballAlg
#EXE=echo

# output directory of ballalg
OUTPUTS="outputs/"

# set to false if not meant to run python verifier
VERIFY=true

# create output dir
mkdir -p $OUTPUTS

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
    OUTPUT="$OUTPUTS$(basename $FILE)"
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
            python verifier.py $FILE $OUTPUT
        fi
    fi

    # print new line if not the last tree
    if [ ${#TESTS[@]} != $(($i + 1)) ]; then
        echo ""
    fi
done
