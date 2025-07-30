#!/bin/bash

checkpoint_iter=""
checkpoint_name=""
line_count=0

# Read from your file
while IFS= read -r line; do
    if (( line_count == 0 )); then
        # Assign the line to checkpoint_iter
        checkpoint_iter=$line
        # Take the absolute value of checkpoint_iter
        if (( checkpoint_iter < 0 )); then
            checkpoint_iter=$(( -checkpoint_iter )) # Multiply by -1 to get the absolute value
        fi
        checkpoint_name=$checkpoint_iter
    elif (( line_count == 1 )); then
        checkpoint_name=quicksave_$line
    fi
    ((line_count++)) # Increment line count for each line read
done < Checkpoints/last_checkpoint

echo "Variable 1 (integer): $checkpoint_iter"
echo "Variable 2 (integer): $checkpoint_name"
