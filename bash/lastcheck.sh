#!/bin/bash

checkpoint_iter=""
checkpoint_name=""
line_count=0

# Read from your file
while IFS= read -r line; do
    if (( line_count == 0 )); then
        # Assign the line to checkpoint_iter
        checkpoint_iter=$line
        checkpoint_name=$checkpoint_iter
    elif (( line_count == 1 )); then
        checkpoint_name=quicksave_$line
    fi
    ((line_count++)) # Increment line count for each line read
done < Checkpoints/last_checkpoint

# make tar ball
echo "Creating lastcheck_$checkpoint_iter.tar"
echo "........................."
tar cvf lastcheck_$checkpoint_iter.tar Checkpoints/last_checkpoint Checkpoints/checkpoint_log Checkpoints/$checkpoint_name
