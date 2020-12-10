#!/bin/bash
while : #true
do
    tail -n1100 $(ls -t | head -n1) | grep Time
    sleep 1
done
