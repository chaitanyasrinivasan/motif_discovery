#!/bin/bash

if ! [ -x "$(command -v git)" ];
then
	echo "yes"
else 
	echo "no"
fi