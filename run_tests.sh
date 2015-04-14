#!/bin/bash

failed=0
total=0
for tst in ./build/*.tst;
do
	echo "Running test $tst..."
	$tst
	if [ $? -eq 0  ]
	then
		echo "  [OK]"
	else
		let failed+=1
		echo "  [ER]"
	fi
	let total+=1
	echo
done

echo "$failed/$total tests failed"
