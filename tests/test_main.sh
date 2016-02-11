#! /bin/bash 

tests_passed=0
tests_failed=0

program="zkc-test"

k_sizes=(13 15)
prefixes=("standard" "with_ns" "blank_line_end" "no_quals" "end_at_plus")


for file_prefix in "${prefixes[@]}"; do

	if [ $file_prefix == "standard" ]; then
		echo "Testing zkc hist for basic files"
		extensions=("fasta" "fastq")

	elif [ $file_prefix == "with_ns" ]; then
		echo "Testing zkc hist for files containing Ns"
		extensions=("fasta" "fastq")

	elif [ $file_prefix == "blank_line_end" ]; then
		echo "Testing zkc hist for files ending in a blank line"
		extensions=("fasta" "fastq")

	elif [ $file_prefix == "no_quals" ]; then
		echo "Testing zkc hist for fastq files ending without quality values"
		extensions=("fastq")

	elif [ $file_prefix == "end_at_plus" ]; then
		echo "Testing zkc hist for fastq files ending after the plus but before the quality values"
		extensions=("fastq")

	fi

	cd $file_prefix
		
	for K in "${k_sizes[@]}"; do

		for extension in "${extensions[@]}"; do

			desired_input=$file_prefix"."$extension

			$program hist -k $K -c $desired_input > stdout.tmp 2> stderr.tmp

			desired_stdout=$file_prefix"."$K"mer_hist.canonical"

			if cmp stdout.tmp $desired_stdout
			then 
				((tests_passed++))
			else
				((tests_failed++))
				echo "Stdout test fails for "$desired_input", k = "$K" canonical"
			fi

			desired_stderr=$file_prefix".stderr"

			if cmp stderr.tmp $desired_stderr
			then 
				((tests_passed++))
			else
				((tests_failed++))
				echo "Stderr test fails for "$desired_input", k = "$K" canonical"
			fi

			rm stdout.tmp stderr.tmp


			$program hist -k $K $desired_input > stdout.tmp 2> stderr.tmp

			desired_stdout=$file_prefix"."$K"mer_hist.not_canonical"

			if cmp stdout.tmp $desired_stdout
			then 
				((tests_passed++))
			else
				((tests_failed++))
				echo "Stdout test fails for "$desired_input", k = "$K" non-canonical"
			fi

			if cmp stderr.tmp $desired_stderr
			then 
				((tests_passed++))
			else
				((tests_failed++))
				echo "Stderr test fails for "$desired_input", k = "$K" non-canonical"
			fi

			rm stdout.tmp stderr.tmp
		done
	done

	if [ $file_prefix == "standard" ]; then
		echo "Finished testing zkc hist for basic files"

	elif [ $file_prefix == "with_ns" ]; then
		echo "Finished testing zkc hist for files containing Ns"

	elif [ $file_prefix == "blank_line_end" ]; then
		echo "Finished testing zkc hist for files ending in a blank line"

	elif [ $file_prefix == "no_quals" ]; then
		echo "Finished testing zkc hist for fastq files ending without quality values"

	elif [ $file_prefix == "end_at_plus" ]; then
		echo "Finished testing zkc hist for fastq files ending after the plus but before the quality values"

	fi

	cd ..

done


printf "\nNumber of tests passed: "$tests_passed"\n"
echo "Number of tests failed: "$tests_failed

if [ $tests_failed == 0 ] 
then
	printf "\nWell done, all tests passed!\n"
fi

printf "\n"
