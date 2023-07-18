for DIR in $(find -maxdepth 1 -type d | sed -e 's/.\/\..*//; s/^\.$//; s/\.\///'); do
	echo $DIR
	cd $DIR
	for FILE in *.md; do
		if [ $FILE == "README.md" ]; then
			OUTPUT=$DIR;
		else
			OUTPUT=$(echo $FILE | sed -e 's/.md//');
		fi

		if [ ! -e $OUTPUT.pdf ] || [ "$FILE" -nt "$OUTPUT.pdf" ]; then
			echo -e  "\tCreating PDF of $FILE as $OUTPUT.pdf in $DIR";
			pandoc $FILE -f markdown --latex-engine=xelatex --columns 100 --smart -s -o $OUTPUT.pdf
		fi
	done
	cd ../
done
