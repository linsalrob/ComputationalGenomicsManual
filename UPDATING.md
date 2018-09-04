# How to update this manual

If you wish to contribute to the manual, please contact Rob. He'll welcome all kinds of help! The easiest way to do so is to clone the git repository, and then make your changes in your clone. Once you are done with that, make a pull request here and Rob can merge your changes with his.

If you are not sure how to do any of that, ask Rob and he'll point you in the right direction to get started.


# Creating PDFs

In each directory we include PDFs of the current state of the README. These are made with pandoc, and you can recreate the PDF at any time using the command (of course, chaning FILE to the appropriate name):

```
pandoc README.md -f markdown --latex-engine=xelatex --columns 100 --smart -s -o FILE.pdf
```

or for all directories:

```
for DIR in $(find -maxdepth 1 -type d | sed -e 's/.\/\..*//; s/^\.$//; s/\.\///'); do
	echo $DIR
	cd $DIR
	pandoc README.md -f markdown --latex-engine=xelatex --columns 100 --smart -s -o $DIR.pdf
	cd ../
done
```

