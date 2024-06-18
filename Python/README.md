# Introduction to Python

[Python](https://www.python.org) is one of the most popular programming languages. There are lots of reasons why, 
here are some of them.

- [It is very easy to use](https://stackoverflow.blog/2021/07/14/getting-started-with-python/)
- There are lots and lots of tools for different things that have already been written and tested by others
    - The scientific python suite includes many interconnected modules
        - [numpy](http://www.numpy.org/) - the numeric Python library that excels at matrix manipulations and common 
          mathematical functions
        - [scipy](https://www.scipy.org/) - the scientific Python library that includes statistics, 
          modeling, and much more
        - [matplotlib](http://matplotlib.org/) - the graphical library for making plots and graphs
        - [seaborn](https://seaborn.pydata.org/) - an extension to matplot lib that makes beautiful graphs easy 
        - [pandas](http://pandas.pydata.org/) - rich data structures (especially tables) for manipulating data
        - [networkx](https://networkx.github.io/) - for network analysis, interrogation, and manipulation
        - [sci-kit learn](http://scikit-learn.org/) - for machine learning
    - There are many, many [domain specific add ons](https://www.scipy.org/topical-software.html) that you
    can leverage when you are writing your code and analyzing your data. For example, some of the topics include:
        -   [Astronomy](https://www.scipy.org/topical-software.html#astronomy)
        -   [Artificial intelligence and machine learning](https://www.scipy.org/topical-software.html#artificial-intelligence-and-machine-learning)
        -   [Bayesian statistics](https://www.scipy.org/topical-software.html#bayesian-statistics)
        -   [Biology](https://www.scipy.org/topical-software.html#biology-including-neuroscience)
        -   [Dynamical systems](https://www.scipy.org/topical-software.html#dynamical-systems)
        -   [Economics and econometrics](https://www.scipy.org/topical-software.html#economics-and-econometrics)
        -   [Electromagnetics and electrical engineering](https://www.scipy.org/topical-software.html#electromagnetics-and-electrical-engineering)
        -   [Geosciences](https://www.scipy.org/topical-software.html#geosciences)
        -   [Molecular modeling](https://www.scipy.org/topical-software.html#molecular-modeling)
        -   [Signal processing](https://www.scipy.org/topical-software.html#signal-processing)
        -   [Symbolic math, number theory, etc.](https://www.scipy.org/topical-software.html#symbolic-math-number-theory-etc)
        -   [Quantum mechanics](https://www.scipy.org/topical-software.html#quantum-mechanics)
- It is one of the most popular langauges on [Stack Overflow](https://stackoverflow.com/questions/tagged/python) with 
  over 1.75 million questions and answers!
- The [PyCharm](https://www.jetbrains.com/pycharm/) development enironment is free!

# The pros and cons of Python

Python is no where near the fastest computer programming language. In fact, it is even slower than some similar
languages [like Perl](https://stackoverflow.com/questions/12793562/text-processing-python-vs-perl-performance), 
but it is much easier to develop in Python than many other languages (e.g. Perl). In a _lot_ of applications 
(e.g. bioinformatics, data science, geoscience), the development time is much more important than the computer
run time. Note that you should still think about [what affects computer run time and how to measure 
that](https://youtu.be/IgeJmTKQlKs), and sometimes you will run into intractable problems, but often the cheapest
solution is to throw more computers at a problem, not to spend more time developing it. _Note_ this is not always 
true, and there are many occassions (e.g. running things on cell phones) where you want to be very careful to 
optimize your applications!


# Learning Python for Bioinformatics

We have adapted Marc Cohen's Google Colab notebooks that teach Python to be more aligned with bioinformatics. This series of Python notebooks will
walk you through the Python basics, and introduce you to more advanced concepts as you progress.

You can access the [first Google Colab notebook here](https://colab.research.google.com/drive/1tw4j91fG8z_et5atakifenjC7Uge2Uoh) 

Or you can jump to a specific lesson. Each of these links opens in Google Colab. You should make a copy of the file and run the code for yourself.

* [Lesson 0 - Introduction and Index](https://colab.research.google.com/drive/1tw4j91fG8z_et5atakifenjC7Uge2Uoh)
* [Lesson 1 - Variables and Types](https://colab.research.google.com/drive/11rcjrwznjZS-qoFoEcCQYF85PCRtN9dr)
  * Variables
  * Naming variables
  * Types of data
  * Numeric Types
  * String Types
  * Using Variables in Python
  * Built in Python functions
* [Lesson 2 - Expressions](https://colab.research.google.com/drive/1Sm7N8Agf0aFj6qbd6GwenGBaqchZYTug)
  * Constants vs. Variables
  * Data Types
  * The Boolean (bool) Type
  * The None Type
  * Comparison Opererators
  * Boolean Operators - and, or, and not
  * Order of Evaluation
  * Python Precedence Rules
  * F strings
* [Lesson 3 Lists](https://colab.research.google.com/drive/1NbWawPfWAQV2x56rG0SvcMNpXI7sn3R0)
  * Lists
  * Creating Lists
  * Sets
  * List Operations
* [Lesson 4 - Dictionaries](https://colab.research.google.com/drive/1IyjNTpdtwaulP_QXbrCYK8Yd2J3MgCQo)
  * Dictionary Operations
  * Rule of thumb for truth value of lists, and dictionaries
* [Lesson 5 - Conditionals](https://colab.research.google.com/drive/1VmGd4AAb1fBKOjmemYIKnPgu58xGE5so)
  * Controlling Program Flow
  * if Statements
  * if Block Structure
  * Python's use of whitespace
  * else Statements
  * elif Statements
  * For loops
* [Lesson 6 - Functions](https://colab.research.google.com/drive/1iLE4rkhVxige0za_H0ehtb9NNzzxQIwu)
  * Defining Functions
  * Docstrings
  * Return Values
* [Lesson 7 - Reading and Writing Files](https://colab.research.google.com/drive/1Uq9ysM5TxMsiS9vElA53Ihl63ONOYDTA)
  * Reading a fasta file
  * Writing to files
  * Reading and writing gzip files
* [Lesson 8 - Modules](https://colab.research.google.com/drive/1ytI6exHPmHDvtzwRsS-zELWlJBOig9L9)
  * The from Statement
  * When to use import vs. from
* [Lesson 9 - Finding and installing other peoples code](https://colab.research.google.com/drive/1JGRJpUPKkkVukyNvtfEJYVVCcdpkyRLZ)
* [Lesson 10 - Translating a DNA sequence](https://colab.research.google.com/drive/1trXzcwT0VnmdnVQY_Wj9b__pXVY8_7GJ)
  * Challenges to try!
* [Lesson 11 - BioPython](https://colab.research.google.com/drive/1N2WL7WDjUQkb7BLWqKYALWCwsUEAVVdf)
  * Translating a DNA sequence
  * Reading a fasta file
  * Reading a fastq file
  * Reading GenBank Files
* [Lesson 12 - Plotting data with Pandas and Seaborn](https://colab.research.google.com/drive/1Qbff17-gZbktQliV6TaFFupDpfUiAnBE)
  * Pandas
  * Seaborn



