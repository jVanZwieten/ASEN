Conventions
===

All of these thoughts essentially spring from me realizing that vectors are sometimes denoted with capital letters, yet programming typically has extensive convention concerning capitalization. What follows is my exploration into defining a convention to make me efficient and consistent writing scientific software.

# standard forms
Keep data to these forms whenever possible, only convert when necessary then convert back when possible.
* datetime
* radians
* row-epoch
* timetables

# variable notation
Coding often encodes much information through variable naming conventions. Written math and science encodes a lot of information that is difficult to express succinctly yet necessary to persist accurately. It is therefore challenging to write in the space of both concerns simultaneously, expressing information pertinent both to programming and mathematical formulae in a accurate, expressive, yet concise way. Conventions should aid in debugging both mathematical operation and application logic.

## verbosity
Programmers often opt to abbreviate a lot in their code, and I've seen that to an egregious extent in MATLAB code. There was a time in programming history when character count could have significant effect on performance and a real burden on development. With modern compilers and development tools that's no longer relevant, and ease of debugging is a far more salient concern. Generally, opt to be more clear and obvious rather than brief with variable names.
However, it's often easier to read complex formulae with very concise variable names. Therefore achieving consistent, dense information encoding can be valuable.
* save long names into short but obvious variables before significant manipulation
* make it easy/local to reference what a short variable is referring to

## domain knowledge
The purpose of maintainable code is not to explain the physics or equation derivation to a novice. It is reasonable to expect that the reader has some domain knowledge to understand what he's reading. The goal is to make it easy to spot an error in an operation if the maintainer were comparing the code to written formula.

## notation
[value][vector modifier][reference frame]_[owner][instance]

### vector modifiers
* V = vector (nx1)
* Hat = unit vector
* M = matrix (nxm)
* d = delta
* Dot = rate (derivative)
* DotDot ... multiple derivative
* T = tilde

If variables a and aV exist, a is the magnitude of aV.