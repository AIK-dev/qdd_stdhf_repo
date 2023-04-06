#!/usr/bin/env bash
#
# CALL: ./regularise-ucF90.sh <source-file.F90>
#
SED=gsed
IN=$1

# Array containing C-preprocessor directives
declare -a Directives=(\
"#include" \
"#define" \
"#undef" \
"#if" \
"#elif" \
"#else" \
"#endif" \
"#ifdef" \
"#ifndef" \
"#error" \
"#line" \
"#pragma" \
)

# Array containing all Fortran 90 standard statement keywords (+ 'import') found in
# Appendix A in https://link.springer.com/content/pdf/bbm%3A978-1-4612-2562-1%2F1.pdf
declare -a Keywords=(\
"ACCESS" \
"ACTION" \
"ADVANCE" \
"ALLOCATABLE" \
"ALLOCATE" \
"ASSIGN" \
"ASSIGNMENT" \
"BACKSPACE" \
"BLANK" \
"BLOCK" \
"CALL" \
"CASE" \
"CHARACTER" \
"CLOSE" \
"COMMON" \
"COMPLEX" \
"CONTAINS" \
"CONTINUE" \
"CYCLE" \
"DATA" \
"DEALLOCATE" \
"DEFAULT" \
"DELIM" \
"DIMENSION" \
"DIRECT" \
"DO" \
"DOUBLE" \
"ELSE" \
"ELSEWHERE" \
"END" \
"ENDFILE" \
"ENTRY" \
"EOR" \
"EQUIVALENCE" \
"ERR" \
"EXIST" \
"EXIT" \
"EXTERNAL" \
"FILE" \
"FMT" \
"FORM" \
"FORMAT" \
"FORMATTED" \
"FUNCTION" \
"IF" \
"IMPLICIT" \
"IMPORT" \
"IN" \
"INOUT" \
"INQUIRE" \
"INTEGER" \
"INTENT" \
"INTERFACE" \
"INTRINSIC" \
"IOLENGTH" \
"IOSTAT" \
"KIND" \
"LEN" \
"LOGICAL" \
"MODULE" \
"NAME" \
"NAMED" \
"NAMELIST" \
"NEXTREC" \
"NML" \
"NONE" \
"NULLIFY" \
"NUMBER" \
"ONLY" \
"OPEN" \
"OPENED" \
"OPERATOR" \
"OPTIONAL" \
"OUT" \
"PAD" \
"PARAMETER" \
"PAUSE" \
"POINTER" \
"POSITION" \
"PRECISION" \
"PRINT" \
"PRIVATE" \
"PROCEDURE" \
"PROGRAM" \
"PUBLIC" \
"READ" \
"READWRITE" \
"REAL" \
"REC" \
"RECL" \
"RECURSIVE" \
"RESULT" \
"RETURN" \
"REWIND" \
"SAVE" \
"SELECT" \
"SEQUENCE" \
"SEQUENTIAL" \
"SIZE" \
"STAT" \
"STATUS" \
"STOP" \
"SUBROUTINE" \
"TARGET" \
"THEN" \
"TO" \
"TYPE" \
"UNFORMATTED" \
"UNIT" \
"USE" \
"WHERE" \
"WHILE" \
"WRITE" \
)

TRIMMED=${IN%.F90}.TRIMMED.F90
CASED=${TRIMMED%.F90}.CASED.F90
REGULARISED=${IN%.F90}.REGULARISED.F90

# Condense excessive spaces around strings to 1 space
tr -s " " < $IN > $TRIMMED

cp $TRIMMED $CASED
# Convert all Fortran standard statement keywords to uppercase
for keyword in ${Keywords[@]}
do
  $SED -i "s/\b${keyword}\b/${keyword}/gi" $CASED
done
# Convert all C-preprocessor directives to lowercase
for directive in ${Directives[@]}
do
  $SED -i "s/${directive}/${directive}/gi" $CASED
done

# Set consistent indentation and spaces around operators
fprettify --indent 2 --whitespace 2 --strict-indent --stdout $CASED > $REGULARISED

# Remove temporary files
rm -v $TRIMMED $CASED

# Uncomment to inspect diff. of endresult with Visual Studio Code
# code --diff $IN $REGULARISED

# DANGER: Only uncomment this if you are sure there are no problems!
rm -v $IN
mv -v $REGULARISED $IN