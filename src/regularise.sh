#!/usr/bin/env bash
#
# CALL: ./regularise.sh  (! execute in $QDD_ROOT/src/qdd/ !)
#
find . -name '*.f' | xargs perl -pi -e 's/(?<!#)endif/end if/ig'
find . -name '*.F' | xargs perl -pi -e 's/(?<!#)endif/end if/ig'
find . -name '*.f90' | xargs perl -pi -e 's/(?<!#)endif/end if/ig'
find . -name '*.F90' | xargs perl -pi -e 's/(?<!#)endif/end if/ig'

find . -name '*.f' | xargs perl -pi -e 's/(?<!#)enddo/end do/ig'
find . -name '*.F' | xargs perl -pi -e 's/(?<!#)enddo/end do/ig'
find . -name '*.f90' | xargs perl -pi -e 's/(?<!#)enddo/end do/ig'
find . -name '*.F90' | xargs perl -pi -e 's/(?<!#)enddo/end do/ig'

find . -name '*.F90' -exec ../regularise-ucF90.sh {} \; 2> regularise-ucF90.err
find . -name '*.f90' -exec ../regularise-lcF90.sh {} \; 2> regularise-lcF90.err