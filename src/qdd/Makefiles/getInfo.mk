### getInfo.mk
# Collection and export of QDD source code status and compilation data

ifneq ($(NO_QDD_INFO), YES)
  GIT_HASH=$(shell git rev-parse HEAD)
  GIT_BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
  ifneq (,$(findstring dirty,$(shell git describe --dirty --always)))
    GIT_STATUS = Dirty
    FILES_CHANGED=$(shell git diff --name-status | cat -T | sed 's@\^I@,@g' | sed 's@/@>@g')
    NFILES_CHANGED=$(shell git diff --name-only | wc -l)
  else
    GIT_STATUS = Clean
    FILES_CHANGED=$(shell git diff --name-only | sed 's@.*/@@' | xargs | sed 's/ /,/g')
    NFILES_CHANGED=$(shell git diff --name-only | wc -l)
  endif
  COMPILE_TIME=$(shell date -u +'%Y/%m/%d %H:%M:%S UTC')
  # COMPILE_TIME=$(shell date -u +'%A %B %d %Y, %H:%M:%S UTC')
  UNAME := $(shell uname -s)
  VERSION_FLAGS= \
    -DGIT_HASH="\"$(GIT_HASH)\"" \
    -DGIT_BRANCH="\"$(GIT_BRANCH)\"" \
    -DGIT_STATUS="\"$(GIT_STATUS)\"" \
    -DFILES_CHANGED="\"$(FILES_CHANGED)\"" \
    -DNFILES_CHANGED="\"$(NFILES_CHANGED)\"" \
    -DCOMPILE_TIME="\"$(COMPILE_TIME)\""
  ifeq ($(UNAME), Linux)
    VERSION_FLAGS+= -D__linux
  endif
  export VERSION_FLAGS
endif