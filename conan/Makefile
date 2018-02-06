SUBDIRS = MPI

.PHONY: exports $(SUBDIRS)

exports: $(SUBDIRS)

$(SUBDIRS):
	(cd $@; conan create . cmsoft/stable)
