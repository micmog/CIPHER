SHELL := /bin/sh

.PHONY: install
install:
	@(cd src; make cipher; cp cipher.exe ../bin; cd -)

.PHONY: clean
clean:
	@(cd src; rm *.exe *.o; cd -)

