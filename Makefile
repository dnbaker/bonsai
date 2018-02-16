.PHONY=all tests clean

all:
	$(MAKE) -C bonsai

clean:
	$(MAKE) clean -C bonsai
